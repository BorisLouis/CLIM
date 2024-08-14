%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\08 - August\Blinking Perovskite\Spot1_dryObj glass side';
file.ext  = '';

info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.thresholdMode = 'auto'; %'auto', 'fixed' or 'None'
threshold = 0.5; %only used if info.thresholdMode is "fixed"
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = false; %this is to use ROI for the whole analysis
%[x y  w h]
ROI = [];
%ROI = [5 71 230 120];
%for intensity extraction
method = 'Mean'; %'Mean'
% For all Data:[5 71 230 120]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI
testROIRadius = 64; %radius of the ROI to find optimal threshold
frame2Process = 1:2000; %number of frame to used for correlation analysis.
minCorr = 0.4;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.7;%maximum correlation to be tested, higher than 0.9 makes little sense

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
%% Loading frames  
data1 = myMovie.loadFrames(1:myMovie.raw.movInfo.maxFrame,ROI);

%% save data 
%myMovie.saveMovie(data1,50);
%dataStorage.nBTiff('driftCorrected.tif',data1,16);

%% Deconvolution

[correctedData] =  myMovie.deconvolve(data1,backgroundThresh,deconvolve);


data2Use = correctedData(:,:,frame2Process);
%% Scanning threshold
if strcmpi(info.thresholdMode,'auto')
    center = [round(size(data2Use,1)/2), round(size(data2Use,2)/2)];
    if info.ROI ==false
        myMovie.info.ROIUsed = [];    
    else
        ROICorrData = data2Use;
        myMovie.info.ROIUsed = ROI;
    end
    try

        ROICorrData = data2Use(center(1)-testROIRadius:center(1)+testROIRadius-1,center(2)-testROIRadius:center(2)+testROIRadius-1,:);
    catch except
        if strcmp(except.identifier, 'MATLAB:badsubscript')
            ROICorrData = data2Use;
        end
    end

    thresh = minCorr:stepCorr:maxCorr;

    [allCorrMask,threshold2Use] = myMovie.findOptimalThreshold(ROICorrData,thresh);
end
%% Calculate final corrMask
switch lower(info.thresholdMode)
    case 'auto'
        corrInfo.thresh = threshold2Use;
    case 'fixed'
        corrInfo.thresh = threshold;
    case 'none'
    corrInfo.thresh = minCorr;
end
%get px correlation
[corrRelation] = myMovie.getPxCorrelation(data2Use);
%get the mask using the optimal threshold
[corrMask] = myMovie.getCorrelationMask(data2Use,corrInfo);

%% clean up mask (+ Silhouette map)

[cleanMask,silMap] = myMovie.cleanCorrMask(data2Use);

%% Old correlation metrics (histogram)
[bestClustEval1,bestRelNum] = myMovie.evalCluster(corrMask,data2Use);


%% Plotting
myMovie.plotContour(data1(:,:,1),corrMask);%raw or clean depending on which we want to use

%% Plot Contour

myMovie.plotContour(corrRelation.corrMap,corrMask);

%% Plot traces
myMovie.plotClusterTraces(data1,4);

%% get image from corrmask with numbered
color= 'colorcube';
[corrMaskIM] = myMovie.getImageFromMask(corrMask,color);

%% Extract intensity traces 
data = myMovie.loadFrames(1:myMovie.raw.maxFrame,ROI);

[traces] = myMovie.getAllTraces(data,correctedData,method);





%% Generate final Output

[corrOutput] = myMovie.generateResults;




%% 
RGBIM = label2rgb(corrOutput.corrMask,color,'k','shuffle');

f1 = figure;
imagesc(RGBIM)
axis image

hold on

centers = regionprops(corrOutput.corrMask,'Centroid');

for i = 1:length(unique(corrOutput.corrMask(corrOutput.corrMask>0)))
     x = centers(i).Centroid(1);
     y = centers(i).Centroid(2);
     
     text(x,y,['.',num2str(i)],'FontSize',8,'Color','white');
end

[folder,file,ext] = fileparts(corrOutput.path);

fileName = [folder filesep 'clusterID.fig'];
saveas(f1,fileName)

fileName = [folder filesep 'ClusterID.svg'];
saveas(f1,fileName,'svg')


%% Get interCluster traces
% allTraces = zeros(size(traces,2),length(traces(1,1).trace));
% for i = 1:size(traces,2)
%     allTraces(i,:) = traces(1,i).trace;
%     
% end
% allTraces = allTraces';
% corrCoeff = corrcoef(allTraces);
% 
% figure
% imagesc(corrCoeff)
% colormap('jet')
% 
% U = triu(corrCoeff);
% U = nonzeros(U);
% tmpData = U(U<1);
% figure
% histogram(tmpData);
% 
% 
% % Let us try to get anti-correlated Crystal
% [val,idx] = mink(corrCoeff(:),20);
% %remove duplicate
% val = val(1:2:end);
% 
% 
% for i=1:length(val)
%     currVal = val(i);
%     idx = find(corrCoeff==currVal);
%     [row,col] = ind2sub(size(corrCoeff),idx(2));
% 
%     figure
%     hold on
%     plot(allTraces(:,row))
%     plot(allTraces(:,col))
%     test = corrcoef(allTraces(:,row),allTraces(:,col));
%     disp(test);
% end



%% Calculate localization on cluster to get the position of the defect
%[allLoc] =myMovie.getAllClusterLocalization(data);

idx =112;
[Localization] = myMovie.getClusterLocalization(data,idx);

pos = sqrt(Localization.x.^2 + Localization.y.^2);
figure
subplot(1,3,1)
scatter(pos,Localization.int,10,'filled')
axis square
box on
title('Intensity vs Position')

subplot(1,3,2)
scatter(pos,Localization.angle,10,'filled')
axis square
box on
title('Position vs angle')

subplot(1,3,3)
scatter(Localization.x,Localization.y,10,'filled')
axis square
box on
title('Intensity vs Angle')

%% test plot
myMovie.plotContour(data);%raw or clean depending on which we want to use
hold on
scatter(Localization.x,Localization.y,5,'filled')

%%

function [relNum1,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,method)
    %Evaluate clusters individually
    [clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,data);
    [clustEval2,relNum2] = corrAnalysis.evalClusters(cleanedCorrMask,data);


    %compare cluster
    relData{1} = relNum1;
    relData{2} = relNum2;
    label{1}   = ['Method' method];
    label{2}   = ['Method',method,'- cleaned'];

    corrAnalysis.compareClusters(relData,label);


end