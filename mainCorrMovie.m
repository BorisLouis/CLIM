%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2023 - data\11 - November\22 - Second day Vacha\perovskite with PEDOT2\Power 2 33ms';
file.ext  = '';

info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0; %0.1 is default, 0 is for no background removal
info.useThreshold = false;%false
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = false; %this is to use ROI for the whole analysis
%      [x y  w h]
ROI = [];
%ROI = [5 71 230 120];
%for intensity extraction
method = 'SilWeigth'; %'Mean'
% For all Data:[5 71 230 120]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI
testROIRadius = 64; %radius of the ROI to find optimal threshold
frame2Process = 500:1500; %number of frame to used for correlation analysis.
minCorr = 0.4;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.8;%maximum correlation to be tested, higher than 0.9 makes little sense

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
%% Loading frames  
data1 = myMovie.loadFrames(frame2Process,ROI);

%% 
%myMovie.saveMovie(data1,50);
%% Deconvolution
if deconvolve
    [correctedData] =  myMovie.deconvolve(data1,backgroundThresh);
else
    correctedData = data1;
end
%% Scanning threshold
if info.useThreshold
    center = [round(size(correctedData,1)/2), round(size(correctedData,2)/2)];
    if info.ROI ==false
        myMovie.info.ROIUsed = [];    
    else
        ROICorrData = correctedData;
        myMovie.info.ROIUsed = ROI;
    end
    try

        ROICorrData = correctedData(center(1)-testROIRadius:center(1)+testROIRadius-1,center(2)-testROIRadius:center(2)+testROIRadius-1,:);
    catch except
        if strcmp(except.identifier, 'MATLAB:badsubscript')
            ROICorrData = correctedData;
        end
    end

    thresh = minCorr:stepCorr:maxCorr;

    [allCorrMask,threshold2Use] = myMovie.findOptimalThreshold(ROICorrData,thresh);
end
%% Calculate final corrMask
if info.useThreshold
    corrInfo.thresh = threshold2Use;
else
    corrInfo.thresh = minCorr;
end
%get px correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData);
%get the mask using the optimal threshold
[corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);


%% clean up mask (+ Silhouette map)

[cleanMask,silMap] = myMovie.cleanCorrMask(correctedData);

%% Old correlation metrics (histogram)
[bestClustEval1,bestRelNum] = myMovie.evalCluster(corrMask,correctedData);


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
data = myMovie.loadFrames(1:1000,ROI);

[traces] = myMovie.getAllTraces(correctedData,method);





%% Generate final Output

[corrOutput] = myMovie.generateResults;


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