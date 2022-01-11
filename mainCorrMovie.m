%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\22-23 - All Measurement\22.11.21\Small grain batch1_2\Place1\Air';
file.ext  = '.spe';

info.runMethod  = 'run';%
info.driftCorr = true;
info.ROI = false;%this is to use ROI for the whole analysis
%      [x y  w h]
ROI = [64,64,96,96]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI

frame2Process = 1:6000;

minCorr = 0.4;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.9;%maximum correlation to be tested, higher than 0.9 makes little sense

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%% Loading frames  
data1 = myMovie.loadFrames(frame2Process,ROI);


%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1);

%% Scanning threshold
if info.ROI ==false
    ROICorrData = correctedData(ROI(2):ROI(2)+ROI(4)-1,ROI(1):ROI(1)+ROI(3)-1,:);
else
    ROICorrData = correctedData;
end

thresh = minCorr:stepCorr:maxCorr;

[allCorrMask,threshold2Use] = myMovie.findOptimalThreshold(ROICorrData,thresh);

%% Calculate final corrMask

corrInfo.thresh = threshold2Use;
%get px correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData,corrInfo);
%get the mask using the optimal threshold
[corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);


%% clean up mask (+ Silhouette map)
% THIS STEP IS VERY LONG
tic
[cleanMask,silMap] = myMovie.cleanCorrMask(data1);
toc

%% Old correlation metrics (histogram)
%compare the two clusters
%[~,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');
[bestClustEval1,bestRelNum] = corrAnalysis.evalClusters(corrMask,correctedData);

relData{1} = bestRelNum;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);

%% Plotting
myMovie.plotContour(data1);%raw or clean depending on which we want to use

%% Plot Contour

myMovie.plotContour(corrRelation.corrMap,corrMask);

%% Plot traces
myMovie.plotClusterTraces(data1,4);


%% Extract intensity traces 
data = myMovie.loadFrames(1:30000);
[traces] = myMovie.getAllTraces(data);


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

%% get image from corrmask
color= 'colorcube';
[corrMaskIM] = myMovie.getImageFromMask(corrMask,color);

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