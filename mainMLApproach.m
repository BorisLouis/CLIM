%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Keep in mind that the deconvolution will affect the background, the
% best would be to make an ROI that removes the background to conteract
% this
%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\30 - Big algorithm evaluation\Small Grain';
file.ext  = '.spe';

info.runMethod  = 'load';
info.driftCorr = true;
info.ROI = true;

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.thresh = 0.6;%correlation threshold Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%
ROI = [26,100,64,64];
data1 = myMovie.loadFrames(frame2Process,ROI);


%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1);
%% Get pixels correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData,corrInfo);


%% Approach 1 - DBScan
testData = reshape(correctedData,[size(correctedData,1)*size(correctedData,1), size(correctedData,3)]);

id = 1:size(testData,1);

% idx = dbscan(double(testData),0.1,3,'Distance','Correlation');
% 
% corrMask = zeros(size(correctedData,1),size(correctedData,2));
% 
% corrMask(id) = idx;
% 
% figure
% imagesc(corrMask)
% 
% %% Approach 1 - Evaluation
% [clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,correctedData);
% 
% relData{1} = relNum1;
% label{1}   = ['Method' '-pseudoClust'];
% corrAnalysis.compareClusters(relData,label);
% 
% myMovie.plotContour(data1,corrMask);
%% Approach 2 - KMeans

clustCenters = imregionalmax(corrRelation.corrMap);

nClusters = sum(clustCenters(:));

distanceMap = 1-corrcoef(double(testData'));


idx = kmeans(double(testData),nClusters,'Distance','correlation');

corrMask = zeros(size(correctedData,1),size(correctedData,2));

corrMask(id) = idx;

figure
imagesc(corrMask)
colormap('colorcube')
axis image

%% Approach 2 - Evaluation
[clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,correctedData);

relData{1} = relNum1;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);

myMovie.plotContour(corrRelation.corrMap,corrMask);