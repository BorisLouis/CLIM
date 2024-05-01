%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2023 - data\12 - December 23\5 - Vacha Lab - Dry Objective - Marked Samples contact layers\on PCBM1_made on 1204\Mov11 - Marked 30%';
file.ext  = '';
corrInfo.thresh = 0.7;
info.runMethod = 'run';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.useThreshold = true;%false
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = false; %this is to use ROI for the whole analysis
%[x y  w h]
ROI = [];
%ROI = [5 71 230 120];
%for intensity extraction
method = 'SilWeigth'; %'Mean'
% For all Data:[5 71 230 120]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI
testROIRadius = 64; %radius of the ROI to find optimal threshold
frame2Process = 1000:1500; %number of frame to used for correlation analysis.
minCorr = 0.3;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.8;%maximum correlation to be tested, higher than 0.9 makes little sense

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



%% plot correlation map

%get px correlation
[corrRelation] = myMovie.getPxCorrelation(data2Use);

figure
imagesc(corrRelation.corrMap)
colormap('jet')
caxis([0 1])
axis image

%get the mask using the optimal threshold
[corrMask] = myMovie.getCorrelationMask(data2Use,corrInfo);


%% clean up mask (+ Silhouette map)

[cleanMask,silMap] = myMovie.cleanCorrMask(data2Use);

%% Old correlation metrics (histogram)
[bestClustEval1,bestRelNum] = myMovie.evalCluster(corrMask,data2Use);


color= 'colorcube';
[corrMaskIM] = myMovie.getImageFromMask(corrMask,color);

%% Extract intensity traces 
data = myMovie.loadFrames(1:myMovie.raw.maxFrame,ROI);

[traces] = myMovie.getAllTraces(data,correctedData,method);


%% Generate final Output

[corrOutput] = myMovie.generateResults;
