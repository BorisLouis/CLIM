%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\08 - August\SOFI\large grain';
file.ext  = '.tif';

info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.thresholdMode = 'fixed'; %'auto', 'fixed' or 'None'
threshold = 0.4; %only used if info.thresholdMode is "fixed"
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
frame2Process = 1:1000; %number of frame to used for correlation analysis.
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


%% Plot and save cluster with numbers image
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
