%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2023 - data\11 - November\29 - Vacha\perovskite on glass\Long WL low pow good';
file.ext  = '';

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
data1 = myMovie.loadFrames(frame2Process,ROI);

%% save data 
%myMovie.saveMovie(data1,50);
dataStorage.nBTiff('driftCorrected',data1,16);

%% Deconvolution
if deconvolve
    [correctedData] =  myMovie.deconvolve(data1,backgroundThresh);
else
    correctedData = data1;
end



%% plot correlation map

%get px correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData);

figure
imagesc(corrRelation.corrMap)
colormap('jet')
caxis([0 1])
axis image