%%
clear;
close all;
clc;
%% User input
file.path = 'D:\Documents\Unif\PostDoc\2023 - data\01 - January\17 - Correlation analysis\File 2 (No biais)';
file.ext  = '';

refFrame = 2500; %Here is the frame to which the correlation is calculated

info.runMethod = 'run';%load % load will try to load existing data from previous run
info.driftCorr = false; % WE DONT CORRECT AT THE MOMENT
deconvolve = true; %to deconvolve the correlated signal
info.useThreshold = false;%false
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = false; %this is to use ROI for the whole analysis
%      [x y  w h]
ROI = [];
% For all Data:[5 71 230 120]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI
testROIRadius = 32; %radius of the ROI to find optimal threshold
frame2Process = 1:6000; %number of frame to used for correlation analysis.


%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
%% Loading frames  
data1 = myMovie.loadFrames(frame2Process,ROI);


%%
correlationProfile = zeros(size(data1,3),1);
for i = 1:size(data1,3)
   currentIm = data1(:,:,i);
   
   currentCorr = corr2(data1(:,:,refFrame),currentIm);
   
   correlationProfile(i) = currentCorr;
   
end

figure
plot(correlationProfile)

