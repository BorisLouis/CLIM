%%
clear;
close all;
clc;
%% User input

file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\10 - October\Resolution limit\FWHM3-2';
file.ext  = '.tif';
info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = false; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.thresholdMode = 'None'; %'auto', 'fixed' or 'None'
info.threshold = 0.4; %only used if info.thresholdMode is "fixed"
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
info.frame2Process = 1:1000; %number of frame to used for correlation analysis.
info.minCorr = 0.4;%Minimum correlation we want to have
info.stepCorr = 0.05; %Correlation difference between different tested threshold
info.maxCorr = 0.7;%maximum correlation to be tested, higher than 0.9 makes litt

%% Create an experiment object
myExperiment = Core.correlationExperiment(file,info);


%% Retrieve the movie in the input folder

myExperiment.retrieveMovies();

%% Get correlation mask

myExperiment.getCorrMask();

%% Save results
fileName = [myExperiment.path filesep 'Experiment.mat'];

save(fileName,'myExperiment');




%%

e=figure(1);
ax1 = gca;



for i =1:size(myExperiment.corrMasks,2)
    
    subplot(2,3,i)
    imagesc(myExperiment.corrMasks(i).mask(1:30,1:30))
    axis image
    colormap('hot')
end
   


f=figure(2);
ax2 = gca;
for i =1:size(myExperiment.corrMasks,2)
    
    subplot(2,3,i)
    imagesc(myExperiment.corrMasks(i).data(1:30,1:30))
    axis image
    colormap('hot')
    drawnow;
   
end
   


