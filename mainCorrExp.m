%%
clear;
close all;
clc;
%% User input

file.path = 'D:\Documents\Unif\PhD\Papers\13 - IntensityCorrletation\Simulations\Blinking\Blink Amp 125';
file.ext  = '.tif';

info.runMethod  = 'load';
info.frame2Process = 1:1000;
info.ROI = false;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)
info.corrInfo = corrInfo;
info.driftCorr = false;

%% Create an experiment object
myExperiment = Core.correlationExperiment(file,info);


%% Retrieve the movie in the input folder

myExperiment.retrieveMovies();

%% Get correlation mask

myExperiment.getCorrMask();

%% Save results
fileName = [myExperiment.path filesep 'Experiment.mat'];

save(fileName,'myExperiment');


