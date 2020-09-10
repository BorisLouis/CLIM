%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2020-Data\09 - Sep\FilmBlinking\testData';
file.ext  = '.spe';

info.runMethod  = 'load';
frame2Process = 1:2000;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

data = myMovie.loadFrames(frame2Process);

%% Data Processing
[corrMask] = myMovie.getCorrelationMask(data,corrInfo);