%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Keep in mind that the deconvolution will affect the background, the
% best would be to make an ROI that removes the background to conteract
% this
%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\12 - December\12 - Code improve attempt\BigGrain';
file.ext  = '.spe';

info.runMethod  = 'run';
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



%% Get stat based correlation

distMap = corrAnalysis.getDistanceMapFromPxList(corrRelation.indPx,correctedData);

%indPx can have holes

corrMask = corrAnalysis.corrClusteringV4(correctedData,corrRelation,distMap);















