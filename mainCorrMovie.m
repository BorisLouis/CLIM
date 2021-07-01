%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\06 - June\BlinkingEx\4Cluster';
file.ext  = '.tif';

info.runMethod  = 'load';
info.driftCorr = true;
info.ROI = false;

frame2Process = 1:300;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)==> 0.6 == 0.4 Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%

data = myMovie.loadFrames(frame2Process);

%% Get pixels correlation
[corrR] = myMovie.getPxCorrelation(data,corrInfo);
%

[corrMask,cleanedCorrMask] = myMovie.getCorrelationMask(data,corrInfo);
%

%compare the two clusters
[relNum1,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');

%% Plotting
myMovie.plotContour(data,'raw');%raw or clean depending on which we want to use


%% Plot traces
myMovie.plotClusterTraces(data,4);


%% Extract intensity traces 
[traces] = myMovie.getAllTraces(data);


%% get image from corrmask
color= 'jet';
[corrMaskIM] = myMovie.getImageFromMask(corrMask,color);

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