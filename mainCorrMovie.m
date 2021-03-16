%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\03 - Mar\15 - Comparison ML Pseudo-Clustering';
file.ext  = '.spe';

info.runMethod  = 'load';
info.driftCorr = true;
info.ROI = true;

frame2Process = 1:6000;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.4;%correlation threshold (smaller is more correlation)==> 0.6 == 0.4 Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%

data = myMovie.loadFrames(frame2Process);

%% Get pixels correlation
[listCorrPx,inds] = myMovie.getPxCorrelation(data,corrInfo);
%

[corrMask,cleanedCorrMask] = myMovie.getCorrelationMask(data,corrInfo);
%

%compare the two clusters
[relNum1,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');

%% Plotting
myMovie.plotContour(data,'raw');

%% Plot traces
myMovie.plotClusterTraces(data,3);


%% Extract intensity traces 
[traces] = myMovie.getAllTraces(data);


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