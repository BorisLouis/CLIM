%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2020-Data\09 - Sep\FilmBlinking\mov1';
file.ext  = '.spe';

info.runMethod  = 'load';
info.driftCorr = true;
info.ROI = true;

frame2Process = 1:3000;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.6;%correlation threshold (smaller is more correlation)==> 0.6 == 0.4 Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%

data = myMovie.loadFrames(frame2Process);

% vidFile = VideoWriter('rawMov.mp4','MPEG-4');
% vidFile.FrameRate = 100;
% open(vidFile);
% figure
% for i = 1:10:size(data,3)
%    imagesc(data(:,:,i));
%    colormap('hot')
%    caxis([0 max(data(:))]);
%    axis image
%    drawnow;
%    im = getframe;
%    writeVideo(vidFile,im);
%    clf
% end
% close(vidFile)


%% Get pixels correlation
[listCorrPx,inds] = myMovie.getPxCorrelation(data,corrInfo);

[corrMask,cleanedCorrMask] = myMovie.getCorrelationMask(data,corrInfo);

%compare the two clusters
[relNum1,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');




%% Correlation clustering

[corrMask,cleanedCorrMask] = myMovie.getCorrClustMask(data,corrInfo);

%compare the two clusters
[relNum3,relNum4] = compare2Cluster(corrMask,cleanedCorrMask,data,'corrClust');



%% ML Data Processing
MLOptions.clust2Test = [2,10];
MLOptions.GPU = true;
MLOptions.replicate = 10;
MLOptions.dist = false; %use dist between point as well as correlation

[MLCorrMask,cleanedMLCorrMask] = myMovie.getMLCorrelationMask(data,MLOptions);

[relNum5,relNum6] = compare2Cluster(MLCorrMask,cleanedMLCorrMask,data,'KMClust');


%% Hierarchical Clustering
[hierMask,corrMask,cleanedHierMask,cleanedCorrMask] = myMovie.getHierarchicalMask(data,MLOptions);

[relNum7,relNum8] = compare2Cluster(corrMask,cleanedCorrMask,data,'HierClust');


%% 
myMovie.showHierarchicalMask;

%% check mask

[meanTraces] = myMovie.checkMask(data,2);

%% clean ML Mask

[cleanMask, hierarchical] = corrAnalysis.cleanMLCorrMask(data,MLCorrMask,0.4);

%% clean mask

[cleanMask] = corrAnalysis.cleanCorrMask(data,myMovie.corrMask,0.9);

cleanMask = imfill(cleanMask,'holes');
se =strel('disk',1);
cleanMask = imclose(cleanMask,se);

figure
imagesc(cleanMask)
colormap('jet')
axis image




%% Correlation CLustering

listCorrPx = myMovie.listCorrPx;
inds = myMovie.indCorrPx;

distMap = corrAnalysis.getDistanceMapFromPxList(inds,data);
%get binary distance map (1 is correlated 0 in not correlated)
binaryDistMap = distMap<1-corrInfo.thresh;


%%






%% Plotting
myMovie.plotContour(data);

%% Plot traces
myMovie.plotTraces(data,3);

%% Extract intensity traces 
[traces] = myMovie.getIntensityTrace(data);
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
