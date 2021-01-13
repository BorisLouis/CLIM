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
corrInfo.thresh = 0.6;%correlation threshold (smaller is more correlation)

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


%% 
[listCorrPx,inds] = myMovie.getPxCorrelation(data,corrInfo);

[corrMask] = myMovie.getCorrelationMask(data,corrInfo);


%% ML Data Processing
MLOptions.clust2Test = [2,10];
MLOptions.GPU = true;
MLOptions.replicate = 10;
MLOptions.dist = false; %use dist between point as well as correlation

[MLCorrMask] = myMovie.getMLCorrelationMask(data,MLOptions);


%%
[hierarchical] = myMovie.getHierarchicalMask(data,MLOptions);


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




%% Test Hierarchical clustering
[distanceMap]      = corrAnalysis.getDistanceMapFromPxList(inds,data);

testHClust = linkage(distanceMap);



%%
figure
dendrogram(testHClust)

%%






%% Plotting
myMovie.plotContour(data);

%% Plot traces
myMovie.plotTraces(data,3);

%% Extract intensity traces 
[traces] = myMovie.getIntensityTrace(data);
%%


