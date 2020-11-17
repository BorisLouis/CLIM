%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2020-Data\10 - Oct\Sudipta\nonBlinking Data\movOD3';
file.ext  = '.spe';

info.runMethod  = 'load';
info.driftCorr = true;

frame2Process = 1:1000;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)

%% Loading data
myMovie = Core.CorrClusterMovie(file,info,driftCorr);

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


%% Data Processing
[corrMask] = myMovie.getCorrelationMask(data,corrInfo);

%% Plotting
myMovie.plotContour(data);

%% Plot traces
myMovie.plotTraces(data,3);

%% Extract intensity traces 
[traces] = myMovie.getIntensityTrace(data);
%%

