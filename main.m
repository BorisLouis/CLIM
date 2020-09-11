%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2020-Data\09 - Sep\FilmBlinking\mov1';
file.ext  = '.spe';

info.runMethod  = 'load';
frame2Process = 1:100;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)
driftCorr = true;
%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

data = myMovie.loadFrames(frame2Process,driftCorr);

vidFile = VideoWriter('rawMov.mp4','MPEG-4');
vidFile.FrameRate = 100;
open(vidFile);
figure
for i = 1:size(data,3)
   imagesc(data(:,:,i));
   colormap('hot')
   caxis([0 max(data(:))]);
   axis image
   drawnow;
   im = getframe;
   writeVideo(vidFile,im);
   clf
end
close(vidFile)


%% Data Processing
[corrMask] = myMovie.getCorrelationMask(data,corrInfo);

%% Plotting
maxIm = max(data,[],3);

figure
hold on
imagesc(maxIm)
axis image
colormap('hot')
for i = 1:max(corrMask(:))
    corrMaskCopy = corrMask;
    
    corrMaskCopy(corrMask~=i) = 0;
    
    contour = bwboundaries(corrMaskCopy);
    
    plot(contour{1}(:,2),contour{1}(:,1),'w','LineWidth',2)
    
    
end
axis ij


