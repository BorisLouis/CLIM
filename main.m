%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2020-Data\10 - Oct\Sudipta\nonBlinking Data\movOD3';
file.ext  = '.spe';

info.runMethod  = 'load';
frame2Process = 1:1000;
corrInfo.r = 2; %radius for checking neighbor
corrInfo.thresh = 0.3;%correlation threshold (smaller is more correlation)
driftCorr = true;
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

%% Extract intensity traces 
idx = 1;
[row,col] = find(corrMask==idx);
trace = zeros(length(row),size(data,3));
for i = 1:length(row)
   
    trace(i,:) = data(row(i),col(i),:);
    
end

