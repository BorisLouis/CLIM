%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\08 - August\Blinking Perovskite\2024-08-07 - nonCS';
file.ext  = '';

pxSize = 0.2;%in um
info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.useThreshold = true;%false
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = true; %this is to use ROI for the whole analysis
%[x y  w h]
%ROI = [];
ROI = [125 75 150 150];
%for intensity extraction
method = 'Mean'; %'Mean'
frame2Process = 1:2000; %number of frame to used for correlation analysis.

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
filename = dir(file.path);
idx = contains({filename.name},'.mat');

load([filename(idx).folder filesep filename(idx).name]);
%% Loading frames  
data1 = myMovie.loadFrames(1:myMovie.raw.movInfo.maxFrame,ROI);
%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1,backgroundThresh,deconvolve);

%% show the data
figure(1)
clf
subplot(1,2,1)
imagesc(squeeze(mean(correctedData,3)));
caxis([400 6000])
axis image
subplot(1,2,2)
imagesc(corrOutput.corrMap(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1)+ROI(3)));
caxis([0.5,1])
colormap('jet')

axis image


%% Look at pixel decorrelation
%CHANGE UNIT TO MICROMETER INSTEAD OF PIXEL
plotROI = [41 10 60 60];
data2Use = correctedData;

%radius = 30; %px

idx    =  [30,85];
%TCP: (X,Y)  High Corr 1) (61,35),(93,42), (84,22), (91 27)
%TCP: (X,Y)  low Corr 1) (85,35), (85,30)


%Big: [58,202] R2 [59,88] R3[38,110]
%Porous: R1[52,139]R2[30,93]  R3[80,87]
%Bigger grain [111 135];%[157,45];%[157 45];%[116,146];% [111 135]

[x,y] = meshgrid(1:size(data2Use(:,:,1)),1:size(data2Use(:,:,2)));
linIdx = sub2ind(size(data2Use(:,:,1)),y,x);
mainIdx = sub2ind(size(data2Use(:,:,1)),idx(1),idx(2));

data2Process = data2Use;

rd2Proc = reshape(data2Process,size(data2Process,1)*size(data2Process,2),size(data2Process,3));
rlinIdx = reshape(linIdx,size(data2Process,1)*size(data2Process,2),1);

correlationMap = corrcoef(double(rd2Proc)');
%correlationMap2 = corr(double(rd2Proc)','Type','Kendall');

corrData = correlationMap(:,rlinIdx==mainIdx);

rCorrData = reshape(corrData,size(data2Process,1),size(data2Process,2));

tmpIm = zeros(size(data2Use,1),size(data2Use,2));
tmpIm(linIdx) = rCorrData;

% plots
figure
imagesc(1:size(rCorrData,1)*pxSize,1:size(rCorrData,2)*pxSize,rCorrData);
colormap('jet')
caxis([0 1])
xlabel('Position in \mum')
ylabel('Position in \mum')
colorbar
axis image
xlim([plotROI(1)*pxSize plotROI(1)*pxSize+plotROI(3)*pxSize])
ylim([plotROI(2)*pxSize plotROI(2)*pxSize+plotROI(4)*pxSize])

figure
imagesc(1:size(rCorrData,1)*pxSize,1:size(rCorrData,2)*pxSize,rCorrData);
colormap('jet')
caxis([-1 1])
xlabel('Position in \mum')
ylabel('Position in \mum')
colorbar
xlim([plotROI(1)*pxSize plotROI(1)*pxSize+plotROI(3)*pxSize])
ylim([plotROI(2)*pxSize plotROI(2)*pxSize+plotROI(4)*pxSize])


figure
surf(rCorrData)

figure 
imagesc(imfuse(squeeze(mean(data2Process,3)),rCorrData))
xlim([plotROI(1) plotROI(1)+plotROI(3)])
ylim([plotROI(2) plotROI(2)+plotROI(4)])


figure 
imagesc(imfuse(corrOutput.corrMap(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1)+ROI(3)),rCorrData))
xlim([plotROI(1) plotROI(1)+plotROI(3)])
ylim([plotROI(2) plotROI(2)+plotROI(4)])


%% Posiion without micrometer
figure
imagesc(rCorrData);
colormap('jet')
caxis([0 1])
% xlabel('Position in \mum')
% ylabel('Position in \mum')
colorbar
axis image



%% Plot selected intensity traces
idx1 = [14 42];%[Y X]
idx2 = [21 10];
idx3 = [80 46];
idx4 = [51 51];%[Y X]
idx5 = [85 25];%[Y X]
idx6 = [35 56];

figure
plot(squeeze(data2Process(idx1(1),idx1(2),:)))
hold on
plot(squeeze(data2Process(idx2(1),idx2(2),:)))
plot(squeeze(data2Process(idx3(1),idx3(2),:)))
plot(squeeze(data2Process(idx4(1),idx4(2),:)))
plot(squeeze(data2Process(idx5(1),idx5(2),:)))
plot(squeeze(data2Process(idx6(1),idx6(2),:)))

legend({'R1','R2','R3','Ref','R4','R5'})


ylim ([1000 4000])
axis square
box on
