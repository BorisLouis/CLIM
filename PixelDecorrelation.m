%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\22-23 - All Measurement\22.11.21\Big grain batch1_2\place1\N2';
file.ext  = '';

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
ROI = [5 71 230 120];
%for intensity extraction
method = 'Mean'; %'Mean'
% For all Data:[5 71 230 120]; %this will be use for scanning threshold and/or the whole analysis based on info.ROI
testROIRadius = 64; %radius of the ROI to find optimal threshold
frame2Process = 1:2000; %number of frame to used for correlation analysis.
minCorr = 0.4;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.7;%maximum correlation to be tested, higher than 0.9 makes little sense

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
%% Loading frames  
data1 = myMovie.loadFrames(1:myMovie.raw.movInfo.maxFrame,ROI);

%% save data 
%myMovie.saveMovie(data1,50);
%dataStorage.nBTiff('driftCorrected.tif',data1,16);

%% Deconvolution

[correctedData] =  myMovie.deconvolve(data1,backgroundThresh,deconvolve);

%%
%CHANGE UNIT TO MICROMETER INSTEAD OF PIXEL
data2Use = correctedData;
pxSize = 200;

radius = 27; %px

idx    =  [58,202];%[59,88];
%Big: [58,202] R2 [59,88] R3[38,110]
%Porous: R1[52,139]R2[30,93]  R3[80,87]
%Bigger grain [111 135];%[157,45];%[157 45];%[116,146];% [111 135]

[x,y] = meshgrid(idx(2)-radius:idx(2)+radius,idx(1)-radius:idx(1)+radius);
linIdx = sub2ind(size(data2Use(:,:,1)),y,x);
mainIdx = sub2ind(size(data2Use(:,:,1)),idx(1),idx(2));

data2Process = data2Use(idx(1)-radius:idx(1)+radius,idx(2)-radius:idx(2)+radius,:);

rd2Proc = reshape(data2Process,size(data2Process,1)*size(data2Process,2),size(data2Process,3));
rlinIdx = reshape(linIdx,size(data2Process,1)*size(data2Process,2),1);

correlationMap = corrcoef(double(rd2Proc)');
%correlationMap2 = corr(double(rd2Proc)','Type','Kendall');

corrData = correlationMap(:,rlinIdx==mainIdx);

rCorrData = reshape(corrData,size(data2Process,1),size(data2Process,2));

tmpIm = zeros(size(data2Use,1),size(data2Use,2));
tmpIm(linIdx) = rCorrData;


figure
imagesc(1:size(rCorrData,1)*pxSize/1000,1:size(rCorrData,2)*pxSize/1000,rCorrData);
colormap('jet')
caxis([0 1])
xlabel('Position in \mum')
ylabel('Position in \mum')
colorbar
axis image


figure
surf(rCorrData)

figure 
imagesc(imfuse(data2Use(:,:,1),tmpIm))


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
