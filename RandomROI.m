clear
clc 
close all

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\Big Grain\Ambient\DriftCorr';
file.ext  = '.mat';

ROISize = 100;
nSample = 50;

%% load data

folder2Data = dir(file.path);
idx = contains({folder2Data.name},file.ext);

folder2Data = folder2Data(idx);

mov = load([folder2Data(1).folder filesep folder2Data(1).name]);
mov = mov.corrData;
%% get bounding box around data
figure 
subplot(1,2,1)
imagesc(mov(:,:,1))
colormap('gray')
BW = imbinarize(uint16(mov(:,:,1)),0.015);
subplot(1,2,2)
imagesc(BW);
se = strel('disk',10);
BW = imclose(BW,se);

box = regionprops(BW,'boundingBox');
box = round(box.BoundingBox);
%%
traceData = struct();
traceData(nSample).trace = [];
for i = 1: nSample
    
    % generate random ROI
    x = randi([box(1) box(1)+box(3)-ROISize],1);
    y = randi([box(2) box(2)+box(4)-ROISize],1);
    
    x = x:1:x+ROISize;
    y = y:1:y+ROISize;
    
    currTrace = squeeze(mean(mean(mov(y,x,:),1),2))';
    
    traceData(i).trace = currTrace;
    traceData(i).box   = [x(1),y(1), ROISize, ROISize];
end


%% save Data

filename = [file.path filesep 'traceData-ROI-' num2str(ROISize)];
save(filename,'traceData');


