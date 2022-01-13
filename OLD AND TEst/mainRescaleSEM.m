clear;
close all;
clc;
%% User input 
pxSizeRef = 200; %in nm;
refDim    = [248, 301];
zoomMultiple = 3; %for zoom in SEM we make the image bigger to be matched with a zoomed in fluorescence image
folder2Analyze = 'F:\Boris - Lund\2021\05 - Mai\28.5.21\Big grain with pores\SEM';

%% Load Data
folder2Process = dir(folder2Analyze);

imageList =  folder2Process(contains({folder2Process.name},'.tif'));
infoList  =  folder2Process(contains({folder2Process.name},'.txt'));
assert(numel(imageList)==numel(infoList),'number of tif files and txt file are not equal, something is wrong');


%% loop through images, rescale and save new image

for i=1:numel(imageList)
    [~,fName,~] = fileparts(imageList(i).name);
    currentIm = [imageList(i).folder filesep imageList(i).name];
    currentIm = imread(currentIm);
    
    currentInfo = [infoList(i).folder filesep infoList(i).name];
    currentInfo = readtable(currentInfo);
    idx = strcmp(currentInfo.InstructName,'PixelSize');
    pxSize = currentInfo.SU8000(idx);
    
    scalingFactor = pxSize*max(size(currentIm))/(pxSizeRef*max(refDim));
    if scalingFactor<0.4
        %when the scaling factor get really small (zoom in SEM) we multiply
        %the size by 3 so we can match
        scalingFactor = scalingFactor*zoomMultiple;
    end
    imRescaled = imresize(currentIm,scalingFactor);
    
    imwrite(imRescaled,[imageList(i).folder filesep fName '-rescaled.png']);
    
    
end