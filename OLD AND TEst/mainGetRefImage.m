clear;
close all;
clc;
%%
zoomMultiple = 3; %for zoom in SEM we make the image bigger to be matched 
folder2Analyze = 'F:\Boris - Lund\2021\05 - Mai\31.5.21\BigGrain - Aymens Sample';

folder2Process = dir(folder2Analyze);

imageList =  folder2Process(contains({folder2Process.name},'.spe'));
cmap = colormap(hot);
close;
for i=1:numel(imageList)
    currentIm = [imageList(i).folder filesep imageList(i).name];
    currentIm = Load.Movie.spe.getFrame(currentIm);
    [~,filename,~] = fileparts(imageList(i).name);
    meanIm = mean(currentIm,3);
    maxIm  = max(currentIm,[],3);
    
    scaledMeanIm = imresize(meanIm,zoomMultiple);
    scaledMaxIm  = imresize(maxIm,zoomMultiple);
    
    %normalization
    normMeanIm = meanIm -min(meanIm(:));
    normMeanIm = normMeanIm./max(normMeanIm(:))*255;
    
    normMaxIm = maxIm -min(maxIm(:));
    normMaxIm = normMaxIm./max(normMaxIm(:))*255;
    
    scaledNormMeanIm = scaledMeanIm -min(scaledMeanIm(:));
    scaledNormMeanIm = scaledNormMeanIm./max(scaledNormMeanIm(:))*255;
    
    scaledNormMaxIm = scaledMaxIm -min(scaledMaxIm(:));
    scaledNormMaxIm = scaledNormMaxIm./max(scaledNormMaxIm(:))*255;
    
    tmpMean = ind2rgb(uint8(normMeanIm),cmap);
    imwrite(tmpMean,cmap, [imageList(i).folder filesep filename '-mean.png']);
    
    tmpMax = ind2rgb(uint8(normMaxIm),cmap);
    imwrite(tmpMax, cmap, [imageList(i).folder filesep filename '-max.png']);
    
    tmpMean = ind2rgb(uint8(scaledNormMeanIm),cmap);
    imwrite(tmpMean,cmap, [imageList(i).folder filesep filename '-meanZoom.png']);
    
    tmpMax = ind2rgb(uint8(scaledNormMaxIm),cmap);
    imwrite(tmpMax, cmap, [imageList(i).folder filesep filename '-maxZoom.png']);
    
    
end
disp('=======>DONE<=======');



%%