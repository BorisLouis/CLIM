%% Simulate correlated data

%% Simulation input
sizeIm = 64;
nFrames = 250;

nParticles = 4;
corrThreshold = 0.7;%smaller is more selective here (0 is perfect correlation) 0.7 is 0.3 pearson coefficient
model.name = 'gaussian';
model.sigma_x = 3;
model.sigma_y = 3;
r = 1; %radius for checking correlation
signThreshold = 0.8; %correlation of pixel to other cluster needs to be smaller than 0.8*correlation of current cluster

%% Simulations
data = zeros(sizeIm,sizeIm,nFrames);
coord = [28,28;36,36;28,36;36,28];
%coord = [32,32;96,96;96,32;32,96];
[X,Y] = meshgrid(1:sizeIm,1:sizeIm);
%coord = zeros(nParticles,2);
for i = 1:nParticles
    
    x0 = coord(i,1);
    y0 = coord(i,2);
    c  = 0;
    BaseInt = 5e4;
    secondInt = 2*BaseInt;
    int = BaseInt;
    for j = 1 : nFrames
        
        num = rand(1);
        if num > 0.7
            if int==BaseInt
                int=secondInt;
            else
                int = BaseInt;
            end
        else
        end
        
        PSF = Sim.getPSF(X,Y,x0,y0,model);
        sPSF = Sim.samplePSF(PSF,int,false);
        
        data(:,:,j) = data(:,:,j) + double(sPSF); 


    end
    coord(i,:) = [x0,y0];   
end

noise = randn(size(data));

finalData = data + ones(size(data))*100 + noise*20;
%% savedata as Movie
%% save a movie as example
% vidFile = VideoWriter('rawMov.mp4','MPEG-4');
% vidFile.FrameRate = 10;
% open(vidFile);
% figure
% for i = 1:size(finalData,3)
%    imagesc(finalData(:,:,i));
%    colormap('hot')
%    caxis([0 max(finalData(:))]);
%    axis image
%    drawnow;
%    im = getframe;
%    writeVideo(vidFile,im);
%    clf
% end
% close(vidFile)


%% get pixel correlation method 1
[listCorrPx]  = corrAnalysis.getCorrRelation(finalData,r,corrThreshold);

%% get pixel correlation method 2
[listCorrPx,corrSum]  = corrAnalysis.getCorrRelation2(finalData,r,corrThreshold);

%% get pixel correlation method 3
[listCorrPx,corrSum3]  = corrAnalysis.getCorrRelation3(finalData,r,corrThreshold);

%%
%listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
inds    = (1:length(listCorrPx))';

%#4 Clean data by keeping only pixel that have correlation
%relation
idx2Delete = cellfun(@isempty,listCorrPx);
listCorrPx(idx2Delete) =[];
inds(idx2Delete) = [];

%calculate distancemap
[n,p] = ind2sub(size(finalData),inds);
pxIntList = zeros(length(n),size(finalData,3));
for i =1:length(n)

    pxIntList(i,:) = finalData(n(i),p(i),:);

end
distanceMap = 1-corrcoef(pxIntList');

%perform pseudo-clustering
dim = size(finalData);
[corrMask] = corrAnalysis.corrClustering(listCorrPx,inds,distanceMap,dim,corrThreshold);


%% clean up mask
%TODO: Code function to clean up the mask by testing significance of
%correlation found(is it significantly higher than the correlation to other
%cluster?)

[cleanCorrMask] = corrAnalysis.cleanCorrMask(finalData,corrMask,signThreshold);


%% SR test
BW = imbinarize(corrSum);
BW = bwlabel(BW);

%% ML test for corrmask
%get linear indices for all pixel in the image
inds    = (1:length(listCorrPx))';
%#4 Clean data by keeping only pixel that have correlation
%relation
idx2Delete = cellfun(@isempty,listCorrPx);
listCorrPx(idx2Delete) =[];
inds(idx2Delete) = [];

[distanceMap]      = corrAnalysis.getDistanceMapFromPxList(inds,finalData);

[clust,clustEval]  = corrAnalysis.clusterCorrelatedPixel(distanceMap,[1,10]);

clust2Use = clustEval.OptimalK;

[MLCorrMask]       = corrAnalysis.getMaskFromMLCluster(clust,inds,clust2Use,size(finalData(:,:,1)));


figure
imagesc(MLCorrMask)
axis image
colormap('jet')

[cleanCorrMask] = corrAnalysis.cleanCorrMask(finalData,MLCorrMask,corrThreshold);

figure
imagesc(cleanCorrMask)
axis image
colormap('jet')
%% Try Hierarchical clustering
inds    = (1:length(listCorrPx))';
%#4 Clean data by keeping only pixel that have correlation
%relation
idx2Delete = cellfun(@isempty,listCorrPx);
inds(idx2Delete) = [];

[distanceMap]      = corrAnalysis.getDistanceMapFromPxList(inds,data);
y = squareform(distanceMap);
testHClust = linkage(distanceMap);

figure
dendrogram(testHClust)



%% Try ML from Scratch
MLData = reshape(finalData,[size(finalData,1)*size(finalData,2),size(finalData,3)]);

%normalize data
minData = min(MLData,[],2);
MLData = MLData-minData;
MLData = MLData./max(MLData,[],2);

nClust = 10;
clust = zeros(size(MLData,1),nClust);
dist = cell(nClust,1);
for i=1:nClust
    [clust(:,i),~,~,dist{i}] = kmeans(MLData,i,'Distance','correlation','emptyaction','drop',...
        'replicate',10);
end
va = evalclusters(MLData,clust,'CalinskiHarabasz');

%clean based on distance
distToClean = dist{4,1};
distThresh = distToClean<corrThreshold;
idx2delete = sum(distThresh,2);
idx2delete = idx2delete<1;

clust(idx2delete,4) = 0;


MLCorrMask = zeros(size(finalData,1),size(finalData,2));
for i = 1:length(MLData)
    MLCorrMask(i) = clust(i,4);
end

figure
subplot(1,2,1)
imagesc(MLCorrMask)
title('ML output')
axis image
colormap('jet')

[cleanCorrMask] = corrAnalysis.cleanCorrMask(finalData,MLCorrMask,signThreshold);

subplot(1,2,2)
imagesc(cleanCorrMask)
title('Cleaned ML output')
axis image
colormap('jet')

%Ask Pavel if it is better to pre-process data so to handle noise or if he
%thinks that pure ML approach can be done?








