%% test splitting on a clear case scenario
clear
close all
clc

%% Load data
folder = 'D:\Documents\Unif\PhD\2020-Data\12 - Dec\Film Blinking\HierarchicalCluster';

dataFile = [folder filesep 'data.mat'];
distMapFile = [folder filesep 'DistMapToSeparate.mat'];
indFile     = [folder filesep 'inds.mat'];

load(dataFile)
load(distMapFile)
load(indFile)

%%
guessFromCorrMat = 574;
figure
imagesc(distanceMap)


theoreticalClust = zeros(size(data(:,:,1)));
theoreticalClust(inds(inds<inds(guessFromCorrMat))) = 1;
theoreticalClust(inds(inds>=inds(guessFromCorrMat))) = 2;

figure
imagesc(theoreticalClust)


%% try different clustering
clust = spectralcluster(distanceMap,2,'Distance','precomputed');

%%

clust = dbscan(distanceMap,2,10,'Distance','precomputed');



%%
[clust,~,sumD,D] = kmeans(distanceMap,2,'replicate',3);

%% create a mask
MLCorrMask = zeros(size(data(:,:,1)));
for i = 1:length(inds)
    MLCorrMask(inds(i)) = clust(i);
end
 figure
imagesc(MLCorrMask)



%% is new clust significantly better

clust1 = MLCorrMask ==1;
clust2 = MLCorrMask ==2;
inds1 = find(clust1);
inds2 = find(clust2);

distMap1 = corrAnalysis.getDistanceMapFromPxList(inds1,data);
distMap2 = corrAnalysis.getDistanceMapFromPxList(inds2,data);

disp(['The median of the cluster distance map ' num2str(median(distanceMap(:)))])
disp(['The median of the subcluster 1 ' num2str(median(distMap1(:)))])
disp(['The median of the subcluster 2 ' num2str(median(distMap2(:)))])


%% further clustering
clust1Lab = bwlabel(clust2);
subClust1 = clust1Lab ==2;
inds = find(subClust1);

subClustDistMap = corrAnalysis.getDistanceMapFromPxList(inds,data);

for i = 2:5
    [subClust(:,i-1),~,sumD,D] = kmeans(subClustDistMap,i,'replicate',3);
end

 clustEval = evalclusters(subClustDistMap,subClust,'CalinskiHarabasz');
 
 MLCorrMask = zeros(size(data(:,:,1)));
for i = 1:length(inds)
    MLCorrMask(inds(i)) = subClust(i,1);
end
 figure
imagesc(MLCorrMask)
 

sClust1 = MLCorrMask ==1;
sClust2 = MLCorrMask ==2;
sInds1 = find(sClust1);
sInds2 = find(sClust2);

sDistMap1 = corrAnalysis.getDistanceMapFromPxList(sInds1,data);
sDistMap2 = corrAnalysis.getDistanceMapFromPxList(sInds2,data);

disp(['The median of the cluster distance map ' num2str(median(subClustDistMap(:)))])
disp(['The median of the subcluster 1 ' num2str(median(sDistMap1(:)))])
disp(['The median of the subcluster 2 ' num2str(median(sDistMap2(:)))])


subClustDistMap(subClustDistMap==0) = nan;
sDistMap1(sDistMap1==0) = nan;
sDistMap2(sDistMap2==0) = nan;

testMainClust = mean(max(subClustDistMap,[],2)-min(subClustDistMap,[],2))
testClust1 = mean(max(sDistMap1,[],2)-min(sDistMap1,[],2))
testClust2 = mean(max(sDistMap2,[],2)-min(sDistMap2,[],2))


 

