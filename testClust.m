%% Simulate correlated data

%% Simulation input
sizeIm = [512 512];
nFrames = 300;
data = zeros(sizeIm(2),sizeIm(1),nFrames);
nParticles = 100;
corrThreshold = 0.3;%smaller is more selective here (0 is perfect correlation)
model.name = 'gaussian';
model.sigma_x = 5;
model.sigma_y = 5;
r = 2; %radius for checking correlation


%% Simulations
%coord=[56,56;72,72;56,72;72,56];
% coord = [16,32;48,32];
[X,Y] = meshgrid(1:sizeIm(1),1:sizeIm(2));
%coord = zeros(nParticles,2);
for i = 1:nParticles
    
%     x0 = coord(i,1);
%     y0 = coord(i,2);
    x0 = randperm(sizeIm(1),1);
    y0 = randperm(sizeIm(2),1);
    c  = 0;
    BaseInt = 500 + i * 100;
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
        
        PSF = getPSF(X,Y,x0,y0,int,model);
        
        
        data(:,:,j) = data(:,:,j) + PSF; 


    end
    coord(i,:) = [x0,y0];   
end

noise = randn(size(data));

finalData = data + ones(size(data))*100 +noise*20;




%%
frameRate = 10;
filename = 'cover.gif';
% Capture the plot as an image 
h = figure('Position',[10 10 768 768]);
for i = 1:size(finalData,3)
      imagesc(finalData(:,:,i));
      colormap('hot')
      caxis([0 max(finalData(:))]);
      
      axis tight
      set(gca,'XColor', 'none','YColor','none')
      a = gca;
      a.XTickLabel = [];
      a.YTickLabel = [];
      
      drawnow;
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      
      % Write to the GIF File 
      
      if i == 1

          imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

      else

          imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

      end
end

%% save a movie as example
vidFile = VideoWriter('rawMov.mp4','MPEG-4');
vidFile.FrameRate = 10;
open(vidFile);
figure
for i = 1:size(finalData,3)
   imagesc(finalData(:,:,i));
   colormap('hot')
   caxis([0 max(finalData(:))]);
   axis image
   set(gca,'XColor', 'none','YColor','none')
   drawnow;
   im = getframe;
   
   writeVideo(vidFile,im);
   clf
end
close(vidFile)

%% Normalization of data
% #1 Normalize the data
data2Cluster = finalData-min(finalData,[],3);
data2Cluster = data2Cluster./max(data2Cluster,[],3);
% # Question: Will the relative amplitude of blinking affect the data?
% Check normalization !!!

%% Clustering

%High degree: If the coefficient value lies between ± 0.50 and ± 1, then it
%is said to be a strong correlation. Moderate degree: If the value lies 
%between ± 0.30 and ± 0.49, then it is said to be a medium correlation. 
%Low degree: When the value lies below + . 29, then it is said to be a small
%correlation. ==>https://www.statisticssolutions.com/pearsons-correlation-coefficient/

%pre-process the data to find pixel with correlation relation
[corrRel]  = corrAnalysis.getCorrRelation(data2Cluster,r,corrThreshold);
listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
inds    = (1:length(listCorrPx))';

idx2Delete = cellfun(@isempty,listCorrPx);

listCorrPx(idx2Delete) =[];

inds(idx2Delete) = [];

[n,p] = ind2sub(size(corrRel),inds);
pxIntList = zeros(length(n),size(data2Cluster,3));
for i =1:length(n)
   
    pxIntList(i,:) = data2Cluster(n(i),p(i),:);
    
end

% get distance map
dist = pdist2([n,p],[n,p]);
%normalize the distance scale
dist = dist./max(dist(:));
distanceMap = 1-corrcoef(pxIntList');

% %#1 Try clustering the data together.
% testMask = zeros(size(data2Cluster,1),size(data2Cluster,2));
% T1 = clusterdata(distanceMap,10);
% for i = 1: length(T1)
%     testMask(inds(i)) = T1(i);
%     
% end
% 
% figure
% imagesc(testMask);
% 
% 
% %#2 DBScan is very similar to what I implemented as pseudo clustering so it
% %fails for similar reasons.
% idx = dbscan(distanceMap,24,50,'Distance','precomputed');
% unique(idx)
% 
%#3 improved pseudo clustering


%% Pseudo clustering ==> Generate correlation mask
dim = size(data2Cluster);
distanceMap = corrAnalysis.getDistanceMapFromPxList(inds,data2Cluster);
[corrMask] = corrClustering(listCorrPx,inds,distanceMap,dim);

figure
imagesc(corrMask)
hold on
scatter(coord(:,1),coord(:,2),20,repmat([1 0 0],size(coord,1),1),'filled')
contour = bwboundaries(corrMask);
        
%% Movie + contour        
vidFile = VideoWriter('contourMov.mp4','MPEG-4');
vidFile.FrameRate = 10;
open(vidFile);
figure
for i = 1:size(finalData,3)
   imagesc(finalData(:,:,i));
   hold on
   for j = 1:length(contour)
       boundary = contour{j};
       plot(boundary(:,2),boundary(:,1),'w','LineWidth',2)
   end
   colormap('hot')
   caxis([0 max(finalData(:))]);
   axis image
   drawnow;
   im = getframe;
   writeVideo(vidFile,im);
   clf
end
close(vidFile)
        


%%
function [neighbor] = findNeighbor(idx,dim,r)

    iIdx = idx(1)-r:idx(1)+r;
    jIdx = idx(2)-r:idx(2)+r;
    
    iIdx = iIdx(iIdx>=1 & iIdx<=dim(1));
    jIdx = jIdx(jIdx>=1 & jIdx<=dim(2));
   % neighbor = [iIdx(:),jIdx(:)];
    neighbor = combvec(iIdx,jIdx)';    
    
%     mem = ismember(neighbor,idx);
%     id = sum(mem,2);
%     
%     neighbor(id==2,:) = [];
end


function [ psf ] = getPSF( Xgrid, Ygrid, xpos, ypos, A, model )
%getPSF calculates the PSF of an emitter using the desired model
%   Detailed explanation goes here

if strcmp(model.name, 'gaussian')
    sigma_x = model.sigma_x;
    sigma_y = model.sigma_y;
    

    f1 = ((Xgrid - xpos).^2) / (2*(sigma_x.^2));
    f2 = ((Ygrid - ypos).^2) / (2*(sigma_y.^2));
    %A  = 1/(2*pi*sigma_x*sigma_y);  
    psf = A * exp(-(f1+f2));

end

end

function [corrRel] = getCorrRelation(data2Cluster,r,corrThreshold)
    corrRel  = cell(size(data2Cluster,1),size(data2Cluster,2));
    for i = 1:size(data2Cluster,1)
        for j = 1:size(data2Cluster,2)
            %take a pixel
            currentPxCoord = [i,j];
            %find its neighbor
            neighbor = findNeighbor(currentPxCoord,size(data2Cluster),r);

            corr = zeros(size(neighbor,1),1);

            %calculate correlation
            data1 = squeeze(data2Cluster(currentPxCoord(1),currentPxCoord(2),:));

            for k = 1:size(neighbor,1)
                data2 = squeeze(data2Cluster(neighbor(k,1),neighbor(k,2),:));

                tmpCorr = corrcoef(data1,data2);

                corr(k) = 1-tmpCorr(2,1);

            end
            if all(corr>corrThreshold)

            else

                idx = neighbor(corr<corrThreshold,:);
                mem = ismember(idx,currentPxCoord);
                id = sum(mem,2);
                idx(id==2,:) = [];

                if ~isempty(idx)
                   %convert to indices for simplicity later
                   corrRel{currentPxCoord(1),currentPxCoord(2)} = sub2ind(size(corrRel),idx(:,1),idx(:,2));

                end

            end
        end
    end
        
end

function [corrMask] = corrClustering(listCorrPx,inds,distanceMap,dim)
    safeCount = 2*dim(1)*dim(2);
    count = 1;
    group = 1;
    corrMask = zeros(dim(1),dim(2));
    treatedIdx = NaN(length(inds),1);
    indsCopy = inds;
    while ~isempty(listCorrPx) 

        currList = listCorrPx{1};
        currIndex = inds(1);
        currList(ismember(currList,treatedIdx,'row'),:) =[];
        if isempty(currList)
            listCorrPx(1) = [];
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;
            inds(1) = [];
            %potential fix, give pixel that should have been treated the same
            %group as the closest pixel.
        else
            while ~isempty(currList)
                % add the new pixel to the group and keep track of it
                corrMask(currIndex) = group;
                idx = find(isnan(treatedIdx(:,1)),1);
            %    [treatedIdx(idx,1),treatedIdx(idx,2)] = ind2sub(dim(1:2),currIndex);
                treatedIdx(idx) = currIndex;
                % remove it from the list
                listCorrPx(inds==currIndex) = [];
                inds(inds==currIndex) = [];

                %update the currentlist to add the new data
                currIndex = currList(1);
                currList(1,:) = [];
                %get list of element correlated to the current element
                list2Add  = listCorrPx{inds==currIndex}; 
                %remove already treated cases and cases already in the list
                list2Add(or(ismember(list2Add,treatedIdx,'row'),ismember(list2Add,currList,'row'))) = []; 
                
                %find pixel that are already inside the cluster to check
                %for correlation
                currCluster = find(corrMask==group);
                
                try
                    r = randperm(length(currCluster,10));
                    currCluster = currCluster(r);
                catch
                    
                end
                
                [list2Add] = checkCorrList(list2Add,currCluster,distanceMap,indsCopy);
                              
                currList  = [currList;list2Add];
                
                if any(ismember(currList,treatedIdx))
                   disp('oups') 
                end

                if isempty(currList)
                    disp(['Found ' num2str(group) ' group(s).']);
                end

                count = count+1;
            end
            %mark the last case as treated:
            corrMask(currIndex) = group;
            idx = find(isnan(treatedIdx(:,1)),1);
            treatedIdx(idx) = currIndex;
            
            % remove it from the list
            listCorrPx(inds==currIndex) = [];
            inds(inds==currIndex) = [];

            %if exit the first while loop, the first group is complete, we need to
            %increment the group number as we are supposed to have treated all
            %cases.
            group = group+1;
        
        end

        count=count+1;

        if count >= safeCount
            error('something went wrong')
        end
    end   
end

function [ind2Add] = checkCorrList(ind2Add,currCluster,distanceMap,indsCopy)
%function to enforce consistency within the cluster by checking that
%to-be-added element are correlated with a random subset of the cluster
nInd = length(ind2Add);
nClust = length(currCluster);
newList = ind2Add;
if ~isempty(currCluster)
    %reshape input to compare correlation
    id2Add  = find(ismember(indsCopy,newList));
    newList = repmat(id2Add,1,length(currCluster))';
    newList = newList(:);
    currCluster = find(ismember(indsCopy,currCluster));
    currCluster = repmat(currCluster,nInd,1);
       
    %combine input to extract value from distance map
    subs = [newList,currCluster];    
    inds = sub2ind(size(distanceMap),subs(:,1),subs(:,2));
    corrVal = distanceMap(inds);
    %reshape corrVal to test for each cluster
    corrVal = reshape(corrVal,nClust,nInd);
    %we are more relaxed since pixel far apart might be tested
    corrTest = corrVal<0.5;
    %we want new data point to be correlated to 80% of the subset from the
    %existing group
    idx2delete = sum(corrTest,1)./size(corrTest,1)<0.7;
    if sum(idx2delete)>0
        disp('oh');
    end
    %delete the data not matching the requested threshold
    ind2Add(idx2delete) = [];    
    ind2Add = unique(ind2Add);
   
end

end



