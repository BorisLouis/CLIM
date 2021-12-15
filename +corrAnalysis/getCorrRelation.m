function [corrRelation] = getCorrRelation(data2Cluster,r,neigh)
    %function to find correlation relation between a each pixel of
    %an image and its neighbor pixels
    
    corrRel  = cell(size(data2Cluster,1),size(data2Cluster,2));
    corrVal  = cell(size(data2Cluster,1),size(data2Cluster,2));
    
    corrMap  = zeros(size(data2Cluster,1),size(data2Cluster,2));
   
    %loop through pixels
    for i = 1:size(data2Cluster,1)
        for j = 1:size(data2Cluster,2)
            %take a pixel
            currentPxCoord = [i,j];
            
            currPxIdx     = sub2ind(size(corrRel),currentPxCoord(:,1),currentPxCoord(:,2));
            if or(data2Cluster(i,j,1)==0,isnan(data2Cluster(i,j,1)))
                corrMap(i,j) = 0;
                pValMap(i,j) = 0;
                corrRel{currPxIdx} = [];
                
            else
                %find its neighbor
                neighbor = corrAnalysis.findNeighbor(currentPxCoord,size(data2Cluster),r,neigh);

                corr = zeros(size(neighbor,1),1);
                pVal = zeros(size(neighbor,1),1);
                data1 = squeeze(data2Cluster(currentPxCoord(1),currentPxCoord(2),:));
                %calculate correlation 1-pearson coefficient ==> 0 is
                %correlated 2 is anti correated 1 is uncorrelated
                for k = 1:size(neighbor,1)
                    data2 = squeeze(data2Cluster(neighbor(k,1),neighbor(k,2),:));
                    [tmpCorr] = corrcoef(data1,data2);
                    corr(k) = tmpCorr(2,1);

                end
                %if correlation between pixel is sufficent we keep
                %track of those pixel as being correlated to the
                %current pixel
                [neighborIdx] = sub2ind(size(corrRel),neighbor(:,1),neighbor(:,2));
               
                corr(neighborIdx==currPxIdx) = [];
                neighborIdx(neighborIdx==currPxIdx) = [];

                corrMap(i,j) = nanmean(corr);
                

                corrRel{currPxIdx} = neighborIdx(~isnan(corr));
                corrVal{currPxIdx} = corr(~isnan(corr));
            end
        end
    end
    %we use 0.2 to kill the background and low correlation pixels
    [idx2Delete] = cellfun(@(x) all(x<0.2),corrVal);
    %the smoothing of idx2delete cause inconsistencies in the list of
    %correlated pixel so we dont use it anymore
    %SE = strel('disk',3);
    %idx2Delete = imclose(idx2Delete,SE);
        
    idx2Delete = reshape(idx2Delete,size(corrRel,1)*size(corrRel,2),1);
    %reshape and store data.
    corrRelation.listPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
    corrRelation.listVal = reshape(corrVal,size(corrRel,1)*size(corrRel,2),1);
    corrRelation.meanPx  = cellfun(@mean,corrRelation.listVal);
    
    corrRelation.indPx    = (1:length(corrRelation.listPx))';
        
    %#4 Clean data by keeping only pixel that have correlation
    %relation idx2Delete is now calculated above (~line94)
    %idx2Delete = cellfun(@isempty,corrRelation.listCorrPx);
    corrRelation.listPx(idx2Delete) =[];
    corrRelation.listVal(idx2Delete) = [];
    corrRelation.indPx(idx2Delete) = [];
    corrRelation.meanPx(idx2Delete) = [];  
    corrRelation.corrMap = corrMap;
  
    disp('======> DONE <=======');
    
    %listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
    %sumCorrPx  = cellfun(@sum,listCorrPx);
end