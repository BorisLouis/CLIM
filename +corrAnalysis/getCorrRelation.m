function [corrRelation] = getCorrRelation(data2Cluster,r)
    %function to find correlation relation between a each pixel of
    %an image and its neighbor pixels
    
    corrRel  = cell(size(data2Cluster,1),size(data2Cluster,2));
    corrVal  = cell(size(data2Cluster,1),size(data2Cluster,2));
    %loop through pixels
    for i = 1:size(data2Cluster,1)
        for j = 1:size(data2Cluster,2)
            %take a pixel
            currentPxCoord = [i,j];
            %find its neighbor
            neighbor = corrAnalysis.findNeighbor(currentPxCoord,size(data2Cluster),r);
            
            corr = zeros(size(neighbor,1),1);                 
            data1 = squeeze(data2Cluster(currentPxCoord(1),currentPxCoord(2),:));
            %calculate correlation 1-pearson coefficient ==> 0 is
            %correlated 2 is anti correated 1 is uncorrelated
            for k = 1:size(neighbor,1)
                data2 = squeeze(data2Cluster(neighbor(k,1),neighbor(k,2),:));
                tmpCorr = corrcoef(data1,data2);
                corr(k) = 1-tmpCorr(2,1);

            end
            %if correlation between pixel is sufficent we keep
            %track of those pixel as being correlated to the
            %current pixel
            [neighborIdx] = sub2ind(size(corrRel),neighbor(:,1),neighbor(:,2));
            currPxIdx     = sub2ind(size(corrRel),currentPxCoord(:,1),currentPxCoord(:,2));
            corr(neighborIdx==currPxIdx) = [];
            neighborIdx(neighborIdx==currPxIdx) = [];

%             if all(corr>corrThreshold)
% 
%             else
%                 idx = neighborIdx(corr<corrThreshold,:);
% 
%                 if ~isempty(idx)
%                    %convert to indices for simplicity later
%                    corrRel{currPxIdx} = idx;
%                 end
%             end
            corrRel{currPxIdx} = neighborIdx;
            corrVal{currPxIdx} = corr;
        end
    end
    
    %here we test the thresholds
    disp('stop');
    disp('Searching for optimal threshold');
  
    area = zeros(1,30);
    meanCorr = zeros(1,30);
    threshold = zeros(1,30);
    startThreshold = 1;
    for i = 1:50
       
       
       %find cell with value bigger than threshold
       idx = cellfun(@(x) x<startThreshold,corrVal,'UniformOutput',false);
       test = cellfun(@(x) sum(x)>0,idx);
       area(i) = sum(test(:));
       
       currentData = corrVal(test);
       
       meanCorr(i) = mean(cellfun(@(x) mean(x(x<startThreshold)),currentData));
       
       threshold(i) = startThreshold;
       startThreshold = startThreshold - 0.02;
       
    end
    
    %find minimal threshold
    ipt = findchangepts(diff(area));
    %We take the a little bit after the change point to make sure the bckg
    %is fully killed
    minThresh = threshold(ipt(1)+2);
    
    %extract the part of the data after background removal
    idx2Thresh = threshold<=minThresh;    
    minArea = area(idx2Thresh);
    [~,idx] = min(diff(smooth(minArea)));
    %get the corresponding threshold, the part with bkg was not removed,
    %intentionally because the minimum represent the moment where the area
    %change the most, we want to stop scanning the threshold before that
    %acceleration.
    maxThresh = threshold(1+idx);
    
    tRange = [minThresh maxThresh];
    
    %get index of pixels that are not at least correlated to their neighbor 
    %with min thresh
    
    [idx2Delete] = cellfun(@(x) ~any(x<minThresh),corrVal);
    SE = strel('disk',3);
    idx2Delete = imclose(idx2Delete,SE);
        
    idx2Delete = reshape(idx2Delete,size(corrRel,1)*size(corrRel,2),1);
    %reshape and store data.
    corrRelation.listPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
    corrRelation.listVal = reshape(corrVal,size(corrRel,1)*size(corrRel,2),1);
    corrRelation.sumPx  = cellfun(@sum,corrRelation.listVal);
    corrRelation.tRange = tRange;
    
    corrRelation.indPx    = (1:length(corrRelation.listPx))';
        
    %#4 Clean data by keeping only pixel that have correlation
    %relation idx2Delete is now calculated above (~line94)
    %idx2Delete = cellfun(@isempty,corrRelation.listCorrPx);
    corrRelation.listPx(idx2Delete) =[];
    corrRelation.listVal(idx2Delete) = [];
    corrRelation.indPx(idx2Delete) = [];
    corrRelation.sumPx(idx2Delete) = [];  
    
    disp('======> DONE <=======');
end