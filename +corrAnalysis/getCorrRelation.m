function [listCorrPx,sumCorrPx] = getCorrRelation(data2Cluster,r,startThreshold)
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
    totArea = size(corrRel,1)*size(corrRel,2);
    area = zeros(1,30);
    meanCorr = zeros(1,30);
    threshold = zeros(1,30);
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
    
    while(true)
       %Test threhsold
       
        
        
        
    end
    
    
    
    
    
    %reshape
    listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
    sumCorrPx  = cellfun(@sum,listCorrPx);
end