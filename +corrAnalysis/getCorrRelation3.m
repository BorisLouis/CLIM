function [listCorrPx, corrProd] = getCorrRelation3(data,r,corrThreshold)
    %function to find correlation relation between a each pixel of
    %an image and its neighbor pixels

    corrRel  = cell(size(data,1),size(data,2));
    corrProd  = zeros(size(data,1),size(data,2));
    %loop through pixels
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            
            %take a pixel
            currentPxCoord = [i,j];
            %find its neighbor
            neighbor = corrAnalysis.findNeighbor(currentPxCoord,size(data),r);
              
            %calculate correlation 1-pearson coefficient ==> 0 is
            %correlated 2 is anti correated 1 is uncorrelated
            
            tmpTrace = zeros(size(neighbor,1),size(data,3));
            for k = 1:size(neighbor,1)
                tmpTrace(k,:) = squeeze(data(neighbor(k,1),neighbor(k,2),:));
            end
            corrMat = corrcoef(tmpTrace');
            
            idx = and(ismember(neighbor(:,1),currentPxCoord(1)),ismember(neighbor(:,2),currentPxCoord(2)));
            corr = 1-corrMat(:,idx);
           
            %if correlation between pixel is sufficent we keep
            %track of those pixel as being correlated to the
            %current pixel
            [neighborIdx] = sub2ind(size(corrRel),neighbor(:,1),neighbor(:,2));
            currPxIdx     = sub2ind(size(corrRel),currentPxCoord(:,1),currentPxCoord(:,2));
            corr(neighborIdx==currPxIdx) = [];
            neighborIdx(neighborIdx==currPxIdx) = [];
            
            %calculate correlation metric 
            corrProd(currPxIdx) =  abs(prod(unique(corrMat))/length(unique(corrMat)));
            %if none are correlated enough, we do not do anything
            if all(corr>corrThreshold)

            else
                idx = neighborIdx(corr<corrThreshold,:);

                if length(idx)>1
                   %test that pixel correlated to centre pixel are correlated
                   %together
                   corrIdx = find(corr<corrThreshold);
                   %get all combination of 2
                   idx2CorrMat = nchoosek(corrIdx,2);

                   idx2CorrMat = sub2ind(size(corrMat),idx2CorrMat(:,1),idx2CorrMat(:,2));
                   %we convert corrMat to distance here
                   corrVal2Test = 1-corrMat(idx2CorrMat);

                   %find which pixel are not correlated to the others
                   if any(corrVal2Test>corrThreshold)
                   else
                       %store indices to correlated pixels
                       corrRel{currPxIdx} = idx;
                   end
                   
               else
                   corrRel{currPxIdx} = idx;
               end
                
            end
        end
    end
    
    %clean corrRel
    listCorrPx = reshape(corrRel,size(corrRel,1)*size(corrRel,2),1);
    
    for i=1:length(listCorrPx)
       if ~isempty(listCorrPx{i})
           
           idx2Test = listCorrPx{i};
           
           idx2Delete = cellfun(@isempty,listCorrPx(idx2Test));
           
           listCorrPx{i}(idx2Delete) = [];
           
       end
        
    end   
end