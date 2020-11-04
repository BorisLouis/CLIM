function [corrRel] = getCorrRelation(data2Cluster,r,corrThreshold)
    %function to find correlation relation between a each pixel of
    %an image and its neighbor pixels

    corrRel  = cell(size(data2Cluster,1),size(data2Cluster,2));
    %loop through pixels
    for i = 1:size(data2Cluster,1)
        for j = 1:size(data2Cluster,2)
            %take a pixel
            currentPxCoord = [i,j];
            %find its neighbor
            neighbor = Core.CorrClusterMovie.findNeighbor(currentPxCoord,size(data2Cluster),r);

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

            if all(corr>corrThreshold)

            else
                idx = neighborIdx(corr<corrThreshold,:);

                if ~isempty(idx)
                   %convert to indices for simplicity later
                   corrRel{currPxIdx} = idx;
                end
            end
        end
    end
end