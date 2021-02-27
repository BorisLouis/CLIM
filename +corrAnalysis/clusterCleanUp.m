function [cleanCorrMask] = clusterCleanUp(corrMask,clusters,distanceMap,threshold)
    nCluster = max(corrMask(:));
    
    sizes = cellfun(@size,clusters,'UniformOutput',false);
    sizes = cell2mat(sizes');
    
    [~,idx] = sort(sizes(:,1));
    
    clusters = clusters(idx);
    
    for i = 1:nCluster
        currCluster = clusters{i};
        nPx = size(currCluster,1);
        px2Move = cell(100,1);
        n = 1;
        for j = 1:nPx
            currentPx = currCluster(j,1);
            idx2Px    = currCluster(j,2);
            pxVec = repmat(idx2Px,nPx,1);
            % get mean correlation to cluster
            try
            idx   = sub2ind(size(distanceMap),pxVec,currCluster(:,2));
            catch
                disp('outch')
            end
            correlation2Cluster = mean(distanceMap(idx));
            
            %get mean correlation to cluster
            corr2OtherClusters = zeros(1,nCluster);
            for k = 1:nCluster
               if ~isempty(clusters{k})
                   cluster2Test = clusters{k}(:,2);
                   nPx2 = size(cluster2Test,1);
                   pxVecTest = repmat(idx2Px,nPx2,1);
                   idx2Test   = sub2ind(size(distanceMap),pxVecTest,cluster2Test);

                   subDistMap = distanceMap(idx2Test);

                   corr2OtherClusters(k) = mean(nonzeros(subDistMap));
               else
                   corr2OtherClusters(k) = 10;
               end
            end
            
            %check if pixel was better correlated to other cluster
            checkRes = correlation2Cluster>=min(corr2OtherClusters);
            
            if or(checkRes,nPx==1)
                if nPx == 1 && min(corr2OtherClusters)> threshold
                    clusters{i} =[];
                else
                    %add pixel to be moved to the list
                    [~,idx] = min(corr2OtherClusters);

                    px2Move{n} = [currentPx, idx2Px, idx];
                    n = n+1;
                end
            else
                
            end
        end
        idx2Del = cellfun(@isempty,px2Move);
        px2Move(idx2Del) = [];
 
        % Move pixel from cluster
        for j = 1:length(px2Move)
            currPx = px2Move{j};
            if currPx(3) ~= i
                clusters{currPx(3)} = [clusters{currPx(3)}; currPx(1:2)];
            end
            
        end
        
    end
        %TODO clean up singleton cluster
        idx2Delete = cellfun(@isvector,clusters);
        clusters(idx2Delete) = [];
        idx2Delete = cellfun(@isempty,clusters);
        clusters(idx2Delete) = [];
        
        newCorrMask = zeros(size(corrMask));
        
        for i = 1 : length(clusters)
            currCluster = clusters{i}(:,1);
            tmpMask = zeros(size(corrMask));
            tmpMask(currCluster) = i;
            tmpMask = imfill(tmpMask,'holes');
            newCorrMask(tmpMask>0) = i;
            
        end
        
        %clear up 'holes' in the numbers
        nClusters = max(newCorrMask(:));
        cleanCorrMask = zeros(size(corrMask));
        for i = 1:nClusters

            BWCopy = newCorrMask==nClusters - (i-1);
            
            if ~isempty(BWCopy)
                cleanCorrMask(BWCopy) = max(cleanCorrMask(:))+1;
            end
        end

        
        
        
end