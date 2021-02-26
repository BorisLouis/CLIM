function [cleanCorrMask] = clusterCleanUpMemEff(corrMask,clusters,data)
    threshold = 0.6;
    nCluster = max(corrMask(:));
    
    sizes = cellfun(@size,clusters,'UniformOutput',false);
    sizes = cell2mat(sizes');
    
    [~,idx] = sort(sizes(:,1));
    
    clusters = clusters(idx);
    
    for i = 1:nCluster
        currCluster = clusters{i};
        nPx = size(currCluster,1);
        px2Move = cell(size(currCluster,1),1);
        
        %get distance map for currently treated cluster
        currId = currCluster(:,1);
        currDistMap = corrAnalysis.getDistanceMapFromPxList(currId,data);
        
        n = 1;
        for j = 1:nPx
            currentPx = currCluster(j,1);
            idx2Px    = currCluster(j,2);
%   
            correlation2Cluster = mean(nonzeros(currDistMap(currId==currentPx,:)));
            
            %get mean correlation to cluster
            corr2OtherClusters = zeros(1,nCluster);
            for k = 1:nCluster
               if ~isempty(clusters{k})
                   cluster2Test = clusters{k}(:,1);
                   cluster2Test = [cluster2Test;currentPx];
                   subDistMap   = corrAnalysis.getDistanceMapFromPxList(cluster2Test,data);
%             
                   corr2OtherClusters(k) = mean(nonzeros(subDistMap(end,:)));
                   
               else
                   
                   corr2OtherClusters(k) = 10;
                   
               end
            end
            
            %check if pixel was better correlated to other cluster
            checkRes = correlation2Cluster> min(corr2OtherClusters);
            
            if or(checkRes,nPx==1)
                %if single pixel we check if another cluster is decently
                %correlated to it (compare with initial threshold)
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