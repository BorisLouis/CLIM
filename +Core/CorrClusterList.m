classdef CorrClusterList < handle
    % CorrClusterSummary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nCluster
        clusters
        statistics
        
    end
    
    methods
        function obj = CorrClusterList(data,inds)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            clusters = cell(size(inds));
            for i = 1:size(inds)
                
               [row,col] = ind2sub(size(data),inds(i));
               
               clusters{i} = Core.CorrCluster(inds(i),squeeze(data(row,col,:)));
                
                
            end
            
            obj.clusters = clusters;
            obj.nCluster = length(obj.clusters);
            
        end
        
        function updateState(obj)
            nClust = obj.nCluster;
            
            stats = table(zeros(nClust,1),zeros(nClust,1),zeros(nClust,1),'VariableNames',{'minCorr','meanCorr','nPx'});
            for i = 1:nClust
                currClust = obj.clusters{i};
                stats(i,:) = currClust.getStats;               
                
                
            end
            
            obj.statistics = table(mean(stats.minCorr),mean(stats.meanCorr),mean(stats.nPx),'VariableNames',{'minCorr','meanCorr','nPx'});
            
        end
        
        
        
        function merge(obj,indA,indB,distanceMap)
            % Merge 2 clusters 
            clustA = obj.clusters{indA};
            clustB = obj.clusters{indB};
            
            
            
            % get s
            
            inds = clustB.getInds;
            trace = clustB.getTrace;
            
            %add ClustB into ClustA
            clustA.addPixel(inds,trace,distanceMap);
            
            
                        
        end
    end
end

