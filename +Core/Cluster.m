classdef Cluster < handle
    
    properties (SetAccess = 'private')
        pxList
        corrVal
        stats
        stdCrit = 0.1;
    end
    
    methods
        function obj = Cluster(pxList,corrVal)
            assert(and(length(pxList)==1,length(corrVal)==1),'Initialisation of cluster expected to start with a single data point');
            % store the inital point in the cluster
            obj.pxList = pxList;
            obj.corrVal = corrVal;
            % initialize statistical value for the cluster
            stat.mean = corrVal;
            stat.std = 0;
            
            obj.stats = stat;
        end
        
        function updateStats(obj)
            corrV = obj.corrVal;
            
            stat.mean = mean(corrV);
            stat.std  = std(corrV);
            
            obj.stats = stat;
        end
        
        function addPxToList(obj,pxList,valList)
            currPxList  = obj.pxList;
            currValList = obj.corrVal;
            newPxList   = [currPxList;pxList];
            newValList  = [currValList;valList];
            
            obj.pxList = newPxList;
            obj.corrVal = newValList;
            
            obj.updateStats;
            
        end
        
        
        
        
        
    end
    
end
