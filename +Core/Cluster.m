classdef Cluster < handle
    
    properties (SetAccess = 'private')
        pxList
        corrVal
        stats
        stdCrit = 0.1;
        ID
        method = 'std'
    end
    
    methods
        function obj = Cluster(pxList,corrVal,id)
            assert(length(pxList)==1,'Initialisation of cluster expected to start with a single data point');
            % store the inital point in the cluster
            obj.pxList = pxList;
            obj.corrVal = corrVal;
            % initialize statistical value for the cluster
            stat.mean = corrVal;
            stat.std = 0;
            
            obj.stats = stat;
            obj.ID    = id;
        end
        
        function updateStats(obj)
            corrV = obj.corrVal;
            
            stat.mean = mean(corrV);
            stat.std  = std(corrV);
            
            obj.stats = stat;
        end
        
        function [addedPx] = addPxToList(obj,pxList,valList)
            currPxList  = obj.pxList;
            currValList = obj.corrVal;
            
            if length(currPxList) == 1
                %if there is only 1 px we add all the neighboring pixels    
                newPxList   = [currPxList;pxList];
                newValList  = [currValList;valList];
                addedPx = pxList;
                
            else
                %if there is more pixels then we check statistic to add or
                %not
                switch obj.method
                    case 'std'
                        currMean = obj.stats.mean;
                        currStd  = 3*obj.stats.std;
%                         if currStd< obj.stdCrit
%                             currStd = obj.stdCrit;
%                         end
                        %threshold is equal to mean +std
                        threshold = currMean+currStd;
                        idx2Add = valList<threshold;
                        addedPx = pxList(idx2Add);

                        newPxList = [currPxList;addedPx];
                        newValList = [currValList; valList(idx2Add)];
                        
                    case ''
                        fe=5;
                        
                end

            end
            obj.pxList = newPxList;
            obj.corrVal = newValList;
            
            obj.updateStats;
            
        end
        
        
        function print(obj)
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['Found cluster number ' num2str(obj.ID)]);
            disp(['Contains: ' num2str(length(obj.pxList)) ' pixels']);
            disp(['Avg Corr: ' num2str(obj.stats.mean,'%.2f')]);
            disp(['Sdev: '     num2str(obj.stats.std,'%.2f')]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            
        end
        
        
    end
    
end
