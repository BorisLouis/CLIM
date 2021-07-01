classdef CorrStatus < handle
    
    properties (SetAccess = 'private')
        listPx 
        listVal   
        sumPx     
        inds
        treatedPx = [];
    end
    
    methods
        function obj = CorrStatus(listCorrPx,listVal,sumPx,inds)
            obj.listPx = listCorrPx;
            obj.listVal = listVal;
            obj.sumPx   = sumPx;
            obj.inds    = inds;
        end
        
        function deletePxFromList(obj,idx2Delete)
            assert(length(idx2Delete) == length(obj.listPx),'list of index to delete needs to be the same size as the list to be deleted from');
            obj.treatedPx = [obj.treatedPx obj.inds(idx2Delete)];
            obj.listPx(idx2Delete) = [];
            obj.listVal(idx2Delete) = [];
            obj.sumPx(idx2Delete) = [];
            obj.inds(idx2Delete) =[];
            
        end
        
        function [res] = isFinished(obj)
            
            if isempty(obj.listPx)
                res = true;
            else
                res = false;
            end            
        end
        
        
    end
    
end
    