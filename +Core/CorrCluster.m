classdef CorrCluster <handle
    % CorrClusterSummary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nPx
        avgTrace
        corrRel
        inds
    end
    
    methods
        function obj = CorrCluster(inds,trace)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.inds = inds;
            obj.avgTrace = trace;
            obj.nPx = length(inds);
            obj.corrRel = 1;
            
        end
        
        function inds = getInds(obj)
            inds = obj.inds;
        end
        
        function avgTrace = getTrace(obj)
            
           avgTrace = obj.avgTrace; 
        end
        
        function [stats] =  getStats(obj)
            map = tril(obj.corrRel);
            
            stats = table(min(map(map>0)),mean(map(map>0)),obj.nPx,'VariableNames',{'minCorr','meanCorr','nPx'});           
            
        end
        
        function addPixel(obj,inds,trace,distMap)   
            currTrace = obj.avgTrace;
            
            obj.inds = sort([obj.inds; inds]);
            %weighed average
            obj.avgTrace = (obj.nPx*currTrace + length(inds)*trace)/(obj.nPx+length(inds));
            
            idx = nchoosek(obj.inds,2);
            
            ind = sub2ind(size(distMap),idx(:,1),idx(:,2));
            
            obj.corrRel = distMap(ind);
            
        end
        
    end
end

