classdef CorrCluster
    % CorrClusterSummary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nPx
        meanCorr
        minCorr
        avgTrace
        inds
    end
    
    methods
        function obj = CorrCluster(inds,trace)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.inds = inds;
            obj.avgTrace = trace;
            obj.minCorr = -inf;
            obj.meanCorr = -inf;
            obj.nPx = length(inds);
            
        end
        
        function merge(obj,cluster)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
            
        end
    end
end

