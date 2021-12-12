classdef CorrClusterList
    % CorrClusterSummary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nCluster
        meanCorr
        stdMeanCorr
        minCorr
        stdMinCorr
        InterClusterCorr
        clusters
    end
    
    methods
        function obj = CorrCluster(clusters)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.clusters = clusters;
            obj.minCorr = -inf;
            obj.meanCorr = -inf;
            obj.InterCLusterCorr = +inf;
            
            obj.nPx = length(inds);
            
        end
        
        function merge(obj,cluster)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
            
        end
    end
end

