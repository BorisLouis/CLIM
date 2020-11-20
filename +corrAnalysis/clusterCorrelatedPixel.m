function [clust,clustEval] = clusterCorrelatedPixel(distanceMap, varargin)

    %parse user input
    narginchk(1,inf)

    %Parse input and give default values
    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('replicate', 5, @(x) isnumeric(x) && x>0);
    params.addParameter('GPU', false, @(x) x == round(x) && or(x==1, x == 0));
    params.addParameter('clust2Test', [1,10], @(x) isvector(x),length(x)==2);
    
    params.parse(varargin{:});

    %Extract values from the inputParser
    clust2Test =  params.Results.clust2Test;
    replicate =  params.Results.replicate;
    GPU = params.Results.GPU;   
    
    
    clust = zeros(size(distanceMap,1),max(clust2Test));
    
    if GPU
        distanceMap = gpuArray(distanceMap);
    end
    
    for i=clust2Test(1):clust2Test(2)
    clust(:,i) = kmeans(distanceMap,i,'emptyaction','drop',...
            'replicate',replicate);
    end
    
    clustEval = evalclusters(distanceMap,clust,'CalinskiHarabasz');
    
end