function [evalClust] = testNumberOfMLCluster(distanceMap,inds,data,varargin)

    %parse user input
    narginchk(1,inf)

    %Parse input and give default values
    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('replicate', 1, @(x) isnumeric(x) && x>0);
    params.addParameter('GPU', false, @(x) x == round(x) && or(x==1, x == 0));
    params.addParameter('clust2Test', 5, @(x) isnumeric(x) && length(x)==1);
    params.addParameter('deltaClust',2,  @(x) isnumeric(x) && length(x)==1);
    params.parse(varargin{:});

    %Extract values from the inputParser
    clust2Test =  params.Results.clust2Test;
    replicate =  params.Results.replicate;
    GPU = params.Results.GPU;   
    deltaClust = params.Results.deltaClust;
    
    nClust = length(clust2Test);
    clust = zeros(size(distanceMap,1),nClust);
    
    if GPU
        try
            distanceMap = gpuArray(distanceMap);
            clust = gpuArray(clust);
        catch
        end
    end
    
    DCV = 1;
    prevCV = 100;
    counter = 1;
    CV = zeros(1,100);
    avg = zeros(1,100);
    stds = zeros(1,100);
    nClustTested = zeros(1,100);
    while DCV > 0
        if counter >200
            break;
        end
        nClustTested(counter) = clust2Test;
        clust = kmeans(distanceMap,clust2Test,'emptyaction','drop',...
            'replicate',replicate);
        
        MLCorrMask = zeros(size(data(:,:,1)));
        for i = 1:length(inds)
            MLCorrMask(inds(i)) = clust(i);
        end

        [~,relNum] = corrAnalysis.evalClusters(MLCorrMask,data);
        
        CV(counter) = relNum.varCoeff;
        avg(counter) = relNum.meanCorr;
        stds(counter)         = relNum.std;
        
        
        
        DCV = abs(CV(counter) - prevCV);
        
        prevCV = CV(counter);
        
        %update counter
        counter = counter+1;
        clust2Test = clust2Test+deltaClust;
        
        
    end
    
    evalClust.CV = CV;
    evalClust.avg = avg;
    evalClust.std = stds;
    evalClust.nClust = nClustTested;


end