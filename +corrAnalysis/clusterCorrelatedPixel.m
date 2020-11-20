function [clust,clustEval] = clusterCorrelatedPixel(distanceMap, clust2Test)

    %TODO: Increase input so GPU or other parameter can be selected outside
    %of the main function

    if length(clust2Test)<2
        clust2Test = [1 clust2Test];
    elseif length(clust2Test) == 2
    else
        error('unexpected format for cluster to test');
    end
    
    clust = zeros(size(distanceMap,1),max(clust2Test));
    for i=clust2Test(1):clust2Test(2)
    clust(:,i) = kmeans(distanceMap,i,'emptyaction','singleton',...
            'replicate',5);
    end
    clustEval = evalclusters(distanceMap,clust,'CalinskiHarabasz');
    
end