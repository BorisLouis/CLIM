function [MLCorrMask] = getMaskFromMLCluster(clust,inds,clust2Use,dim)
    
    MLCorrMask = zeros(dim);

    for i = 1:length(inds)
        MLCorrMask(inds(i)) = clust(i,clust2Use);
    end
    
end