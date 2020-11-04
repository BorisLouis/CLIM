function [neighbor] = findNeighbor(idx,dim,r)
    %function to find neighboring pixel in a given radius of a
    %central pixel given by idx. The function also makes sure that
    %we do not have indices outside the images.
    iIdx = idx(1)-r:idx(1)+r;
    jIdx = idx(2)-r:idx(2)+r;

    iIdx = iIdx(iIdx>=1 & iIdx<=dim(1));
    jIdx = jIdx(jIdx>=1 & jIdx<=dim(2));

    neighbor = combvec(iIdx,jIdx)';    
end