function [distanceMap] = getDistanceMapFromPxList(inds,data)
    
    %calculate distancemap
    [n,p] = ind2sub(size(data),inds);
    pxIntList = zeros(length(n),size(data,3));
    for i =1:length(n)

        pxIntList(i,:) = data(n(i),p(i),:);

    end
    distanceMap = 1-corrcoef(pxIntList');