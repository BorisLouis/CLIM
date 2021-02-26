function [distanceMap] = getDistanceMapFromPxList(inds,data)
    
     %calculate distancemap
    [n,p] = ind2sub(size(data),inds);
    pxIntList = zeros(length(n),size(data,3),'single');
    for i =1:length(n)

        pxIntList(i,:) = single(data(n(i),p(i),:));

    end
    
    clear data inds n p
    distanceMap = single(1-corrcoef(pxIntList'));