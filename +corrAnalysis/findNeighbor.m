function [neighbor] = findNeighbor(idx,dim,r,neigh)
    %function to find neighboring pixel in a given radius of a
    %central pixel given by idx. The function also makes sure that
    %we do not have indices outside the images.
    
    switch neigh
        case 8
            iIdx = idx(1)-r:idx(1)+r;
            jIdx = idx(2)-r:idx(2)+r;

            iIdx = iIdx(iIdx>=1 & iIdx<=dim(1));
            jIdx = jIdx(jIdx>=1 & jIdx<=dim(2));

            neighbor = combvec(iIdx,jIdx)';
        case 4
            neighbor = [idx(1)-r,idx(2);
                        idx(1)+r, idx(2);
                        idx(1),idx(2);
                        idx(1),idx(2)-r;
                        idx(1),idx(2)+r;
                        ];
            
            idx2Del1 = neighbor(:,1)>=1 & neighbor(:,1)<=dim(1);
            idx2Del2 = neighbor(:,2)>=1 & neighbor(:,2)<=dim(2);
            
            idx2Del = and(idx2Del1,idx2Del2);
            
            neighbor(~idx2Del,:) = [];
            
        otherwise
            error('inconsistent number of neighbor chosen, please use 4 or 8')
    end
end