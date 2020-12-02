function [corrMask, hierarchical] = cleanMLCorrMask(data,corrMask,cleanThresh)

    data = double(data);
    hierarchical = cell(size(corrMask));
    for i = 1:max(corrMask(:))
         
        %get the requested cluster
        subMask = corrMask == i;
        subMask = bwareaopen(subMask,8);
        subMask = imfill(subMask,8,'holes');

        currMask = bwlabel(subMask);

        idx = max(currMask(:));

        for j = 1:idx
            subMask = currMask==j;
            % get indices of traces 
            [row,col] = find(subMask);
            r = repmat(row,1,size(data,3));
            r = reshape(r',size(data,3)*length(row),1);
            c = repmat(col,1,size(data,3));
            c = reshape(c',size(data,3)*length(col),1);

            f = repmat((1:size(data,3))',length(col),1);
            idx = sub2ind(size(data),r,c,f);

            %get the data
            tmpTrace = data(idx);
            tmpTrace = reshape(tmpTrace,size(data,3),length(row));

            %check correlation
            tmpCorr = corrcoef(double(tmpTrace));
            %test for pixel with mostly correlation <0.3
            test = tmpCorr<0.3;
            test = sum(test,1);

            %calculate the percentage of uncorrelated pixel in each
            %col
            test = test/length(tmpCorr);
            % if 40% of the correlation is lower than the threshold we
            
            idx2Delete = find(test>cleanThresh);
            idx2Delete = [row(idx2Delete), col(idx2Delete)];
            idx2Delete = sub2ind(size(subMask),idx2Delete(:,1),idx2Delete(:,2));
            % delete the pixels
            subMask(idx2Delete) = 0;
            
            %store the cleaned mask
            corrMask(corrMask==i) =0;
            corrMask(subMask>0) = i;
            %Keep track of clusters hierarchy
            nElem = ones(length(subMask(subMask>0)),1);
            idxMat = [nElem*i,nElem*j];
           
            hierarchical(subMask>0) = mat2cell(idxMat,nElem);
            
        end

        
    end





end