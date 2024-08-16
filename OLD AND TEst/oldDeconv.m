       
%             %remove background
%             medAutoCorr = zeros(size(signal,1),size(signal,2));
%             for i = 1:size(signal,1)
%                 for j = 1:size(signal,2)
%                 
%                     currData = double(squeeze(signal(i,j,:)));
%                     AC = autocorr(currData);
%                     medAutoCorr(i,j) = AC(2);
%                     
%                 end
%             end
%             
%             deleteMask = abs(medAutoCorr)>decThresh;
%             %clean up the binary image mask
%             SE = strel('disk',4);
%             deleteMask = imopen(deleteMask,SE);
%             deleteMask = imfill(deleteMask,'holes');
%             %keep only the largest area
%             d = regionprops(deleteMask,'Area','PixelIdxList');
%             [~,idx] = max([d.Area]);
%             delMask = zeros(size(deleteMask));
%             delMask(d(idx).PixelIdxList) = 1;
%             
%             %repeeat it in z and multiply it
%             delMask = repmat(delMask,1,1,size(signal,3));
%             
%             bkgCorrData = double(delMask).*double(signal);