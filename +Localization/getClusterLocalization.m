 function [loc,contour] = getClusterLocalization(data,mask,idx)
            images = regionprops(mask,'Image','BoundingBox');
            
            bBox = round(images(idx).BoundingBox);
            currentImage = images(idx).Image;
            
            contour = bwboundaries(currentImage);
            contour = contour{1,1};
            
            contour(:,1) = contour(:,1)+bBox(2);
            contour(:,2) = contour(:,2)+bBox(1);
            
            data2Fit = double(data(bBox(2):bBox(2)+bBox(4)-1,bBox(1):bBox(1)+bBox(3)-1,:));
            
            [X,Y]= meshgrid(bBox(1):bBox(1)+bBox(3)-1,bBox(2):bBox(2)+bBox(4)-1);
            domain(:,:,1) = X;
            domain(:,:,2) = Y;
            
            loc = struct('x',zeros(size(data2Fit,3),1),'y',zeros(size(data2Fit,3),1),...
                'int',zeros(size(data2Fit,3),1),'angle',zeros(size(data2Fit,3),1));
            
            for i = 1:size(data2Fit,3)
                currentFrame = data2Fit(:,:,i);
                currentFrame(~currentImage) = 0;
                
                %fit the region with a Gaussian
                [gPar, fit] = Gauss.Gauss2D_Fit(currentFrame,domain);
%                 figure
%                 surf(currentFrame)
%                 hold on 
%                 surf(fit)
                
                loc.x(i) = gPar(5);
                loc.y(i) = gPar(6);
                loc.int(i) = gPar(1);
                loc.angle(i) = gPar(7);
                
            end 
        end
        
        
       