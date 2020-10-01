function [correctedStack,Drift]=CorrelationDrift(im_In,scalingFactor,correlationInfo)

corrSz = correlationInfo.corrSz; 
driftPeriod = correlationInfo.driftPeriod;
maxDrift = correlationInfo.maxDrift;
frame1 = im_In(:,:,1);
%get cropping boundaries
[lb,ub] = getCropPar(frame1,corrSz);

%cropping
imCropped = im_In(lb(2): ub(2),lb(1):ub(1),:);

correctedStack = zeros(size(im_In));
Drift=zeros(size(imCropped,3),2);

mid = round(size(imCropped,3)/2);
imRef = stdfilt(mean(imCropped(:,:,mid+1:mid+driftPeriod),3));
h = waitbar(0,'Calculating drift');
for i=1:floor(size(imCropped,3)/driftPeriod)
    corrMatrix = normxcorr2(imRef,...
        stdfilt(mean(imCropped(:,:,(i-1)*driftPeriod+1:i*driftPeriod),3)));
    
    corrMatCenter=corrMatrix(...
    ceil(size(corrMatrix,1)/2)-10:ceil(size(corrMatrix,1)/2)+10,...
    ceil(size(corrMatrix,2)/2)-10:ceil(size(corrMatrix,2)/2)+10);

    [row,col,~,~,~] = Gauss.phasor(corrMatCenter);
%     
%     [X,Y] = meshgrid(1:size(corrMatCenter,2),1:size(corrMatCenter,1));
%     x0 = size(corrMatCenter,2)/2;
%     y0 = size(corrMatCenter,1)/2;
%     domain(:,:,1) = X;
%     domain(:,:,2) = Y;
%     % Gaussian fit
%     [gPar] = Gauss.MultipleFitting(corrMatCenter,x0,y0,domain,1,0);    
% 
%     DriftA((i-1)*driftPeriod+1:i*driftPeriod,1) = -(size(corrMatCenter,1)+1)/2+gPar(6);
%     DriftA((i-1)*driftPeriod+1:i*driftPeriod,2) = -(size(corrMatCenter,2)+1)/2+gPar(5);
%     
    Drift((i-1)*driftPeriod+1:i*driftPeriod,1) = row;
    Drift((i-1)*driftPeriod+1:i*driftPeriod,2) = col;
    
    
    waitbar(i/floor(size(imCropped,3)/driftPeriod),h,['Calculating drift ' num2str(i) '/' num2str(floor(size(imCropped,3)/driftPeriod))])
end
close(h)

h = waitbar(0,'Starting drift correction');
for i=1:size(imCropped,3)
    im2Correct = imresize(im_In(:,:,i),scalingFactor);
    [tmp] = PreProcess.correctImageDrift(im2Correct,...
        -round(Drift(i,:)*scalingFactor));
    
    correctedStack(:,:,i) = imresize(tmp,1/scalingFactor);
    waitbar(i/size(imCropped,3),h,['Drift correction ' num2str(i) '/' num2str(size(imCropped,3))])
end
close(h);

end

function [lb, ub] = getCropPar(im,corrSz)
    %find center of data
%     BW = imbinarize(uint16(im));
%     cent = regionprops(BW,'centroid');
%     tmp = [cent.Centroid];
%     coord(:,1) = tmp(1:2:end);
%     coord(:,2) = tmp(2:2:end);
%     xIdx = round(mean(coord(:,1)));
%     yIdx = round(mean(coord(:,2)));
%     
    maxProjY = max(im,[],2);
    maxProjX = max(im,[],1);
    
    maxProjY(maxProjY<1.2*min(maxProjY)) = 0;
    maxProjX(maxProjX<1.2*min(maxProjX)) = 0;
    
    [~,locY] = findpeaks(maxProjY);
    [~,locX] = findpeaks(maxProjX);
    xIdx = round(mean(locX));
    yIdx = round(mean(locY));
   
    %Crop the image to save computing time in the correlation
    lb = [xIdx - corrSz + 1, yIdx - corrSz + 1];
    ub = [xIdx + corrSz - 1, yIdx + corrSz - 1];
    %Check that the cropping occurs within the frame
    if lb(1)<1
        lb(1) = 1;
    end

    if lb(2)<1
        lb(2) = 1;
    end

    if ub(1)>size(im,1)
        ub(1) = size(im,1);
    end

    if ub(2)>size(im,2)
        ub(2) = size(im,2);
    end

end