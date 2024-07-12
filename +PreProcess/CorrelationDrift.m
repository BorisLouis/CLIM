function [correctedStack,Drift]=CorrelationDrift(im_In,scalingFactor,correlationInfo)

corrSz = correlationInfo.corrSz; 
driftPeriod = correlationInfo.driftPeriod;
maxDrift = correlationInfo.maxDrift;
frame1 = im_In(:,:,1);
%get cropping boundaries
[lb,ub] = getCropPar(frame1,corrSz);

%cropping
if or(isnan(lb),isnan(ub))
    imCropped = im_In;
else
    imCropped = im_In(lb(1): ub(1),lb(2):ub(2),:);

end

correctedStack = zeros(size(im_In));
Drift=zeros(size(imCropped,3),2);

mid = round(size(imCropped,3)/2);
imRef = mean(imCropped(:,:,1:1+driftPeriod),3);
fft_ref=fft2(imRef); 

h = waitbar(0,'Calculating drift');


for i=1:floor(size(imCropped,3)/driftPeriod)
    currentFrame = mean(imCropped(:,:,(i-1)*driftPeriod+1:i*driftPeriod),3);
    
%     corrMatrix = normxcorr2(imRef,currentFrame);
%     
     fft_frame=fft2(currentFrame);
%     prod=fft_ref.*conj(fft_frame);
%     cc=ifft2(prod);
%     [maxRow,maxCol]=find(fftshift(cc)==max(max(cc)));
    %testing
    usfac = 10;
    [output, Greg] = PreProcess.dftregistration(fft_ref,fft_frame,usfac);
    
    maxRow = output(3);
    maxCol = output(4);
    
     if i > 1
        if abs(maxRow)>2-Drift((i-1),1)
            maxRow = Drift((i-1),1);
        end
    
        if abs(maxCol)>2-Drift((i-1),2)
            maxCol = Drift((i-1),2);
        end
     end

    Drift((i-1)*driftPeriod+1:i*driftPeriod,1) = maxRow;
    Drift((i-1)*driftPeriod+1:i*driftPeriod,2) = maxCol;

%     
    
    waitbar(i/floor(size(imCropped,3)/driftPeriod),h,['Calculating drift ' num2str(i) '/' num2str(floor(size(imCropped,3)/driftPeriod))])
end
close(h)

h = waitbar(0,'Starting drift correction');
usfac2 = 5;
corrDrift(:,1) = smooth(Drift(:,1),50);
corrDrift(:,2) = smooth(Drift(:,2),50);
corrDrift = round(corrDrift*usfac2);
for i=1:size(imCropped,3)
    im2Correct = imresize(im_In(:,:,i),usfac2);
    [tmp] = PreProcess.correctImageDrift(im2Correct,...
        +corrDrift(i,:));
    
    correctedStack(:,:,i) = imresize(tmp,1/usfac2);
    waitbar(i/size(imCropped,3),h,['Drift correction ' num2str(i) '/' num2str(size(imCropped,3))])
end
%remove circ Drift.
minDrift = round(min(Drift));
maxDrift = round(max(Drift));

correctedStack([1:minDrift(1),end+minDrift(1)+1:end],:,:) = 0;
correctedStack(:,[1:minDrift(2),end+minDrift(2)+1:end],:) = 0;
correctedStack([1:maxDrift(1),end+maxDrift(1)+1:end],:,:) = 0;
correctedStack(:,[1:maxDrift(2),end+maxDrift(2)+1:end],:) = 0;




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
    lb = [yIdx - corrSz + 1, xIdx - corrSz + 1];
    ub = [yIdx + corrSz - 1, xIdx + corrSz - 1];
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