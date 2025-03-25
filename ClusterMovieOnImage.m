

filename = [myMovie.pathRes filesep 'CorrMapMovie'];
toPlot = 'corrMap';%clusters or corrMap
frameRate = 20;
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780  0.1840],[1 1 1]};
ROI = [50 125];
nFrames = 200;
Fig = figure;
scaleBar = 10;
transparancy = 1;
for i = 1:nFrames %loop over frames
    %plot original image
    if i>50
        if mod(i,10) ==0 %every 10 frames we upgrade transparency
           transparancy = transparancy-0.07; 
        end
    end
   
    im2Plot = data(ROI(1):ROI(1)+ROI(2),ROI(1):ROI(1)+ROI(2),i);
    imagesc(cat(2,im2Plot,fliplr(im2Plot)))
    colormap('gray');
    map = colormap('gray')*2^16;
    
%     cFrame = ind2rgb(data1(:,:,i),'gray');    
%    % cFrame = cat(3, data(:,:,i),data(:,:,i),data(:,:,i));
%     imagesc(cFrame)
    hold on
        
    frame2 =label2rgb(corrMask(ROI(1):ROI(1)+ROI(2),ROI(1):ROI(1)+ROI(2)),'colorcube','k','shuffle');
    corrMap = myMovie.corrRelation.corrMap(ROI(1):ROI(1)+ROI(2),ROI(1):ROI(1)+ROI(2));
    corrMap = corrMap-0.5;
    corrMap = corrMap./(1.1*max(corrMap));
    frame3 = ind2rgb(uint8(corrMap*255),colormap('gray'));

    im = imagesc(cat(2,frame2,fliplr(frame3)*255));
    im.AlphaData = 1-transparancy;
    axis image;
    %add scale  bar
    scaleBarPx = scaleBar/200*1000;
    x = size(im2Plot,2)-scaleBarPx-(0.05*size(im2Plot,2)):size(im2Plot,2)-0.05*size(im2Plot,2);
    y = ones(1,length(x))*size(im2Plot,1)-0.05*size(im2Plot,2);
    text(mean(x),mean(y)-0.05*size(im2Plot,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
    plot(x,y,'-w','LineWidth',3);
    
    set(gca,'visible','off');
    set(gcf,'color','w');
    drawnow;

    frame = getframe(Fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

    else

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

    end
    
    hold off
    clf;
end