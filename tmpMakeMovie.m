

filename = [myMovie.pathRes filesep 'clusterMovie'];
frameRate = 20;
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780  0.1840],[1 1 1]};

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
    imagesc(data(:,:,i))
    colormap('gray');
%     map = colormap('gray')*2^16;
%     
%     %cFrame = ind2rgb(uint16(data1(:,:,i)),'gray');    
%     cFrame = cat(3, data(:,:,i),data(:,:,i),data(:,:,i));
%     imagesc(cFrame)
    hold on
    
    frame2 =label2rgb(corrMask,'colorcube','k','shuffle');
    im = imagesc(frame2);
    im.AlphaData = 1-transparancy;
    
    axis image;
    
       
    %add scale  bar
    scaleBarPx = scaleBar/200*1000;
    x = size(cFrame,2)-scaleBarPx-(0.05*size(cFrame,2)):size(cFrame,2)-0.05*size(cFrame,2);
    y = ones(1,length(x))*size(cFrame,1)-0.05*size(cFrame,2);
    text(mean(x),mean(y)-0.05*size(cFrame,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
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