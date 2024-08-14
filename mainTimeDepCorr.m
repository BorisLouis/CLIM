%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for temporal anaylsis of correlation of fluctuation           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\08 - August\Blinking Perovskite\Spot1_dryObj glass side';
file.ext  = '';
info.runMethod = 'run';
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.ROI = true; %this is to use ROI for the whole analysis
%[x y  w h]
ROI = [75 75 100 100];
%ROI = [5 71 230 120];

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
%% Loading frames  
data1 = myMovie.loadFrames(1:myMovie.raw.movInfo.maxFrame,ROI);

%% save data 
%myMovie.saveMovie(data1,50);
%dataStorage.nBTiff('driftCorrected.tif',data1,16);

%% Deconvolution

[correctedData] =  myMovie.deconvolve(data1,backgroundThresh,deconvolve);


%% time window correlation map analysis

nFrames = 100;%number of frame to used for correlation analysis.
lagTime  = 20; %this means that startign points will be 500 frames appart

procFrame = length(correctedData) - nFrames;
idx = 1:lagTime:procFrame+1;

h = waitbar(1/length(idx),sprintf('%d/%d time point',1,length(idx)));
tic
allCorrMap = zeros(size(correctedData,1),size(correctedData,2),length(idx));
for i = 1:length(idx)
    
    frame2Process = idx(i):idx(i)+nFrames-1;
    data2Use = correctedData(:,:,frame2Process);
    
    [corrRelation] = myMovie.getPxCorrelation(data2Use);
    allCorrMap(:,:,i) = corrRelation.corrMap;
    
    waitbar(i/length(idx),h, sprintf('%d/%d time point',i,length(idx)));
end
toc
close(h)

allDiffMap = diff(allCorrMap,1,3);
%% make GIF
data = allCorrMap;
filename = [myMovie.pathRes filesep 'tempDepCorr.gif'];
frameRate = 5;
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780  0.1840],[1 1 1]};

Fig = figure;
scaleBar = 10;

for i = 1:length(idx)
 	imagesc(data(:,:,i))
    colormap('jet');
    caxis([0.5 1])
    %caxis([0 0.1])
    colorbar
    axis image;
    hold on
    %add scale  bar
    scaleBarPx = scaleBar/200*1000;
    x = size(data,2)-scaleBarPx-(0.05*size(data,2)):size(data,2)-0.05*size(data,2);
    y = ones(1,length(x))*size(data,1)-0.05*size(data,2);
    text(mean(x),mean(y)-0.05*size(data,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
    
    plot(x,y,'-w','LineWidth',3);
    text(15,5,[num2str(idx(i)) ' fr'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
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

%% DiffMap
allDiffMap = diff(allCorrMap,1,3); % frame by frame difference
%data = allDiffMap;
data = allCorrMap-squeeze(mean(allCorrMap,3));
filename = [myMovie.pathRes filesep 'diffTempDepCorr.gif'];
frameRate = 5;
color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780  0.1840],[1 1 1]};

Fig = figure;
scaleBar = 10;

for i = 1:length(idx)
 	imagesc(data(:,:,i))
    colormap('jet');
    caxis([0 0.2])
    %caxis([0 0.1])
    colorbar
    axis image;
    hold on
    %add scale  bar
    scaleBarPx = scaleBar/200*1000;
    x = size(data,2)-scaleBarPx-(0.05*size(data,2)):size(data,2)-0.05*size(data,2);
    y = ones(1,length(x))*size(data,1)-0.05*size(data,2);
    text(mean(x),mean(y)-0.05*size(data,1),[num2str(scaleBar) ' µm'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
    
    plot(x,y,'-w','LineWidth',3);
    text(15,5,[num2str(idx(i)) ' fr'],'HorizontalAlignment','center','Color','white','fontWeight','bold','fontSize',14);
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


