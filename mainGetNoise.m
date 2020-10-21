%%
file.path = 'D:\Documents\Unif\PhD\2020-Data\10 - Oct\Sudipta\Noise - Laser Fluctuations';
file.ext  = '.spe';
info.runMethod  = 'load';

myMovie = Core.Movie(file,info);

%% load all frames

nFrames = myMovie.raw.maxFrame;
data = zeros(myMovie.raw.movInfo.Length,myMovie.raw.movInfo.Width,nFrames);
for i = 1:nFrames
    
    data(:,:,i) = myMovie.getFrame(i);
        
end


%% Get Background
figure
imagesc(data(:,:,1));
h=drawrectangle();
posBkg = round(h.Position);
close
%% extract distribution from bkg

bkg = data(posBkg(2):posBkg(2)+posBkg(4)-1,posBkg(1):posBkg(1)+posBkg(3)-1,:);

[N,edges] = histcounts(bkg(:),1000);

edges = edges(1:end-1) + mean(diff(edges));

figure
bar(edges,N)

% get cumlative distribution
bkgData = Dist.getCDF(bkg(:));
bkgData.center = edges;
bkgData.weight = N;
bkgData.Avg = mean(bkg(:));

figure
plot(bkgData.x,bkgData.y)

% % %test Sampling ==========> OK !!!
% nSamples = 20000;
% %testData = zeros(1,nSamples);
% 
% testData = datasample(edges,nSamples,'weight',N);
% % for i = 1:nSamples
% % %     n = rand(1);
% % %     
% % %     [~,idx] = min(abs(bkgCDF.y - n));
% % %     testData(i) = bkgCDF.x(idx(1));
% %     
% %     
% % end
% 
% [testN,testEdges] = histcounts(testData);
% testEdges = testEdges(1:end-1) + mean(diff(testEdges));
% figure
% bar(edges,N/max(N))
% hold on
% bar(testEdges,testN/max(testN));

%% Get Laser fluctuations
figure
imagesc(data(:,:,1));
h=drawrectangle();
posLas = round(h.Position);
close

%% Extract distibution from Laser
las = data(posLas(2):posLas(2)+posLas(4)-1,posLas(1):posLas(1)+posLas(3)-1,:);

[N,edges] = histcounts(las(:),1000);

edges = edges(1:end-1) + mean(diff(edges));

figure
bar(edges,N)

% get cumlative distribution
lasData = Dist.getCDF(las(:));
lasData.Avg = mean(las(:));
lasData.center = edges;
lasData.weight = N;
figure
plot(lasData.x,lasData.y)


%% Saving

fileName = [file.path filesep 'noiseDist.mat'];
save(fileName,'bkgData');

fileName = [file.path filesep 'laserDist.mat'];
save(fileName,'lasData');




