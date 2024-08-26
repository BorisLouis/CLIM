%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for localization analysis on clusters                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

%% User input
file.path = 'D:\Documents\Unif\PostDoc\2024 - Data\08 - August\Blinking Perovskite\2024-08-07 - nonCS';
file.ext  = '';

pxSize = 0.2;%in um
info.runMethod = 'load';%load % load will try to load existing data from previous run
info.driftCorr = true; % true to correct for drift, false to not
deconvolve = true; %to deconvolve the correlated signal
backgroundThresh = 0.01; %0.1 is default, 0 is for no background removal
info.useThreshold = true;%false
info.doPlot = false;% default-false, do plot will generate a movie of the clustering
%procedure as it goes.
info.ROI = false; %this is to use ROI for the whole analysis
%[x y  w h]
%ROI = [];
ROI = [125 75 150 150];
%for intensity extraction
method = 'Mean'; %'Mean'
frame2Process = 1:5000; %number of frame to used for correlation analysis.

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);
myMovie.correctDrift;
    
filename = dir(file.path);
idx = contains({filename.name},'.mat');

load([filename(idx).folder filesep filename(idx).name]);




%% Loading frames  
data1 = myMovie.loadFrames(1:myMovie.raw.movInfo.maxFrame,ROI);
%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1,backgroundThresh,deconvolve);

%% show clusters
open([file.path filesep 'clusterID.fig'])


%% Run localization for every cluster:
% !!!!!!!!!!!!!!!!!!!! Time consuming!!!!!!!!!!!!!!!!!
%[allLoc] =myMovie.getAllClusterLocalization(corrOutput.corrMask,data);



%% Calculate localization on a specific cluster
idxToCluster = 1104;
[Loc,contour] = Localization.getClusterLocalization(correctedData,corrOutput.corrMask,idxToCluster);


figure
subplot(1,2,1)
scatter(Loc.x,Loc.y,10,'filled')
title('Localization and cluster boundary')
hold on
plot(contour(:,2),contour(:,1),'k','Linewidth',2)
axis square
axis ij
box on
subplot(1,2,2)
scatter(Loc.x,Loc.y,10,'filled')
title('Zoom on the Localization ')
axis image
box on
axis ij


pos = sqrt(Loc.x.^2 + Loc.y.^2);
figure
subplot(1,3,1)
scatter(pos,Loc.int,10,'filled')
axis square
box on
title('Intensity vs Position')

subplot(1,3,2)
scatter(pos,Loc.angle,10,'filled')
axis square
box on
title('Position vs angle')

subplot(1,3,3)
scatter(Loc.x,Loc.y,10,'filled')
axis square
box on
title('Intensity vs Angle')

%% test plot
myMovie.plotContour(data);%raw or clean depending on which we want to use
hold on
scatter(Loc.x,Loc.y,5,'filled')