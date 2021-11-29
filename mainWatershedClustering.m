%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Watershed based clustering of intensity correlation in images           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Keep in mind that the deconvolution will affect the background, the
% best would be to make an ROI that removes the background to conteract
% this
%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\29 - Algorithm comparison\BiggerGrain';
file.ext  = '.spe';

info.runMethod  = 'load';
info.driftCorr = true;
info.ROI = true;

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.thresh = 0.6;%correlation threshold Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%
ROI = [96,96,64,64];
data1 = myMovie.loadFrames(frame2Process,ROI);


%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1);
%% Get pixels correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData,corrInfo);


%% Watershed clustering
figure
subplot(1,2,1)
imagesc(corrRelation.corrMap);
axis image


% watershed test
corrMask = watershed(1-corrRelation.corrMap);
subplot(1,2,2)
imagesc(corrMask);
axis image

figure
imagesc(imfuse(corrRelation.corrMap,corrMask))

%% Cluster evaluation
[clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,correctedData);

relData{1} = relNum1;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);


%% 
myMovie.plotContour(data1,corrMask);
%%
% % merge local maxima
% LocMax = imregionalmax(Ds); %find local maxima
% LocMax2 = imdilate(LocMax, strel('disk',10));%merge local maxima
% LocMax3 = bwmorph(LocMax2,'thin',8);
% LocMax3(~im) = 0;
% %    figure, imshowpair(LocMax3,LocMax2,'falsecolor')
% 
% %% perform watershed
% imD = -Ds;
% imD(~im) = Inf;
% imD(LocMax3) = min(imD(:));
% ws = watershed(imD);
% ws(~im) = 0;
% ws = bwareaopen(logical(ws),2);%remove pores smaller than 2px
% ws = bwlabel(logical(ws));