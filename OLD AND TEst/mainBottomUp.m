%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Keep in mind that the deconvolution will affect the background, the
% best would be to make an ROI that removes the background to conteract
% this
%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\11 - November\30 - Big algorithm evaluation\Bigger Grain';
file.ext  = '.spe';

info.runMethod  = 'run';
info.driftCorr = true;
info.ROI = true;

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.thresh = 0.6;%correlation threshold Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%
ROI = [26,100,64,64];
data1 = myMovie.loadFrames(frame2Process,ROI);


%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1);
%% Get pixels correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData,corrInfo);

%% calculate correlation matrix

[n,p] = ind2sub(size(correctedData),corrRelation.indPx);
pxIntList = zeros(length(n),size(correctedData,3),'single');
for i =1:length(n)

    pxIntList(i,:) = single(correctedData(n(i),p(i),:));

end

corrMat = corrcoef(pxIntList');


%% 
meanPx = corrRelation.meanPx;
inds = corrRelation.indPx;
[corrMask] = corrAnalysis.corrClusteringBottomUp(meanPx,inds,corrMat,correctedData);

%%
figure
subplot(1,2,1)
imagesc(corrRelation.corrMap)
axis image
subplot(1,2,2)
imagesc(corrMask)
axis image


