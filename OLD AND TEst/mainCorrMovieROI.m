%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Keep in mind that the deconvolution will affect the background, the
% best would be to make an ROI that removes the background to conteract
% this
%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\12 - December\12 - Code improve attempt\BigGrain';
file.ext  = '.spe';

info.runMethod  = 'run';
info.driftCorr = true;
info.ROI = true;
ROI = [26,100,64,64];

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.neighbor = 8; %4 or 8 (8 means that diagonal are taken too)

minCorr = 0.4;%Minimum correlation we want to have
stepCorr = 0.05; %Correlation difference between different tested threshold
maxCorr = 0.9;%maximum correlation to be tested, higher than 0.9 makes little sense

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%
data1 = myMovie.loadFrames(frame2Process,ROI);


%% Deconvolution
[correctedData] =  myMovie.deconvolve(data1);

%% Get pixels correlation
[corrRelation] = myMovie.getPxCorrelation(correctedData,corrInfo);

%% Scanning threshold

thresh = minCorr:stepCorr:maxCorr;

allCorrMask = zeros(size(correctedData,1),size(correctedData,2),length(thresh));
threshold = zeros(length(thresh),1);
treatedArea = zeros(length(thresh),1);

rearrangeData = reshape(correctedData,[size(correctedData,1)*size(correctedData,2),size(correctedData,3)]);
for i = 1:length(thresh)
    corrInfo.thresh = thresh(i);
    [corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);
    
    label = reshape(corrMask,[size(corrMask,1)*size(corrMask,2),1]);
    
    %Calculate silhouette of clusters
    tmpRearrangeData = rearrangeData;
    silDist = silhouette(tmpRearrangeData,label,'Correlation');
    silLabel = label;
    sDist = silDist;
    sDist(label==0) = [];
    label(label==0) = [];
    sil(i) = mean(sDist);
       
    %[clustEval1,relNum1(i)] = corrAnalysis.evalClusters(corrMask,correctedData);
    threshold(i) = corrInfo.thresh;
    
    %calculate % of data treated   
    treatedArea(i) = sum(corrMask(:)>0)./numel(corrMask);
   
    allCorrMask(:,:,i) = corrMask;
   
end


%% find optimal threshold

optMetric = sil(:).*treatedArea(:);
figure
hold on
plot(threshold,optMetric)
xlabel('Correlation threshold')
ylabel('Mean Silhouette coefficient');
axis square
box on

guess.sig = 0.3;
guess.mu = threshold(optMetric==max(optMetric));
[FitPar,fit] = Gauss.gauss1D(optMetric,threshold,guess);

plot(threshold,fit);

threshold2Use = FitPar(2);

%% Calculate final corrMask
corrInfo.thresh = threshold2Use;
[corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);

label = reshape(corrMask,[size(corrMask,1)*size(corrMask,2),1]);

%Calculate silhouette of clusters
tmpRearrangeData = rearrangeData;
silDist = silhouette(tmpRearrangeData,label,'Correlation');
silLabel = label;
sDist = silDist;
sDist(label==0) = [];
label(label==0) = [];
bestSil = mean(sDist);



%%
%compare the two clusters
%[~,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');
[bestClustEval1,bestRelNum] = corrAnalysis.evalClusters(corrMask,correctedData);

relData{1} = bestRelNum;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);

%% Plotting
myMovie.plotContour(data1);%raw or clean depending on which we want to use

%% 

myMovie.plotContour(corrRelation.corrMap,corrMask);

%% Plot traces
myMovie.plotClusterTraces(data1,4);


%% Extract intensity traces 
data = myMovie.loadFrames(1:36000);

[traces] = myMovie.getAllTraces(data);


%% Get interCluster traces
% allTraces = zeros(size(traces,2),length(traces(1,1).trace));
% for i = 1:size(traces,2)
%     allTraces(i,:) = traces(1,i).trace;
%     
% end
% allTraces = allTraces';
% corrCoeff = corrcoef(allTraces);
% 
% figure
% imagesc(corrCoeff)
% colormap('jet')
% 
% U = triu(corrCoeff);
% U = nonzeros(U);
% tmpData = U(U<1);
% figure
% histogram(tmpData);
% 
% 
% % Let us try to get anti-correlated Crystal
% [val,idx] = mink(corrCoeff(:),20);
% %remove duplicate
% val = val(1:2:end);
% 
% 
% for i=1:length(val)
%     currVal = val(i);
%     idx = find(corrCoeff==currVal);
%     [row,col] = ind2sub(size(corrCoeff),idx(2));
% 
%     figure
%     hold on
%     plot(allTraces(:,row))
%     plot(allTraces(:,col))
%     test = corrcoef(allTraces(:,row),allTraces(:,col));
%     disp(test);
% end

%% get image from corrmask
color= 'colorcube';
[corrMaskIM] = myMovie.getImageFromMask(corrMask,color);

%% Calculate localization on cluster to get the position of the defect
%[allLoc] =myMovie.getAllClusterLocalization(data);

idx =112;
[Localization] = myMovie.getClusterLocalization(data1,idx);

pos = sqrt(Localization.x.^2 + Localization.y.^2);
figure
subplot(1,3,1)
scatter(pos,Localization.int,10,'filled')
axis square
box on
title('Intensity vs Position')

subplot(1,3,2)
scatter(pos,Localization.angle,10,'filled')
axis square
box on
title('Position vs angle')

subplot(1,3,3)
scatter(Localization.x,Localization.y,10,'filled')
axis square
box on
title('Intensity vs Angle')

%% test plot
myMovie.plotContour(data1,'raw');%raw or clean depending on which we want to use
hold on
scatter(Localization.x,Localization.y,5,'filled')

%%

function [relNum1,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,method)
    %Evaluate clusters individually
    [clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,data);
    [clustEval2,relNum2] = corrAnalysis.evalClusters(cleanedCorrMask,data);


    %compare cluster
    relData{1} = relNum1;
    relData{2} = relNum2;
    label{1}   = ['Method' method];
    label{2}   = ['Method',method,'- cleaned'];

    corrAnalysis.compareClusters(relData,label);


end