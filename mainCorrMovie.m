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

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.thresh = 0.6;%correlation threshold Pearson coefficient
corrInfo.neighbor = 4; %4 or 8 (8 means that diagonal are taken too)

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



%% get correlation mask from deconvolve data
corrInfo.thresh = 0.8; 
[corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);
[clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,correctedData);

relData{1} = relNum1;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);

%% Scanning threshold

thresh = 0.3:0.05:0.95;
allCorrMask = zeros(size(correctedData,1),size(correctedData,2),length(thresh));
threshold = zeros(length(thresh),1);
treatedArea = zeros(length(thresh),1);

for i = 1:length(thresh)
    corrInfo.thresh = thresh(i);
    [corrMask] = myMovie.getCorrelationMask(correctedData,corrInfo);
    
    %clean CorrMask
    a = regionprops(corrMask,'MajorAxisLength','MinorAxisLength');
    
    idx = find(and([a.MinorAxisLength]<3,[a.MajorAxisLength]./[a.MinorAxisLength]>2));
    cleanCorrMask = corrMask;
    cleanCorrMask(ismember(cleanCorrMask,idx)) = 0;
    
    [clustEval1,relNum1(i)] = corrAnalysis.evalClusters(corrMask,correctedData);
    
    [clustEval1Clean,relNum1Clean(i)] = corrAnalysis.evalClusters(cleanCorrMask,correctedData);
    threshold(i) = corrInfo.thresh;
    %calculate % of data treated
    binarizedMask = corrMask>0;

    treatedArea(i) = sum(binarizedMask(:))./numel(binarizedMask);
    
    allCorrMask(:,:,i) = cleanCorrMask;    

end




%% max number of cluster
relClusters = [relNum1.nClusters]./max([relNum1.nClusters]);
optMetric = treatedArea.*relClusters';
figure
plot(optMetric)

thresholdToUse = threshold(optMetric==max(optMetric));


%% clean the cluster based on statistical significance




%% Find best threshold
%calculate number relative number of cluster
relClusters = [relNum1.nClusters]./max([relNum1.nClusters]);


% Metric based on average correlation within cluster, intercluster
% correlation and number of clusters
corrMetric1= ([relNum1.meanCorr].*[relNum1.minCorr])./([relNum1.stdInterClusterCorr].*relClusters);

corrMetric2 = ([relNum1.meanCorr].*[relNum1.minCorr].*treatedArea')./([relNum1.stdInterClusterCorr].*relClusters);

corrMetric3 = ([relNum1.meanCorr].*[relNum1.minCorr].*treatedArea')./([relNum1.stdInterClusterCorr]);

figure
subplot(1,3,1)
plot(threshold,corrMetric1)
xlabel('Threshold')
ylabel('Goodness Metric')
box on
axis square
subplot(1,3,2)
plot(threshold,corrMetric2)
xlabel('Threshold')
ylabel('Goodness Metric')
box on
axis square
subplot(1,3,3)
plot(threshold,corrMetric3)
xlabel('Threshold')
ylabel('Goodness Metric')
box on
axis square

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

myMovie.plotContour(corrRelation.corrMap,cleanCorrMask);

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