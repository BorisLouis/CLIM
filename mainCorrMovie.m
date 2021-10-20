%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Main code for analysis of intensity fluctuation correlation in images   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input
file.path = 'D:\Documents\Unif\PhD\2021-Data\10 - October\20 - Film Blinking\Big Grain Ambiant';
file.ext  = '.spe';

info.runMethod  = 'run';
info.driftCorr = true;
info.ROI = false;

frame2Process = 1:6000;
corrInfo.r = 1; %radius for checking neighbor
corrInfo.thresh = 0.4;%correlation threshold (smaller is more correlation)==> 0.6 == 0.4 Pearson coefficient

%% Loading data
myMovie = Core.CorrClusterMovie(file,info);

myMovie.correctDrift;

    
%%

data1 = myMovie.loadFrames(frame2Process);

meanData1 = smooth(squeeze(mean(mean(data1,1),2)));
% 
testData = [squeeze(data1(128,128,:));squeeze(data1(128,128,:))];

figure
subplot(1,3,1)
plot(testData)
subplot(1,3,2)
plot(meanData1)
subplot(1,3,3)
[a,r] = deconv(testData,[meanData1;meanData1]);

plot(r+mean(testData))


%test deconvolution on the whole image
%Need to improve speed here
meanData1 = smooth(squeeze(mean(mean(data1,1),2)));
correctedData = data1;

for i =1:size(data1,1)
    for j=1:size(data1,2)
        currentData = squeeze(data1(i,j,:));
        [a,r] = deconv(currentData,meanData1);
        cleanData = r+mean(currentData);
        correctedData(i,j,:) = cleanData;
    end
end

data2Use = correctedData;

%% Get pixels correlation
[listCorrPx,inds] = myMovie.getPxCorrelation(data2Use,corrInfo);
%

[corrMask,cleanedCorrMask] = myMovie.getCorrelationMask(data2Use,corrInfo);
%
%%
%compare the two clusters
%[~,relNum2] = compare2Cluster(corrMask,cleanedCorrMask,data,'V1');
[clustEval1,relNum1] = corrAnalysis.evalClusters(corrMask,data1);

relData{1} = relNum1;
label{1}   = ['Method' '-pseudoClust'];
corrAnalysis.compareClusters(relData,label);

%% Plotting
myMovie.plotContour(data1,'raw');%raw or clean depending on which we want to use


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