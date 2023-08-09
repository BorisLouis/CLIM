clearvars -except corrOutput

%% User input
idx2Plot = 10;


%% Get distance and correlation between clusters

%extract data from loaded data
traces = zeros(length(corrOutput.results(1).trace),length(corrOutput.results));
pos    = zeros(length(corrOutput.results),2);
idx    = zeros(length(corrOutput.results),1);
for i = 1:length(corrOutput.results)
    
    traces(:,i) = corrOutput.results(i).trace;
    clustPos = corrOutput.results(i).clustPos;
    [row,col] = ind2sub(size(corrOutput.corrMap),clustPos);
    pos(i,:) = [mean(row), mean(col)]; % [row col] = [y x]
    idx(i,:) = corrOutput.results(i).clustIdx;
     
    
end

% calculate the correlation matrix and distance matrix

corrMat = corrcoef(traces);

distMat = squareform(pdist(pos));

%% Plot decay example
idx2Plot = 10;
figure
scatter(distMat(idx2Plot,:),corrMat(idx2Plot,:),10,'filled')
axis square
box on

%% Image of the decay for single cluster
idx2Plot = 2;
bwimage = zeros(size(corrOutput.corrMap));
for i = 1:size(corrMat,1)
   bwimage(corrOutput.results(i).clustPos) = corrMat(idx2Plot,i);
   bwimage(bwimage==1) = max(corrMat(idx2Plot,corrMat(idx2Plot,:)<1))+0.1*max(corrMat(idx2Plot,corrMat(idx2Plot,:)<1));
  
end

figure
imagesc(bwimage)
colormap('jet')
caxis([0.05 max(bwimage(:))])
axis image
colorbar
%% fitting exponential
Results = table(zeros(size(distMat,1),1),zeros(size(distMat,1),1),zeros(size(distMat,1),1),...
    'Variablenames',{'ClustIdx','tauCorr','AmpCorr'});
allFit = zeros(size(corrMat));
sortedDMat = zeros(size(corrMat));
sortedCMat = zeros(size(corrMat));

% modfun = @(b,x) b(2).*exp(-(x/b(1))) +b(3);
% F = x(2) .* exp(-(data/x(1))+x(3));
for i = 1:size(distMat,1)
    
    distance = distMat(i,:);
    correlation = corrMat(i,:);
    
    [distance,idx] = sort(distance);
    correlation = correlation(idx);
    
    sortedDMat(i,:) = distance;
    sortedCMat(i,:) = correlation;
    
    [FitPar,fit]= SimpleFitting.expFit(correlation,distance);
%     tauGuess = prctile(distance,20);
%     AGuess = max(correlation)-min(correlation);
%     bkgGuess = min(correlation);
%     beta0 = [tauGuess,AGuess,bkgGuess];
%     [FitPar] = nlinfit(distance,correlation,SimpleFitting.exponential,beta0);
    allFit(i,:) = fit;
    Results.ClustIdx(i) = i;
    Results.tauCorr(i) = FitPar(1);
    Results.AmpCorr(i) = FitPar(2);
   
end


%% Plot the fit
idx = randi(size(corrMat,1),9);

figure
hold on
for i = 1:length(idx)
    id = idx(i);
    subplot(3,3,i)
    scatter(sortedDMat(id,:),sortedCMat(id,:),5,'filled')
    axis square
    box on
    hold on 
    plot(sortedDMat(id,:),allFit(id,:),'r','Linewidth',2)
    
end


%% image of the decay fit

decayImage = zeros(size(corrOutput.corrMap));
ampImage   = zeros(size(corrOutput.corrMap));
for i = 1:size(corrMat,1)
   decayImage(corrOutput.results(i).clustPos) = Results.tauCorr(i);
   ampImage(corrOutput.results(i).clustPos)   = Results.AmpCorr(i);
  
end

figure
subplot(1,2,1)

imagesc(decayImage)
colormap('jet')
caxis([0 30])
axis image
title('Correlation decay lifetime')
colorbar

subplot(1,2,2)
title('Correlation decay amplitude')
imagesc(ampImage)
colormap('jet')
axis image
colorbar
title('Correlation decay amplitude')







