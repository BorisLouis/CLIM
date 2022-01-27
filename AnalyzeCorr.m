clc 
close all

%% Analyze data and make plot
pxSize = 0.2;
pxArea = pxSize^2;


%% get folder to save plots
[folder,file,ext] = fileparts(corrOutput.path);




%% Make map/mask plot

f1 = figure('Position',[250 100 1200 800]);
subplot(1,3,1)
imagesc(corrOutput.Image)
axis image
colormap(gca,'gray')
title('PL image');
subplot(1,3,2)
imagesc(corrOutput.corrMap)
axis image
title('Correlation map');
colormap(gca,'parula')
cb = colorbar;
cb.Position = [0.63 0.355 0.012 0.325];
subplot(1,3,3)
RGBIM = label2rgb(corrOutput.corrMask,normCol,'k','shuffle');
%RGBIM = label2rgb(corrOutput.corrMask,'parula','k','shuffle');
imagesc(RGBIM);
%imagesc(corrOutput.corrMask);
axis image
colormap(gca,'colorcube')
title('Cluster map')

fileName = [folder filesep 'CorrerelationClust.fig'];
saveas(f1,fileName)

fileName = [folder filesep 'CorrerelationClust.svg'];
saveas(f1,fileName,'svg')


%% Contour
f2 = figure('Position',[250 100 1200 800]);
subplot(1,2,1)
hold on
plotContour(f2,corrOutput.Image,corrOutput.corrMask);
colormap(gca,'gray')
title('Mask on raw Data')
subplot(1,2,2)
hold on
plotContour(f2,corrOutput.corrMap,corrOutput.corrMask);
title('Mask on correlation Map')

fileName = [folder filesep 'Contour.fig'];
saveas(f2,fileName)

fileName = [folder filesep 'Contour.svg'];
saveas(f2,fileName,'svg')
%% Different maps

silColors = colormap('jet');
silColors(1,:) = [0,0,0];

f3 = figure('Position',[250 100 1200 800]);

subplot(1,3,1)
hold on
imagesc(corrOutput.corrMap)
colormap(gca,'parula')
cb = colorbar;
cb.Position = [0.35 0.355 0.012 0.325];
title('Correlation Map')
axis image
axis ij
subplot(1,3,2)
hold on
imagesc(corrOutput.corrClustMap)
colormap(gca,'jet')

cb = colorbar;
cb.Position = [0.63 0.355 0.012 0.325];
title('Correlation to cluster');
axis image
axis ij
subplot(1,3,3)

hold on

imagesc(corrOutput.silMap)
cb = colorbar;
cb.Position = [0.91 0.355 0.012 0.325];
caxis([0 max(corrOutput.silMap(:))])
%colormap(gca,silColors)
colormap(gca,'jet')
title('Silhouette Map')
axis image
axis ij


fileName = [folder filesep 'corrMap.fig'];
saveas(f3,fileName)

fileName = [folder filesep 'corrMap.svg'];
saveas(f3,fileName,'svg')

%% Population analysis Pt1 - Size dependence

f4 = figure('Position',[250 100 1200 800]);
subplot(2,2,1)
scatter([corrOutput.results.nPx].*pxArea,[corrOutput.results.meanSil],20,'filled')
axis square
box on
xlim([0 5])
ylim([0 1])
xlabel('Cluster Size (\mum)')
ylabel('Silhouette')
title('Silhouette vs Size')

subplot(2,2,2)
scatter([corrOutput.results.nPx].*pxArea,[corrOutput.results.Intensity],20,'filled')
axis square
xlim([0 5])
box on
xlabel('Cluster Size (\mum)')
ylabel('Avg. Px. Intensity (A.U.)')
title('Intensity vs Size')


subplot(2,2,3)
scatter([corrOutput.results.nPx].*pxArea,[corrOutput.results.meanCorr],20,'filled')
axis square
xlim([0 5])
ylim([0 1])
box on
xlabel('Cluster Size (\mum)')
ylabel('meanCorr')
title('Correlation vs Size')

subplot(2,2,4)
scatter([corrOutput.results.nPx].*pxArea,[corrOutput.results.maxInterClustCorr],20,'filled')
axis square
xlim([0 5])
ylim([0 1])
box on
xlabel('Cluster Size (\mum)')
ylabel('max Inter-cluster corr')
title('maxCluster Intercorr vs Size')

fileName = [folder filesep 'SizeDependence.fig'];
saveas(f4,fileName)

fileName = [folder filesep 'SizeDependence.svg'];
saveas(f4,fileName,'svg')

%% Intensity Dependence

f5 = figure('Position',[250 100 1200 800]);

subplot(1,3,1)
scatter([corrOutput.results.Intensity],[corrOutput.results.meanSil],20,'filled')
axis square
box on

xlabel('Avg. Px. Intensity (A.U.)')
ylabel('Silhouette')
title('Intensity vs silhouette')
ylim([0 1])

subplot(1,3,2)
scatter([corrOutput.results.Intensity],[corrOutput.results.meanCorr],20,'filled')
axis square
box on
xlabel('Avg. Px. Intensity (A.U.)')
ylabel('Mean Correlation')
title('Correlation vs Intensity')
ylim([0 1])

subplot(1,3,3)
scatter([corrOutput.results.Intensity],[corrOutput.results.maxInterClustCorr],20,'filled')
axis square
box on
xlabel('Avg. Px. Intensity (A.U.)')
ylabel('max Inter-cluster corr')
title('maxCluster Intercorr vs Intensity')
ylim([0 1])

fileName = [folder filesep 'IntensityDependence.fig'];
saveas(f5,fileName)

fileName = [folder filesep 'IntensityDependence.svg'];
saveas(f5,fileName,'svg')

%%
function plotContour(f,IM,corrM)

    n = 1;
    for i = unique(corrM(:))'
        currMask = corrM == i;
        if not(all(currMask(:)==0))
            contour{n} = Plotting.bwperimtrace(corrM == i,[1 size(corrM,2)],[1 size(corrM,1)]);
            n=n+1;
        end
    end

    f;
    hold on
    imagesc(IM)
    for i = 1:length(contour)


        %contour = bwboundaries(corrMaskCopy);

        plot(contour{i}{1}(:,1),contour{i}{1}(:,2),'w','LineWidth',2)

    end
    axis ij
    xlim([1 size(IM,2)])
    ylim([1 size(IM,1)]);
    axis image
end
