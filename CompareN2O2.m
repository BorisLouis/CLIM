
%% Compare N2 and O2 Data


%loop through oxygen cluster
n = 1;
for i = 1:length(corrOutputO2.results)
    currentClustPos = corrOutputO2.results(i).clustPos;
    for j=1:length(corrOutputN2.results)
       currN2Clust = corrOutputN2.results(j).clustPos; 
       
       intersection = intersect(currentClustPos,currN2Clust);
       
       if ~isempty(intersection)
          intersec = length(intersection)/length(currentClustPos);
          
          if intersec>0.75
              
              compRes(n).traceO2 = corrOutputO2.results(i).trace;
              compRes(n).traceN2 = corrOutputN2.results(j).trace;
              
              compRes(n).nPxO2 = corrOutputO2.results(i).nPx;
              compRes(n).nPxN2 = corrOutputN2.results(j).nPx;
              
              compRes(n).AmpO2 = corrOutputO2.results(i).Amp;
              compRes(n).AmpN2 = corrOutputN2.results(j).Amp;
              
              compRes(n).SilO2 = corrOutputO2.results(i).meanSil;
              compRes(n).SilN2 = corrOutputN2.results(j).meanSil;
              
              compRes(n).idxO2 = i;
              compRes(n).idxN2 = j;
             
              compRes(n).overlap = intersec;
              n = n+1;
              disp('overlap found');
          end
           
       end
    end
    
end

%% Plot comparison

figure('Position',[200,200,1200,800])

subplot(2,4,1)
idx = find([compRes.overlap]==1);
rIdx = randi([1 length(idx)],1);
plot(compRes(rIdx).traceO2)
hold on
plot(compRes(rIdx).traceN2);
legend({'O2','N2'})
xlabel('Time (frames)')
ylabel('Intensity (A.U.)');
axis square
box on

subplot(2,4,2)
edges = [0:0.05:1];
NO2 = histcounts([compRes.AmpO2],edges);
NN2 = histcounts([compRes.AmpN2],edges);
hold on
bar(edges(1:end-1),NO2)
bar(edges(1:end-1),NN2);
legend({'O2','N2'})
xlabel('Amplitude')
ylabel('Occurence')
axis square
box on

subplot(2,4,3)
edges = [5:5:200];
NO2 = histcounts([compRes.nPxO2],edges);
NN2 = histcounts([compRes.nPxN2],edges);
hold on
plot(edges(1:end-1),NO2)
plot(edges(1:end-1),NN2);
legend({'O2','N2'})
xlabel('Area (px)')
ylabel('Occurences')
axis square
box on

subplot(2,4,4)
edges = [0:0.05:1];
NO2 = histcounts([compRes.SilO2],edges);
NN2 = histcounts([compRes.SilN2],edges);
hold on
plot(edges(1:end-1),NO2)
plot(edges(1:end-1),NN2);
legend({'O2','N2'})
xlabel('Silhouette')
ylabel('occurences')
axis square
box on


subplot(2,4,5)
idx = find([compRes.overlap]==1);
rIdx = randi([1 length(idx)],1);
plot(compRes(rIdx).traceO2)
hold on
plot(compRes(rIdx).traceN2);
legend({'O2','N2'})
xlabel('Time (frames)')
ylabel('Intensity (A.U.)');
axis square
box on

subplot(2,4,6)
scatter([compRes.AmpO2],[compRes.AmpN2],10,'filled')
xlabel('Amplitude O2')
ylabel('Amplitude N2')
axis image
xlim([0.1 0.5])
ylim([0.1 0.5])

box on

subplot(2,4,7)
scatter([compRes.nPxO2]*0.04,[compRes.nPxN2]*0.04,10,'filled')
xlabel('Area O2 (\mum)')
ylabel('Area N2 (\mum)')
axis image
xlim([0 10])
ylim([0 10])
box on

subplot(2,4,8)
scatter([compRes.SilO2],[compRes.SilN2],10,'filled')
xlabel('Silhouette O2')
ylabel('Silhouette N2')
axis image
xlim([0.1 0.5])
ylim([0.1 0.5])

box on
