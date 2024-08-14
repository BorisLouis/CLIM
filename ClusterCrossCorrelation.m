%% correlation analysis


%figure(1)
highTimeCorrList = cell(size(traces,2),1);
bestR = zeros(size(traces,2),length(traces(1).trace));
for i = 1:size(traces,2)
   idx = 1:size(traces,2);
   currTrace = traces(i).trace;
   
   for j = 1:size(traces,2)
       if j == i
           r(:,j) = zeros(6001,1);
       else
           trace2Test = traces(j).trace;

%            r(:,j)  = xcov(currTrace,trace2Test,'normalized');
%            r1(:,j) = normxcorr2(currTrace,trace2Test);
           r(:,j) = crosscorr(currTrace,trace2Test,'NumLags',3000,'NumSTD',3);
           
%            plot(r(:,j),'DisplayName','r')
%            axis square
%            clf
       end

   end
   [val,id] = max(r,[],1);
   
   idx2Best = find(val>0.5);
   
   highTimeCorrList{i,1} = idx2Best;
   highTimeCorrList{i,2} = id(val>0.5) - 3001;
   
end

%% store the data

for i =1:size(highTimeCorrList,1)
    traces(i).correlatedClust = highTimeCorrList{i,1};
    traces(i).lag = highTimeCorrList{i,2}; 
end


%% plot them
idx = 1;
tmpIm = zeros(size(silMap));

correlatedClust = traces(idx).correlatedClust;
lags = traces(idx).lag;

tmpIm(traces(idx).clustPos) = 1;
for i = 1:length(correlatedClust)
   
    currIdx = correlatedClust(i);
    
    tmpIm(traces(currIdx).clustPos) = traces(currIdx).nPx;
        
    
end


figure
imagesc(tmpIm)









