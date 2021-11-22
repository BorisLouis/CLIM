function compareClusters(relData,labels)
    assert(iscell(relData),'data is expected to be a cell containing struct')
    assert(isstruct(relData{1}),'data is expected to be a cell containing struct');
    assert(iscell(labels),'label is expected to be a cell containing a string')
    assert(ischar(labels{1}),'label is expected to be a cell containing a string');
        
    nFields = numel(fieldnames(relData{1}));
    data2Plot = zeros(5,length(relData));
    
    figure
    hold on
    for i = 1 : length(relData)
        currData = relData{i};
        
        data2Plot(1,i) = currData.meanCorr;
        data2Plot(2,i) = currData.medCorr;
        data2Plot(3,i) = currData.minCorr;
        data2Plot(4,i) = currData.maxCorr;
        data2Plot(5,i) = currData.meanInterClusterCorr;
        data2Plot(6,i) = currData.stdInterClusterCorr;
        data2Plot(7,i) = currData.nClusters;
    
    end
    data2Plot(7,:) = data2Plot(7,:)./max(data2Plot(7,:)); 
    X = categorical({'mean','median','min','max','interCorrM','interCorrStd','nClusters'});
    hb = bar(X',data2Plot,'grouped');
    axis square
    box on
    legend(labels)
    ylabel('Metric')

end