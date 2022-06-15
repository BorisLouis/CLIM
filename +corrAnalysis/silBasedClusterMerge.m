function [newCorrMask] = silBasedClusterMerge(corrM,data,clustTrace,meanCorr)
    
    traces = [clustTrace.trace];

    clusterCorr = corrcoef(traces);
    
    label = corrM(corrM>0);
    %px2Move = table(zeros(size(label)),zeros(size(label)),zeros(size(label)),'VariableNames',{'currClust','Sil','bestClust'});
    lab = unique(label);
    
    for i = lab
        
        
        
    end

    end



