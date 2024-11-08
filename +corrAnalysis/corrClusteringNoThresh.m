function [corrMask,frames] = corrClusteringNoThresh(corrRelation,data,thresh,doPlot)
    %The idea of the modifications here is to:
    %1) Use listVal to and threshold to add to cluster(this will allow to
    %test various threshold without having to re-run the first steps)
    %2) Add multiple indices to the cluster at once instead of 1 by 1 to
    %reduce the processing time.
    listCorrPx = corrRelation.listPx;
    inds   = corrRelation.indPx;
    meanPx = corrRelation.meanPx;
    listVal = corrRelation.listVal;
    indsCopy = inds;
    %function actually does the clustering.
    safeCount = 2*size(data,1)*size(data,2);
    %initialize the variable
    count = 1;
    group = 1;
    corrMask = zeros(size(data,1),size(data,2));
    treatedIdx = NaN(length(inds),1);
    indsCopy = inds;
    clusters{group} = [];
    %as long as the list of pixel that are correlated is not empty
    %we keep going
    tic
    while ~isempty(listCorrPx) 
        
        %get Index of most correlated pixel m
        [~,idx] = max(meanPx);
        %take the index of the first pixel to be treated
        currIndex = inds(idx);
        %add to a new list the pixel that are correlated with the
        %currently treated pixel
        currList = listCorrPx{idx};
        currVal  = listVal{idx};
        
        currList(currVal<thresh) = [];
       
        tmpList = [currIndex; currList];
        
        %add the pixel to the cluster
        %corrMask(tmpList) = group;
        if isempty(currList)
        
        else
            if ~isempty(clusters{1})
                %check if there is an overlap between the "to-be-added"
                %pixels and an existing cluster
                overlap = zeros(size(clusters));
                for i = 1: length(clusters)
                    overlap(i) = sum(ismember(tmpList,clusters{i}))/length(tmpList);
                end

                if sum(overlap>0)>1
                    %clean the list
                    tmpList(ismember(tmpList,treatedIdx,'row'),:) =[];
                    tmpList(~ismember(tmpList,indsCopy,'row'),:) = [];
                    %loop over the list and check for each pixel which cluster
                    %is a better
                    if ~isempty(tmpList)
                       for i = 1:length(tmpList)
%                           corrPx = listCorrPx{inds==tmpList(i)};
%                           currVal = listVal{inds == tmpList(i)};
                          %need to look  in which cluster are corrPx and check
                          %correlation to see which one the current pixels
                          %(tmpList(i) better belongs too)
                          idx = find(overlap>0);
                          corrScores = zeros(length(idx),1);
                          for j = 1:length(nonzeros(idx))
                              %compare current pixel to random subset of
                              %the overlapping clusters to select the best
                              %matching one
                              
                              clusterToTest = clusters{idx(j)};
                              if length(clusterToTest)>10
                                list2Comp = [tmpList(i);clusterToTest(randperm(length(clusterToTest),10))];
                              else
                                list2Comp = [tmpList(i);clusterToTest];
                              end
                              
                              distanceMap = corrAnalysis.getDistanceMapFromPxList(list2Comp,data);

                              corrVal = 1-distanceMap(1,2:end);
                              
                              %before just checking the local pixels
                              %ida = ismember(corrPx,clusters{idx(j)});
                              %corrScores(j) = mean(currVal(ida));
                              corrScores(j) = mean(corrVal);
                          end

                          [~,idb] = max(corrScores);

                          clusters{idx(idb)} = [clusters{idx(idb)}; tmpList(i)];


                       end
                      idx = find(isnan(treatedIdx(:,1)),1);
                      treatedIdx(idx:idx+length(tmpList)-1) = tmpList;
                    else
                    end
                    
                elseif sum(overlap>0) ==1
                    %if there is only a single pixel then we dont consider
                    %it
                    %find with which one it overlaps
                    id = find(overlap>0);
                    group = id;
                     %keep track of the added pixels
                    tmpList(ismember(tmpList,treatedIdx,'row'),:) =[];
                    tmpList(~ismember(tmpList,indsCopy,'row'),:) = [];
        
                    tmpList = corrAnalysis.enforceClusterConsistency(tmpList,clusters{group},data,thresh);
                 
                    clusters{group} = [clusters{group}; tmpList];
                    idx = find(isnan(treatedIdx(:,1)),1);
                    treatedIdx(idx:idx+length(tmpList)-1) = tmpList;
                
                elseif sum(overlap>0) ==0
                    %if there is no overlap then all pixel are added
                    group = length(clusters)+1;

                    tmpList(ismember(tmpList,treatedIdx,'row'),:) =[];
                    tmpList(~ismember(tmpList,indsCopy,'row'),:) = [];
                    
                    clusters{group} = tmpList;

                    idx = find(isnan(treatedIdx(:,1)),1);
                    treatedIdx(idx:idx+length(tmpList)-1) = tmpList;
                end

            else

                 %keep track of the added pixels
                tmpList(ismember(tmpList,treatedIdx,'row'),:) =[];
                tmpList(~ismember(tmpList,indsCopy,'row'),:) = [];
                clusters{group} = [clusters{group}; tmpList];
                idx = find(isnan(treatedIdx(:,1)),1);
                treatedIdx(idx:idx+length(tmpList)-1) = tmpList;
            end
        end
        listCorrPx(inds==currIndex) = [];
        
        meanPx(inds==currIndex) =[];
        listVal(inds==currIndex) = [];
        inds(inds==currIndex) = [];          
        count=count+1;

        if and(count >= safeCount,~isempty(currList))
            error('something went wrong')
        end
                       
    end
    
    for i = 1:length(clusters)
        corrMask(clusters{i}(:,1)) = i;
    
    end
    toc
    