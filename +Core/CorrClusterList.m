classdef CorrClusterList < handle
    % CorrClusterSummary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nCluster
        clusters
        statistics
        
    end
    
    methods
        function obj = CorrClusterList(data,inds)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            clusters = cell(size(inds));
            for i = 1:size(inds)
                
               [row,col] = ind2sub(size(data),inds(i));
               
               clusters{i} = Core.CorrCluster(inds(i),i,squeeze(data(row,col,:)));
                
                
            end
            
            obj.clusters = clusters;
            obj.nCluster = length(obj.clusters);
            
        end
        
        function updateState(obj)
%             nClust = obj.nCluster;
%             
%             stats = table(zeros(nClust,1),zeros(nClust,1),zeros(nClust,1),'VariableNames',{'minCorr','meanCorr','nPx'});
%             for i = 1:nClust
%                 currClust = obj.clusters{i};
%                 stats(i,:) = currClust.getStats;               
%                 
%                 
%             end
%             
%             obj.statistics = table(mean(stats.minCorr),mean(stats.meanCorr),mean(stats.nPx),'VariableNames',{'minCorr','meanCorr','nPx'});

              obj.nCluster = length(obj.clusters);
            
        end
        
        
        
        function didMerge = merge(obj,indA,indB,distanceMap)
            %indA, indB are in distance map space
            %idA, idB are in data space
            %if all pixel are part of the analysis then they are equal.
            
            didMerge = false;
            doMerge = true;
            
            % Get the two clusters
            clustA = obj.clusters{indA};
            clustB = obj.clusters{indB};
            
            indsA = clustA.getInds;
            indsB = clustB.getInds;
            
            corrRelA = clustA.corrRel;
            corrRelB = clustB.corrRel;
            
            % Check if merging is worth
            avgTrace = (clustA.nPx*clustA.avgTrace + clustB.nPx*clustB.avgTrace)/(clustB.nPx+clustA.nPx);
            
            idx = nchoosek([indsA(:,2);indsB(:,2)],2);
            
            ind = sub2ind(size(distanceMap),idx(:,1),idx(:,2));
            
            newCorrRel = distanceMap(ind);
            
            if or(size(indsA,1)==1,size(indsB,1)==1)
                doMerge = true;
            else
                meanCorrDistA = mean(corrRelA);
                meanCorrDistB = mean(corrRelB);
                meanM = (meanCorrDistA*clustA.nPx + meanCorrDistB*clustB.nPx)/(clustA.nPx+clustB.nPx);
               
                maxCorrDistAB = max(newCorrRel);
                
                sil = (maxCorrDistAB-meanM)/max([maxCorrDistAB,meanM]);
                
                if sil <0.7
                    doMerge = true;
                else 
                    doMerge = false;
                end
            end
            
            
            if doMerge
                %add ClustB into ClustA
                newInds = [indsA;indsB];
                newTrace = avgTrace;
                
                newClust = Core.CorrCluster(newInds(:,1),newInds(:,2),newTrace,newCorrRel);

                obj.clusters(indA) = [];
                obj.clusters(indB) = [];
                obj.clusters{end+1} = newClust;
                
                didMerge = true;
            else
                
                a=5;
                
            end
            
                        
        end
        
        function [corrMask] = getCorrMask(obj,dim)
           
            corrMask = zeros(dim(1),dim(2));
            clustList = obj.clusters;
            
            for i = 1:length(clustList)
                currClust = clustList{i};
                
                corrMask(currClust.inds(:,2)) = i;
                
                
                
            end
            
            
            
        end
                  
    end
end

