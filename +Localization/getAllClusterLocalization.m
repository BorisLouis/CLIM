function [allLoc] = getAllClusterLocalization(obj,data)
           
            nCluster = obj.corrMask.rawNCluster;
           allLoc = cell(nCluster,1);
           for i = 1 :nCluster
               
               loc = obj.getClusterLocalization(data,i);
               allLoc{i} = loc;
               
           end
           obj.clustLoc = allLoc;
        end