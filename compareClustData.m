

for i = 1:length(relNum1)
   relData{i} = relNum1SilClean(i);
   lab{i}   = num2str(thresh(i));
   
    
end

corrAnalysis.compareClusters(relData,lab);

% relData{1} = relNum0_3_2500;
% relData{2} = relNum0_3;
% relData{3} = relNum0_5_2500;
% relData{4} = relNum0_5;
% 
% label{1}   = ['0.3 - 2500'];
% label{2}   = ['0.3 - 6000'];
% label{3}   = ['0.5 - 2500'];
% label{4}   = ['0.5 - 6000'];
% 
% corrAnalysis.compareClusters(relData,label);


% 
% relData{1} = relNum1(11);
% relData{2} = relNum1Clean(11);
% label{1}   = ['Normal'];
% label{2}   = ['Cleaned'];
% corrAnalysis.compareClusters(relData,label);
%% best output
clear relData;
clear lab
met1 = sil.*treatedArea';
[~,idx] = max(met1);
relData{1} = relNum1(idx);
lab{1} = 'Normal';

met2 = cleanSil.*cleanTreatedArea';
[~,idx] = max(met2);
relData{2} = relNum1Clean(idx);
lab{2} = 'Area Clean';

met3 = silClean.*silCleanTreatedArea';
[~,idx] = max(met3);
relData{3} = relNum1SilClean(idx);
lab{3} = 'SilClean';


corrAnalysis.compareClusters(relData,lab);







%%
relData{1} = Pseudo;
relData{2} = Watershed;
relData{3} = Kmean;
%relData{4} = DBScan;

lab{1} = 'Pseudo';
lab{2} = 'Watershed';
lab{3} = 'Kmean';
%label{4} = 'DBScan';

corrAnalysis.compareClusters(relData,lab);
