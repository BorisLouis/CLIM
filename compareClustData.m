

relData{1} = relNum0_3_2500;
relData{2} = relNum0_3;
relData{3} = relNum0_5_2500;
relData{4} = relNum0_5;

label{1}   = ['0.3 - 2500'];
label{2}   = ['0.3 - 6000'];
label{3}   = ['0.5 - 2500'];
label{4}   = ['0.5 - 6000'];

corrAnalysis.compareClusters(relData,label);