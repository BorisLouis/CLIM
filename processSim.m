% Code to process results of simulations
clear;
close;
clc;

%% User input
path = 'D:\Documents\Unif\PhD\Papers\13 - IntensityCorrletation\Simulations\test';
nParticles = 15;
%% Load Data

folder2Analyze = dir(path);
folder2Analyze(1:2) = [];

nFolder = length(folder2Analyze);

%% Process data
delta = 5;
FWHM_pix = 4;
chi2 = 24;
allData = cell(1,nFolder);
for i = 1:nFolder
   
    currPath = [folder2Analyze(i).folder filesep folder2Analyze(i).name filesep 'Experiment.mat'];
    
    %get the saved experiment object
    currExp  = load(currPath);
    currExp  = currExp.myExperiment;
    
    %get movies
    mov = fieldnames(currExp.corrMovies);
    
    data = table(zeros(length(mov),1), zeros(length(mov),1), zeros(length(mov),1),...
    zeros(length(mov),1),zeros(length(mov),1),zeros(length(mov),1),...
    'VariableNames',{'noise','amp','SNRamp','SNR','nPart','detRat'});

    for j = 1:length(mov)
        currMov = currExp.corrMovies.(mov{j});
        frames = currMov.info.frame2Load;
        
        data = currMov.getFrame(frames);
        
        im = data(:,:,1);
        
        [pos, ~, ~] = Localization.smDetection( double(im), delta, FWHM_pix, chi2 );
        pos = round(pos);
        
        ROIs = zeros(size(im));
        iTrace = zeros(size(pos,1),size(data,3));
        
        for k = 1: size(pos,1)
            
            idxRow = pos(k,1)-delta:pos(k,1)+delta;
            idxCol = pos(k,2)-delta:pos(k,2)+delta;          
            iTrace(k,:) = squeeze(sum(sum(data(idxRow,idxCol,:),1),2));        
            ROIs(idxRow,idxCol) = 1;
            
        end
        
        iTrace = iTrace - mean(data(ROIs==0))*(2*delta+1)^2;
        noise = std(double(data(ROIs==0)));
        minState = mean(mink(iTrace,10,2));
        maxState = mean(maxk(iTrace,10,2));
        
        Amp = maxState-minState;
        
        SNRAmp = Amp/noise;
        SNR    = mean(iTrace(:))/noise;
        nPart  = size(iTrace,1);
        detectionRatio = nPart/nParticles;
        
        data(j).noise = noise;
        data(j).amp = amp;
        data(j).SNRamp = SNRamp;
        data(j).SNR = SNR;
        data(j).nPart = nPart;
        data(j).detRat = detRat;
        
    end
    
    allData{i} = data;
    
end
