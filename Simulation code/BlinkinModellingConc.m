%% Simulation experiment #1 - increasing of trap concentration
%In this experiment we assume that the volume of the crystal and the amount
%of charges available to quench is constant. We slowy increase the number
%of traps and thus the concentration
clear
clc 
close all
%% User input
path2Save = 'D:\Documents\Unif\PhD\2022-Data\04 - April\26 Blinking Models';
simParam.baseCounts = 10000;
simParam.trapCapacity = 2000;
simParam.nSim = 50;
simParam.nTraps = 0; %0 is for getting 1 simulation per number of traps, otherwise it is a fix number
simParam.sdCapacity  = 0.1;

%0.05 probability is the standard (=1switch every 20 frames = 1 sec)
simParam.onProb = 0.07; %here on/off time refers to the trap being active or not, so it is reverse on the blinking
simParam.offProb = 0.07;
simParam.sdProb = 0;% % of fluctuation of on/off prob around their nomitive value
%put to 0 to obtain constant on/prob value


simParam.nFrames = 6000;
simParam.bkgCounts = 500;
resolution = 10;
simParam.onProb = simParam.onProb /resolution;
simParam.offProb = simParam.offProb /resolution;

nTrapList = [1:simParam.nSim];
% nTrapList = [1, 2, 3, 4, 5, 8, 10, 15, 20, 25, 30, 50, 100, 150, 200, 250, 500, 1000, 1500, 2000];
% nTrapList = [5, 5, 5, 5, 10, 10, 10, 10, 20, 20, 20, 20, 30, 30, 30, 30, 50, 50, 50, 50];


%%

statsTS = table(zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),...
    zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),'VariableNames',...
    {'width','med','absAmp','relAmp','widthNN','medNN','absAmpNN','relAmpNN','nTraps'});
statsInt = statsTS;
allData = zeros(simParam.nSim,simParam.nFrames,4);
corrOutput = struct();
corrOutput.results = struct();
corrOutput.results(simParam.nSim).traceTs = 0;
corrOutput.results(simParam.nSim).traceInt = 0;

syms x
eqn =  simParam.baseCounts * (simParam.baseCounts/(simParam.baseCounts+x)) == simParam.baseCounts-simParam.trapCapacity ;
S= solve(eqn,x);
simParam.trapCapacity2 = double(S);

symObj = syms;
cellfun(@clear,symObj)

if simParam.nTraps <= 0
    mode = 'multi';
    
else
    mode = 'single';
end

for i = 1:min([length(nTrapList),simParam.nSim])
    switch mode
        case 'multi'
            simParam.nTraps = nTrapList(i);
        case 'single'
        otherwise
            error('Unknown trap simulation mode');
    end
   
    simParam.nFrames = simParam.nFrames*resolution;
    
    [trapSatIntensity,interactionIntensity,state] = Sim.trapSim(simParam);
    simParam.nFrames = simParam.nFrames/resolution;
    
    trace = zeros(simParam.nFrames,1);
    for j = 1:simParam.nFrames
        trace(j) = mean(interactionIntensity((j-1)*resolution+1:(j-1)*resolution+resolution));

    end
    interactionIntensity = trace;
    
    
    trace = zeros(simParam.nFrames,1);
    for j = 1:simParam.nFrames
        trace(j) = mean(trapSatIntensity((j-1)*10+1:(j-1)*10+1+9));

    end
    trapSatIntensity = trace;
    
    noise = true;
    noiseAmp = [50];

    for j = 1:length(noiseAmp)
        tmpTs = trapSatIntensity;
        tmpInt = interactionIntensity;
        if noise
            background = int16(ones(size(tmpTs))*simParam.bkgCounts);
            noise2Add = randn(size(tmpTs))*noiseAmp(j);
            tsIntensity2Save = int16(int16(tmpTs)+background+int16(noise2Add));
            intIntensity2Save = int16(int16(tmpInt)+background+int16(noise2Add));
        else
            tsIntensity2Save = trapSatIntensity;
            intIntensity2Save = interactionIntensity;
        end
    end
    
    %add bkg to initial trace
    trapSatIntensity = trapSatIntensity+double(background);
    interactionIntensity = interactionIntensity+double(background);
   
    
    allData(i,:,1) = tsIntensity2Save;
    allData(i,:,2) = intIntensity2Save;
    allData(i,:,3) = trapSatIntensity;
    allData(i,:,4) = interactionIntensity;
    
    corrOutput.results(i).traceTs = tsIntensity2Save;
    corrOutput.results(i).traceInt = intIntensity2Save;
    corrOutput.results(i).traceTsNN = trapSatIntensity;
    corrOutput.results(i).traceIntNN = interactionIntensity;
    corrOutput.results(i).nTrap = simParam.nTraps;
    
    tsIntensity2Save = double(tsIntensity2Save);     
    statsTS(i,:).width = std(tsIntensity2Save);
    statsTS(i,:).med = mean(tsIntensity2Save);
    statsTS(i,:).absAmp = (mean(tsIntensity2Save(tsIntensity2Save>prctile(tsIntensity2Save,90)))...
        -mean(tsIntensity2Save(tsIntensity2Save<prctile(tsIntensity2Save,10))));
    statsTS(i,:).relAmp = statsTS(i,:).absAmp/mean(tsIntensity2Save(tsIntensity2Save>prctile(tsIntensity2Save,90)));
    
    trapSatIntensity = double(trapSatIntensity);     
    statsTS(i,:).widthNN = std(trapSatIntensity);
    statsTS(i,:).medNN = mean(trapSatIntensity);
    statsTS(i,:).absAmpNN = (mean(trapSatIntensity(trapSatIntensity>prctile(trapSatIntensity,90)))...
        -mean(trapSatIntensity(trapSatIntensity<prctile(trapSatIntensity,10))));
    statsTS(i,:).relAmpNN = statsTS(i,:).absAmpNN*(mean(trapSatIntensity(trapSatIntensity>prctile(trapSatIntensity,90))));
    statsTS(i,:).nTraps = simParam.nTraps;
    
    intIntensity2Save = double(intIntensity2Save);
    statsInt(i,:).width = std(double(intIntensity2Save));
    statsInt(i,:).med = mean(intIntensity2Save);
    statsInt(i,:).absAmp = (mean(intIntensity2Save(intIntensity2Save>prctile(intIntensity2Save,90)))...
        -mean(intIntensity2Save(intIntensity2Save<prctile(intIntensity2Save,10))));
    statsInt(i,:).relAmp = statsInt(i,:).absAmp/(mean(intIntensity2Save(intIntensity2Save>prctile(intIntensity2Save,90))));
    
    interactionIntensity = double(interactionIntensity);
    statsInt(i,:).widthNN = std(double(interactionIntensity));
    statsInt(i,:).medNN = mean(interactionIntensity);
   
    statsInt(i,:).absAmpNN = (mean(interactionIntensity(interactionIntensity>prctile(interactionIntensity,90)))...
        -mean(interactionIntensity(interactionIntensity<prctile(interactionIntensity,10))));
    statsInt(i,:).relAmpNN = statsInt(i,:).absAmpNN/(mean(interactionIntensity(interactionIntensity>prctile(interactionIntensity,90))));
    statsInt(i,:).nTraps = simParam.nTraps;
    corrOutput.statsTS = statsTS;
    corrOutput.statsInt = statsInt;
   
end