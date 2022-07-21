%% Simulation experiment #2 - increase of crystal size
%In this experiment we assume a certain trap concentration
%e.g. 1 trap per 500x500x200nm
%to simulate that we increase the amount of charges proportional to the
%volume and the number of traps

clear
clc 
close all
%% User input
path2Save = 'D:\Documents\Unif\PhD\2022-Data\04 - April\26 Blinking Models';

initCount = 1000;
initialVolume = 0.5*0.5*0.2;

simParam.trapCapacity = 200;
simParam.sdCapacity  = 0.1;
%0.05 probability is the standard (=1switch every 20 frames = 1 sec)
simParam.onProb = 0.1; %here on/off time refers to the trap being active or not, so it is reverse on the blinking
simParam.offProb = 0.1;
simParam.sdProb = 0.5;

simParam.nSim = 50;
simParam.nFrames = 6000;
simParam.bkgCounts = 500;
resolution = 10;
simParam.onProb = simParam.onProb /resolution;
simParam.offProb = simParam.offProb /resolution; 

statsTS = table(zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),...
    zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),'VariableNames',...
    {'width','med','absAmp','relAmp','widthNN','medNN','absAmpNN','relAmpNN'});
statsInt = statsTS;
allData = zeros(simParam.nSim,simParam.nFrames,2);
currentVolume = initialVolume;
simParam.baseCounts = initCount;

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

for i = 1:simParam.nSim
    currentVolume = initialVolume*(i);
    simParam.baseCounts = initCount*(i);
    simParam.nTraps = uint16(1*currentVolume/initialVolume);

    simParam.nFrames = simParam.nFrames*resolution;
    
    [trapSatIntensity,interactionIntensity] = Sim.trapSim(simParam);
    simParam.nFrames = simParam.nFrames/resolution;
    trapSatIntensity = imresize(trapSatIntensity,[1, simParam.nFrames]);
    interactionIntensity = imresize(interactionIntensity,[1, simParam.nFrames]);   
    
    noise = true;
    noiseAmp = [50];

    for j = 1:length(noiseAmp)
        tmpTs = trapSatIntensity;
        tmpInt = interactionIntensity;
        if noise
            background = int32(ones(size(tmpTs))*simParam.bkgCounts);
            noise2Add = randn(size(tmpTs))*noiseAmp(j);
            tsIntensity2Save = int32(int32(tmpTs)+background+int32(noise2Add));
            intIntensity2Save = int32(int32(tmpInt)+background+int32(noise2Add));
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
    statsTS(i,:).med = median(tsIntensity2Save);
    statsTS(i,:).absAmp = (max(tsIntensity2Save)-min(tsIntensity2Save));
    statsTS(i,:).relAmp = (max(tsIntensity2Save)-min(tsIntensity2Save))/max(tsIntensity2Save);
    
    trapSatIntensity = double(trapSatIntensity);     
    statsTS(i,:).widthNN = std(trapSatIntensity);
    statsTS(i,:).medNN = median(trapSatIntensity);
    statsTS(i,:).absAmpNN = (max(trapSatIntensity)-min(trapSatIntensity));
    statsTS(i,:).relAmpNN = (max(trapSatIntensity)-min(trapSatIntensity))/max(trapSatIntensity);
    
    intIntensity2Save = double(intIntensity2Save);
    statsInt(i,:).width = std(double(intIntensity2Save));
    statsInt(i,:).med = median(intIntensity2Save);
    statsInt(i,:).absAmp = (max(intIntensity2Save)-min(intIntensity2Save));
    statsInt(i,:).relAmp = (max(intIntensity2Save)-min(intIntensity2Save))/max(intIntensity2Save);
    
    interactionIntensity = double(interactionIntensity);
    statsInt(i,:).widthNN = std(double(interactionIntensity));
    statsInt(i,:).medNN = median(interactionIntensity);
    statsInt(i,:).absAmpNN = (max(interactionIntensity)-min(interactionIntensity));
    statsInt(i,:).relAmpNN = (max(interactionIntensity)-min(interactionIntensity))/max(interactionIntensity);
        
    corrOutput.statsTS = statsTS;
    corrOutput.statsInt = statsInt;
end