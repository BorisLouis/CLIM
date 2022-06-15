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

simParam.trapCapacity = 900;
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


statsTS = table(zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),'VariableNames',...
    {'nTraps','volume','width','med','absAmp','relAmp'});
statsInt = statsTS;
allData = zeros(simParam.nSim,simParam.nFrames,2);
currentVolume = initialVolume;
simParam.baseCounts = initCount;

corrOutput = struct();
corrOutput.results = struct();
corrOutput.results(simParam.nSim).traceTs = 0;
corrOutput.results(simParam.nSim).traceInt = 0;
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
            background = uint32(ones(size(tmpTs))*simParam.bkgCounts);
            noise2Add = randn(size(tmpTs))*noiseAmp(j);
            tsIntensity2Save = uint32(uint32(tmpTs)+background+uint32(noise2Add));
            intIntensity2Save = uint32(uint32(tmpInt)+background+uint32(noise2Add));
        else
            tsIntensity2Save = trapSatIntensity;
            intIntensity2Save = interactionIntensity;
        end
    end
        
    allData(i,:,1) = tsIntensity2Save;
    allData(i,:,2) = intIntensity2Save;
    corrOutput.results(i).traceTs = tsIntensity2Save;
    corrOutput.results(i).traceInt = intIntensity2Save;
    
    
    tsIntensity2Save = double(tsIntensity2Save);     
    statsTS(i,:).width = std(tsIntensity2Save);
    statsTS(i,:).med = median(tsIntensity2Save);
    statsTS(i,:).absAmp = (max(tsIntensity2Save)-min(tsIntensity2Save));
    statsTS(i,:).relAmp = (max(tsIntensity2Save)-min(tsIntensity2Save))/max(tsIntensity2Save);
    statsTS(i,:).volume = currentVolume;
    statsTS(i,:).nTraps = simParam.nTraps;
    
    intIntensity2Save = double(intIntensity2Save);
    statsInt(i,:).width = std(double(intIntensity2Save));
    statsInt(i,:).med = median(intIntensity2Save);
    statsInt(i,:).absAmp = (max(intIntensity2Save)-min(intIntensity2Save));
    statsInt(i,:).relAmp = (max(intIntensity2Save)-min(intIntensity2Save))/max(intIntensity2Save);
    statsInt(i,:).volume = currentVolume;
    statsInt(i,:).nTraps = simParam.nTraps;
    
    corrOutput.statsTS = statsTS;
    corrOutput.statsInt = statsInt;
end