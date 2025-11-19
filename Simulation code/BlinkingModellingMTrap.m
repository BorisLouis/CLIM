%% Simulation experiment #3 - constant concentration and volume
%In this experiment we assume a certain trap concentration
%e.g. 1 trap per 250x250x250nm grain
%to simulate that, we increase the amount of charges proportional to the
%volume and the number of traps

clear
clc 
close all
addpath('E:\Users\Boris\Documents\MATLAB\CLIM');
%% User input
path2Save = 'C:\Users\Boris\OneDrive - KU Leuven\Documents\Unif\Postdoc\Paper\01 - Film vs individual grains\March 2025\Figure\Simulation - trapN2';

initCount = 4500;%30000
initialVolume = 0.25*0.25; %um^2 assuming constant thickness

nTraps = 10;% 75
%0.05 probability is the standard (=1switch every 20 frames = 1 sec)
trapType = 20;
simParam.trapCapacity = [linspace(3500,500,trapType)];% [ones(trapType,1)*26400];%[linspace(4000,1000,trapType)];%[ones(trapType,1)*1500]; %Oxygen: 4000;4000;4000;4000];%, 1000 ];%5000];%] 5000];
% simParam.onProb = [0.001 0.01 0.05 0.1 0.3 ];%0.3 0.005]; %here on/off time refers to the trap being active or not, so it is reverse on the blinking
% simParam.offProb =[0.003 0.01 0.05 0.1 0.3 ];% 0.3 0.005];

simParam.onProb = [logspace(-6,-0.8,trapType)];% Oxygen0.1 0.1 0.1 0.1];
simParam.offProb = simParam.onProb;

simParam.sdProb = 0;
trapList = [1, 1, 1];
simParam.nSim = 100;
simParam.nFrames = 6000;
simParam.bkgCounts = 500; %without noise, no background
resolution = 10;
simParam.onProb = simParam.onProb /resolution;
simParam.offProb = simParam.offProb /resolution; 

statsTS = table(zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),...
    zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),zeros(simParam.nSim,1),'VariableNames',...
    {'width','med','absAmp','relAmp','widthNN','medNN','absAmpNN','relAmpNN','nTraps'});
statsInt = statsTS;
allData = zeros(simParam.nSim,simParam.nFrames,2);
currentVolume = initialVolume;
simParam.baseCounts = initCount;

%% Distribute onOff efficiencies
randIdx = randi(length(simParam.trapCapacity),nTraps);

trapCapacity = simParam.trapCapacity(randIdx(1,1:end));
onProb = simParam.onProb(randIdx(1,1:end));
offProb = simParam.offProb(randIdx(1,1:end));

simParam.trapCapacity = trapCapacity;
simParam.onProb = onProb;
simParam.offProb = offProb;

corrOutput = struct();
corrOutput.results = struct();
corrOutput.results(simParam.nSim).traceTs = 0;
corrOutput.results(simParam.nSim).traceInt = 0;
assert(and(length(simParam.trapCapacity) == nTraps,length(simParam.onProb) == nTraps),'Something is wrong with the dimensions of the simulation parameters, check that the parameters are consistent with the number of traps requested')

for i = 1:nTraps
    syms x
    eqn =  simParam.baseCounts * (simParam.baseCounts/(simParam.baseCounts+x)) == simParam.baseCounts-simParam.trapCapacity(i) ;
    S= solve(eqn,x);
    simParam.trapCapacity2(i) = double(S);
    symObj = syms;
    cellfun(@clear,symObj)
end

for i = 1:simParam.nSim
    currentVolume = initialVolume;
    simParam.baseCounts = initCount;
    simParam.nTraps = nTraps;

    simParam.nFrames = simParam.nFrames*resolution;
       
    [trapSatIntensity,interactionIntensity,state] = Sim.trapSim2(simParam);
    simParam.nFrames = simParam.nFrames/resolution;
    
    trace = zeros(simParam.nFrames,1);
    for j = 1:simParam.nFrames
        trace(j) = mean(interactionIntensity((j-1)*10+1:(j-1)*10+1+9));

    end
    interactionIntensity = trace;
    
    
    trace = zeros(simParam.nFrames,1);
    for j = 1:simParam.nFrames
        trace(j) = mean(trapSatIntensity((j-1)*10+1:(j-1)*10+1+9));

    end
    trapSatIntensity = trace;
    noise = true;
    noiseAmp = [10];

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
    corrOutput.results(i).volume = currentVolume;
    corrOutput.results(i).charges = simParam.baseCounts;
    
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