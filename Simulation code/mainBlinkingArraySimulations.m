%% Code to simulate different behavior

clear
clc 
close all

%% User input
path2Save = 'D:\Documents\Unif\PhD\2022-Data\01 - Jan\20 - Fluctuation rate';
simParam.sizeIm = 64;
simParam.nFrames = 6000;
simParam.nParticles = 100;

simParam.baseCounts = 3000;
simParam.sdCounts = 0;
simParam.intMod = 1.25;
simParam.sdIntMod  = 0;
simParam.baseProb = 0.05;
simParam.sdProb = 0.05;
simParam.bkgCounts = 500;
simParam.enhancement = 0.5; %in %
simParam.enhTau = 9000;
simType = 'blinking'; %'enhancement', 'bleaching', 'blinking','blinkenha'

delta = 5;
model.name = 'gaussian';
model.sigma_x = 3;
model.sigma_y = 3;

resolution = 10;


%% Simulated intensity profile depending on requested type
%generate two intensity level for each particles with some distribution
%we generate 3 times more frames and then resample to simulate exposure
%time
simParam.nFrames = simParam.nFrames*resolution;
simParam.baseProb = 0.05/resolution;
int = Sim.simIntensity(simParam,simType);
simParam.nFrames = simParam.nFrames/resolution;
intensity = zeros(simParam.nParticles,simParam.nFrames);
for j = 1:simParam.nFrames
    intensity(:,j) = mean(int(:,(j-1)*resolution+1:(j-1)*resolution+resolution),2);

end

% intensity(intensity<simParam.baseCounts) = simParam.baseCounts;
% intensity(intensity>simParam.baseCounts*simParam.intMod) = simParam.baseCounts*simParam.intMod;

%% Simulation
sizeIm = simParam.sizeIm;
nFrames = simParam.nFrames;
[X,Y] = meshgrid(1:sizeIm,1:sizeIm);
data = zeros(sizeIm,sizeIm,nFrames,'uint16');

x = [10:5:55];
y = [10:5:55];

[x,y] = meshgrid(x,y);
PSF = Sim.getPSF(X,Y,simParam.sizeIm/2,simParam.sizeIm/2,model);

for j = 1:nFrames
   idx = sub2ind([size(data,1),size(data,2)],round(y),round(x));
   currConvFrame = data(:,:,j);
   currConvFrame(idx) = intensity(:,j);

   convFrame = conv2(currConvFrame,PSF,'same');
   
   data(:,:,j) = convFrame;

end


%% add camera noise
noise = true;
noiseAmp = [50];

for i = 1:length(noiseAmp)
    tmpData = data;
    if noise
        background = int16(ones(size(tmpData))*simParam.bkgCounts);
        noise2Add = randn(size(tmpData))*noiseAmp(i);
        data2Save = uint16(int16(tmpData)+background+int16(noise2Add));
        
    else
        data2Save = tmpData;
    end



    %% save data
    filename = [path2Save filesep 'mov_' simType '_noiseAmp_' num2str(noiseAmp(i)) '.tif'];
    count = 1;
    while isfile(filename)
        filename = [filename '_' num2str(count),'.tif'];
        count=count+1;
    end

    tiffObj = Tiff(filename,'a');
    tiffObj = dataStorage.writeTiff(tiffObj,data2Save,16);
    tiffObj.close;
end