%% Code to simulate different behavior

clear
clc 
close all
addpath('E:\Users\Boris\Documents\MATLAB\CLIM')
%% User input
path2Save = 'D:\Documents\Unif\PostDoc\2024 - Data\10 - October\Resolution limit';
simParam.sizeIm = 64;
simParam.nFrames = 1000;
simParam.nParticles = 2;

simParam.baseCounts = 3000;
simParam.sdCounts = 0;
simParam.intMod = 1.5;
simParam.sdIntMod  = 0;
simParam.baseProb = 0.05;
simParam.sdProb = 0;
simParam.bkgCounts = 500;
simParam.enhancement = 0.5; %in %
simParam.enhTau = 9000;
simType = 'blinking'; %'enhancement', 'bleaching', 'blinking','blinkenha'

model.name = 'gaussian';
FWHM = 3;


distance = [2 1 0.75, 0.5, 0.25, 0.1];%expresses in FWHM

resolution = 10;%time resolution for the simulation before resampling


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

intensity(intensity<simParam.baseCounts) = simParam.baseCounts;
intensity(intensity>simParam.baseCounts*simParam.intMod) = simParam.baseCounts*simParam.intMod;

%% Simulation
simRFac = 10;
model.sigma_x = FWHM/2.35*simRFac;
model.sigma_y = FWHM/2.35*simRFac;
for i = 1:length(distance)
    
    delta= distance(i)*FWHM*simRFac;
    nFrames = simParam.nFrames;
    data = zeros(simParam.sizeIm,simParam.sizeIm,nFrames,'uint16');
    sizeIm = simParam.sizeIm*simRFac;
    
    [X,Y] = meshgrid(1:sizeIm,1:sizeIm);
    

    x = [12*simRFac 12*simRFac+delta];
    y = [16*simRFac 16*simRFac];

    PSF = Sim.getPSF(X,Y,sizeIm/2,sizeIm/2,model);

    for j = 1:nFrames
        currConvFrame = zeros(simParam.sizeIm*simRFac,simParam.sizeIm*simRFac,'uint16');
        idx = sub2ind([size(currConvFrame,1),size(currConvFrame,2)],round(y),round(x));

        currConvFrame(idx) = intensity(:,j);

        convFrame = conv2(currConvFrame,PSF,'same');

        data(:,:,j) = imresize(convFrame,1/simRFac);
       
    end
    noiseAmp = [50];
    background = int16(ones(size(data))*simParam.bkgCounts);
    noise2Add = randn(size(data))*noiseAmp;
    data2Save = uint16(int16(data)+background+int16(noise2Add));

      %% save data
    filename = [path2Save filesep 'mov_' simType '_FWHM_' num2str(FWHM) '_delta_' num2str(distance(i)) '.tif'];
    count = 1;
    while isfile(filename)
    filename = [filename '_' num2str(count),'.tif'];
    count=count+1;
    end

    tiffObj = Tiff(filename,'a');
    tiffObj = dataStorage.writeTiff(tiffObj,data2Save,16);
    tiffObj.close;

end




 
