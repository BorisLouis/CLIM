clc;
close all;
clear;
%% Loading empirical data
path2CorrFile = ['D:\Documents\Unif\PhD\2020-Data\10 - Oct\Sudipta\Noise - Laser Fluctuations'];
laserFluct = load([path2CorrFile filesep 'laserDist.mat']);
laserFluct = laserFluct.lasData;
camNoise   = load([path2CorrFile filesep 'noiseDist.mat']);
camNoise   = camNoise.bkgData;

%% Simulation input
simParam.sizeIm = 128;
simParam.nFrames = 100;
simParam.nParticles = 15;

model.name = 'gaussian';
model.sigma_x = 1.5;
model.sigma_y = 1.5;

simParam.baseCounts = 1e4;
simParam.sdCounts = 0.25;
simParam.baseProb = 0.1;%number between 0 and 0.5
simParam.intMod = 2;
simParam.sdIntMod  = 0.3;
simParam.baseProb = 0.2;
simParam.sdProb = 0.1;
simType = 'blinking'; %'enhancement', 'Bleaching'

%% Simulated intensity profile depending on requested type
%generate two intensity level for each particles with some distribution
simParam.nFrames = simParam.nFrames*3;
intensity = Sim.simIntensity(simParam,simType);
simParam.nFrames = simParam.nFrames/3;
intensity = imresize(intensity,[simParam.nParticles, simParam.nFrames]);
%% add laser noise
lasNoise = laserFluct.center/ laserFluct.Avg;

sampled = datasample(lasNoise,simParam.nFrames,'weight',laserFluct.weight);

sampled = repmat(sampled,simParam.nParticles,1);

intensity = intensity.*sampled;


%% Simulation
sizeIm = simParam.sizeIm;
nFrames = simParam.nFrames;
[X,Y] = meshgrid(1:sizeIm,1:sizeIm);
data = zeros(sizeIm,sizeIm,nFrames,'uint16');
for i = 1:simParam.nParticles
    x0 = randperm(sizeIm,1);
    y0 = randperm(sizeIm,1);
    PSF = Sim.getPSF(X,Y,x0,y0,model);
    for j = 1:nFrames
       
       sPSF = Sim.samplePSF(PSF,intensity(i,j),false);
       data(:,:,j) = data(:,:,j) + sPSF; 
       
    end
    
end

%% add camera noise
noise = datasample(camNoise.center,simParam.sizeIm^2*simParam.nFrames,'weight',camNoise.weight);

noise = reshape(noise,[simParam.sizeIm,simParam.sizeIm,simParam.nFrames]);

data = data+uint16(noise);

%% make Movie
vidFile = VideoWriter('rawMov.mp4','MPEG-4');
vidFile.FrameRate = 10;
open(vidFile);
figure
for i = 1:size(data,3)
   imagesc(data(:,:,i));
   colormap('hot')
   caxis([0 max(data(:))]);
   axis image
   drawnow;
   im = getframe;
   writeVideo(vidFile,im);
   clf
end
close(vidFile)



