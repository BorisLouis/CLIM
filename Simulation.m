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
sizeIm = 128;
nFrames = 100;
data = zeros(sizeIm,sizeIm,nFrames);
nParticles = 15;

model.name = 'gaussian';
model.sigma_x = 2;
model.sigma_y = 2;

counts = 1e3;
sdCounts = 0.25;
blinkMod = 2;

sim = 'blinking'; %'enhancement', 'Bleaching'



%% Simulation
for i = 1:nParticles
    x0 = randperm(sizeIm,1);
    y0 = randperm(sizeIm,1);
    secondInt = 2*BaseInt;
    BaseInt = counts;
    int = BaseInt;
    
    
    for j = 1:nFrames
       
        
        
    end
    
end
