function Intensity = simIntensity(simParam,simType)
    baseCounts = simParam.baseCounts;
    sdCounts = simParam.sdCounts;
    blinkMod = simParam.intMod;
    sdBlink  = simParam.sdIntMod;
    baseProb = simParam.baseProb;
    sdProb   = simParam.sdProb;
    nParticles = simParam.nParticles;
    nFrames = simParam.nFrames;
    
    switch lower(simType)
        case 'blinking'
            I0 = baseCounts-sdCounts*baseCounts + rand(1,nParticles)*2*sdCounts*baseCounts;
            I0 = I0';
            blinkingF = blinkMod-sdBlink+rand(1,nParticles)*2*sdBlink;
            I1 = I0.*blinkingF;
            %give a switching probability to each particles with some distribution
            switchProb = baseProb -sdProb + rand(1,nParticles)*2*sdProb;

            %throw dices to know if particles will switch or not
            Intensity = rand(nParticles,nFrames);
            Intensity(:,1) = I0;

            prevIntensity = Intensity(:,1);
            for i=2:nFrames
               idxToSwitch = Intensity(:,i)<=switchProb(:);

               %The one that do not switch get same intensity as before
               Intensity(~idxToSwitch,i) = prevIntensity(~idxToSwitch);

               %The one that switch, we need to check what was the previous level
               idxI0 = prevIntensity == I0(:);

               idxI0ToSwitch = logical(idxToSwitch.*idxI0);
               idxI1ToSwitch = logical(idxToSwitch.*~idxI0);

               %switch the intensity level that needs switch
               Intensity(idxI0ToSwitch,i) = I1(idxI0ToSwitch);
               Intensity(idxI1ToSwitch,i) = I0(idxI1ToSwitch);  

               prevIntensity = Intensity(:,i);
            end
        otherwise
            error('Simulation option requested does not exist or is not implemented yet');
    end
end