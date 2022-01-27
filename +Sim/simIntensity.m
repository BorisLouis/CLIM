function Intensity = simIntensity(simParam,simType)
    baseCounts = simParam.baseCounts;
    sdCounts = simParam.sdCounts;
    intMod = simParam.intMod;
    sdIntMod  = simParam.sdIntMod;
    baseProb = simParam.baseProb;
    sdProb   = simParam.sdProb;
    nParticles = simParam.nParticles;
    nFrames = simParam.nFrames;
    
    switch lower(simType)
        
        case 'blinking'
            
            I0 = baseCounts-sdCounts*baseCounts + rand(nParticles,1)*2*sdCounts*baseCounts;
            
            blinkingF = intMod-sdIntMod+rand(nParticles,1)*2*sdIntMod;
            I1 = I0.*blinkingF;
            %give a switching probability to each particles with some distribution
            switchProb = baseProb -sdProb + rand(nParticles,1)*2*sdProb;

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
            
        case 'enhancement'
            
            x  = 1:nFrames;
            x  = repmat(x,nParticles,1);
            x0 = 1 + round(rand(nParticles,1)*8);
            I0 = baseCounts-sdCounts*baseCounts + rand(nParticles,1)*2*sdCounts*baseCounts;
            %get intensity modifier and calculate I1 the final intensity
            %level
            intModifier = intMod-sdIntMod+rand(nParticles,1)*2*sdIntMod;
            I1 = I0.*intModifier;
            %get lifetime
            tau = baseProb -sdProb + rand(nParticles,1)*2*sdProb;
            tau = tau*nFrames;
               
            Intensity = I0 + (I1-I0).* (1 - exp(-((x-x0)./tau))); 
            
        case 'bleaching'
            x  = 1:nFrames;
            x  = repmat(x,nParticles,1);
            x0 = 1 + round(rand(nParticles,1)*8);
            I0 = baseCounts-sdCounts*baseCounts + rand(nParticles,1)*2*sdCounts*baseCounts;
            
            intModifier = intMod-sdIntMod+rand(nParticles,1)*2*sdIntMod;
            %for bleaching we divide by the intensity modifier
            I1 = I0./intModifier;
            
            tau = baseProb -sdProb + rand(nParticles,1)*2*sdProb;
            tau = tau*nFrames;
               
            Intensity = I0 + (I1-I0).* (1 - exp(-((x-x0)./tau)));
            
        case 'blinkenha'
            I0 = baseCounts-sdCounts*baseCounts + rand(nParticles,1)*2*sdCounts*baseCounts;
            
            blinkingF = intMod-sdIntMod+rand(nParticles,1)*2*sdIntMod;
            I1 = I0.*blinkingF;
            %give a switching probability to each particles with some distribution
            switchProb = baseProb -sdProb + rand(nParticles,1)*2*sdProb;

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
            
            x  = 1:nFrames;
            x  = repmat(x,nParticles,1);
            x0 = 1 + round(rand(nParticles,1)*8);
            I0 = 1;
            %get intensity modifier and calculate I1 the final intensity
            %level
            I1 = I0 + simParam.enhancement;
            %get lifetime
            tau = simParam.enhTau;               
            enhancement = I0 + (I1-I0).* (1 - exp(-((x-x0)./tau)));
            
            for i = 1:size(Intensity,1)
               
                Intensity(i,:) = Intensity(i,:).*enhancement(i,:);
                
            end
            
        otherwise
            error('Simulation option requested does not exist or is not implemented yet');
    end
end