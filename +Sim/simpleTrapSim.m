function Intensity = simpleTrapSim(simParam)
%store useful parameters
nTraps = simParam.nTraps;

baseCounts = simParam.baseCounts;
trapCapacity = simParam.trapCapacity;
nFrames = simParam.nFrames;

onProb = simParam.onProb;
offProb = simParam.offProb;

initialState = round(rand(1,nTraps));
I0 = baseCounts;

I1 = I0 - sum(initialState)*trapCapacity;

%throw dices to know if particles will switch or not
prob = rand(nTraps,nFrames);
Intensity = zeros(1,nFrames);
Intensity(:,1) = I1;

prevState = initialState;
for i=2:nFrames
    
   idx2Switch = prob(:,i) <= onProb;

   currentState = prevState;
   
   currentState(idx2Switch) = ~currentState(idx2Switch);
   
   
   Intensity(:,i) = I0 - sum(currentState)*trapCapacity;
   
   prevState = currentState;
end
              
end