function [satTrapIntensity, interactionIntensity,state] = trapSim(simParam)
%store useful parameters
nTraps = simParam.nTraps;

baseCounts = simParam.baseCounts;
trapCapacity = simParam.trapCapacity;
nFrames = simParam.nFrames;

onProb = simParam.onProb;
offProb = simParam.offProb;

%generate traps
trapList = cell(nTraps,1);
initialState = round(rand(1,nTraps));
capacityList = zeros(1,nTraps);
for i = 1:nTraps
    newTrap = Sim.Trap(trapCapacity,onProb,offProb,initialState(i));
    capacityList(i) = newTrap.capacity;
    trapList{i} = newTrap;
end

I0 = baseCounts;

I1 = I0 - sum(initialState(:).*capacityList(:));

capacityList2 = ones(size(capacityList)) * simParam.trapCapacity2;

I2 = I0 * (I0/(I0+sum(initialState(:).*capacityList2(:))));
%throw dices to know if particles will switch or not
satTrapIntensity = zeros(1,nFrames);
interactionIntensity = zeros(1,nFrames);
satTrapIntensity(:,1) = I1;
interactionIntensity(:,1) = I2;
prevState = initialState;

rng('shuffle')
dice = rand(nFrames,nTraps);
state= zeros(1,nFrames);
state(1) = sum(initialState/nTraps);
for i=2:nFrames
   currentState = prevState;
   for j = 1:nTraps
      currentTrap = trapList{j};
      currentTrap.doSwitch(dice(i,randi(nTraps)));
      currentState(j) = currentTrap.getState; 
      
   end  
   
   satTrapIntensity(:,i) = I0 - sum(currentState(:).*capacityList(:));
   interactionIntensity(:,i) = I0 * (I0/(I0+sum(currentState(:).*capacityList2(:))));
   
   prevState = currentState;
   
    state(i) = sum(currentState/nTraps);
   
end
              
end