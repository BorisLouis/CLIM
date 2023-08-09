function [satTrapIntensity, interactionIntensity,state] = trapSim2(simParam)
%store useful parameters
nTraps = simParam.nTraps;

baseCounts = simParam.baseCounts;
trapCapacity = simParam.trapCapacity;
nFrames = simParam.nFrames;

%generate traps
trapList = cell(nTraps,1);
initialDice = rand(1,nTraps);
probOn = simParam.onProb/(simParam.onProb+simParam.offProb);
initialState = initialDice<probOn;

capacityList = zeros(1,nTraps);
for i = 1:nTraps
    onProb = simParam.onProb(i) + (-1 + 2*rand(1))*simParam.sdProb*simParam.onProb(i);
    offProb = simParam.offProb(i) + (-1 + 2*rand(1))*simParam.sdProb*simParam.offProb(i);
    
    newTrap = Sim.Trap(trapCapacity(i),onProb,offProb,initialState(i));
    capacityList(i) = newTrap.capacity;
    trapList{i} = newTrap;
end

I0 = baseCounts;

I1 = I0 - sum(initialState(:).*capacityList(:));

capacityList2 = ones(size(capacityList)) .* simParam.trapCapacity2;

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
      currentTrap.doSwitch(dice(i,j));%randi(nTraps)));
      currentState(j) = currentTrap.getState; 
      
   end  
%    if all(currentState == 0)
%        disp('stop')
%    end
   satTrapIntensity(:,i) = I0 - sum(currentState(:).*capacityList(:));
   interactionIntensity(:,i) = I0 * (I0/(I0+sum(currentState(:).*capacityList2(:))));
   
   prevState = currentState;
   
    state(i) = sum(currentState/nTraps);
   
end
              
end