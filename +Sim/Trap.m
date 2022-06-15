classdef Trap < handle
    properties
        capacity
        onProb
        offProb
        state
    end
    methods
        %constructor
        function obj = Trap(capacity,onProb,offProb,state)
            obj.capacity = capacity;
            assert(and(offProb>1/1e9,onProb>1/1e9),'Probability too low for rng')
            obj.onProb = onProb;
            obj.offProb = offProb;
            
            obj.state = state;
        end
        
        function state = getState(obj)
            state = obj.state;
            
        end
        
        function [dice] = doSwitch(obj,dice)
            currentState = obj.state;
            %choose probability based on current state
            switch currentState
                case 1
                    switchProb = obj.offProb;
                case 0 
                    switchProb = obj.onProb;
                otherwise
                    error('impossible state')
            end
            
            %toss a coin to know if it switches or not
           
          
            if dice <= switchProb
                obj.state = ~currentState;
            else 
                obj.state = currentState;
            end           
        end
        
        
        
    end
end