% function takes a vector with exactly marked REMs and returns veectors
% with miniepochs marked with 0 (no REMs in miniepoch) and 1(There are REMs in miniepoch)

%miniEpochSize: Time length (miniepochSize = 1 -> 1sec)


function REMsInMiniEpochs = computeMiniEpochsWithREMs(exactlyMarkedREMs, sampling, miniEpochSize)

REMsInMiniEpochs = zeros(ceil(length(exactlyMarkedREMs)/sampling/miniEpochSize), 1);


for i = 1:1:length(exactlyMarkedREMs)    
    
    if exactlyMarkedREMs(i)
       
        REMsInMiniEpochs(floor(i/sampling/miniEpochSize) + 1) = 1;        
                
    end
        
end











