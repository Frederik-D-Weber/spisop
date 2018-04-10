function [sumREM funName remEpochs]= allNightMiniEp1WithREMsFromExactVct(hip, epochLength, exactREMsVct, sampling)


funName = 'allNightR1';
miniEpochSize = 1; % how long our miniepochs should be?
REMsInMiniEpochs = computeMiniEpochsWithREMs(exactREMsVct, sampling, miniEpochSize);


%plot(hip)
%figure
%plot(REMsInMiniEpochs)

miniepPerEpoch = epochLength/miniEpochSize;
%miniepPerEpoch
%length(REMsInMiniEpochs)
%length(exactREMsVct)/300

if round(miniepPerEpoch) ~= miniepPerEpoch
    
    error('Miniepoch size not compatible to epoch size!!!')
     
end

hipLen = length(hip);


if hipLen > length(REMsInMiniEpochs)/miniepPerEpoch
    
   fprintf('Vector is shorter than hipnogram!!!\n');     
   hipLen = floor(length(REMsInMiniEpochs)/miniepPerEpoch);
   
end
%hipLen

remEpochs = 0;
sumREM = 0;

for i = 1:1:hipLen
        
    if hip(i) == 5
             
        sumREM = sumREM + sum( REMsInMiniEpochs((i-1)*miniepPerEpoch + 1 : i*miniepPerEpoch) );
        remEpochs = remEpochs+1;
        
    end        
    
end % for i = 1:1:hipLen
