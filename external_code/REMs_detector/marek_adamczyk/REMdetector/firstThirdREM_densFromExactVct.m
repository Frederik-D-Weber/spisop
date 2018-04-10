function [avREMsNb funName remEpochs] = firstThirdREM_densFromExactVct(hip, epochLength, exactREMsVct, sampling)


funName = '1stThirdR3D';
miniEpochSize = 3; % how long our miniepochs should be?
REMsInMiniEpochs = computeMiniEpochsWithREMs(exactREMsVct, sampling, miniEpochSize);


miniepPerEpoch = epochLength/miniEpochSize;


if round(miniepPerEpoch) ~= miniepPerEpoch
    
    error('Miniepoch size not compatible to epoch size!!!')
     
end

hipLen = length(hip);

if hipLen > length(REMsInMiniEpochs)/miniepPerEpoch
    
   fprintf('Vector is shorter than hipnogram!!!\n');     
   hipLen = floor(length(REMsInMiniEpochs)/miniepPerEpoch);
   
end


[result groupsNames] = returnThirdsfromPSG(hip, epochLength);
placeInHip = floor(result/epochLength);

begin = placeInHip(1, 1);
eend = placeInHip(1, 2) ;





remEpochs = 0;
sumREM = 0;

for i = begin:1:eend
    
    
    if hip(i) == 5
             
        sumREM = sumREM + sum( REMsInMiniEpochs((i-1)*miniepPerEpoch + 1 : i*miniepPerEpoch) );
        remEpochs = remEpochs+1;
        
    end
        
    
end % for i = 1:1:hipLen

if remEpochs
   
    avREMsNb = sumREM / remEpochs;
    
else 
    
    avREMsNb = nan;
    
end



