function [sumOfBursts funName remEpochs]= allNightREM1NbOfMidBurstsFromExactVct(hip, epochLength, exactREMsVct, sampling)


funName = 'allNightR1MidBNb';
miniEpochSize = 1; % how long our miniepochs should be?
REMsInMiniEpochs = computeMiniEpochsWithREMs(exactREMsVct, sampling, miniEpochSize);

distanceWithinBurst = 2; % how many miniepochs can be empty btween miniepochs with REMs in order to clissify them as one burst?
minNbOfREMs = 3;

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
%sumInBurst = 0;
thisBurst = 0;

sumOfBursts = 0;

previousMinis = -1000;
lastInBurst = 0;
%FollowMinis = 1000;

for i = 1:1:hipLen
        
    if hip(i) == 5
 
        
        for jj = (i-1)*miniepPerEpoch + 1 : 1 : i*miniepPerEpoch
        
            if REMsInMiniEpochs(jj)
               
                sumREM = sumREM + 1;

                if jj - previousMinis <= distanceWithinBurst + 1 % is it possible, that current REM belongs to the burst in which previous REM is ? (if it is at all) 
                   
                    thisBurst = thisBurst + 1;
                    
                    if lastInBurst == 0
                    
                        thisBurst = thisBurst + 1;
                        
                    end
                    
                    lastInBurst = 1;
                    
                else % new burst (if it is at all)
                    
                    lastInBurst = 0;
                    
                    if thisBurst >= minNbOfREMs % Does previous REM belong to a burst according to defined conditions?, if yes, save this burst data
                                            
%                        sumInBurst = sumInBurst + thisBurst;
                        sumOfBursts = sumOfBursts + 1;
                    end

                    thisBurst = 0;
                        
                end
                
                previousMinis = jj;
                    
            end
              
        end % for jj = (i-1)*miniepPerEpoch + 1 : 1 : i*miniepPerEpoch
        
        
        
%        sumREM = sumREM + sum( REMsInMiniEpochs((i-1)*miniepPerEpoch + 1 : i*miniepPerEpoch) );
        remEpochs = remEpochs+1;
        
    end % if hip(i) == 5
        
end % for i = 1:1:hipLen

    
if thisBurst >= minNbOfREMs % Does previous REM belong to a burst according to defined conditions?, if yes, save this burst data
                                            
    %sumInBurst = sumInBurst + thisBurst;
    sumOfBursts = sumOfBursts + 1;
    
end

