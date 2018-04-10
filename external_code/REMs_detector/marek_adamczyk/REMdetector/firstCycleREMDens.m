function [avREMsNb remEpochs] = firstCycleREMDens(hip, REMsData)



%resLen = min([length(hip) length(REMsData)]);

count = 0;
sumREM = 0;

[result groupsNames] = findSleepCycles(hip, 30);

%result
%plot(hip)
%pause


temp = floor(result/30);
begin = temp(1, 1);
eend = temp(1, 2);


for i = begin:1:eend
   
    if hip(i) == 5
       
      sumREM = sumREM + REMsData(i);
      count = count+1;
        
    end
    
end

if count
   
    avREMsNb = sumREM / count;
    
else 
    
    avREMsNb = nan;
    
end

remEpochs = count;
