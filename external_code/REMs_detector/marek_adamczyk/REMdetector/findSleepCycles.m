%function marks sleep cycles in the following way:
% NREMrREM sleep cycles were deﬁned as each REM sleep episode and the NREM 
%sleep immediately preceding it, going back to the sleep onset (1st NREM/REM 
%sleep cycle) or to the limit of another REM sleep episode (from the 2nd 
%NREM/REM sleep cycle to the end of the night). ﬁrst NREM sleep episode began with
%the ﬁrst epoch of stage 2 and ended with the ﬁrst REM sleep episode. Each REM sleep episode
%began with the ﬁrst epoch of REM sleep and ended after the last epoch of REM sleep that was
%followed by at least 15 min of NREM sleep.

%hip - hipnogram with numbers 
%-1 art
% 0 wake
% 1 stage1
% 2 stage2
% 3 stage3
% 4 stage4
% 5 REM

                       %markreaady - in each row of markready we have a sleep cycle. column 1 -
                       %begin, column 2 - REM phase begin, column 3 - end

%function returns matrix with start and end of time fragments in 
%signal to be computed (in seconds nb). To distinguish between fragments, whose results 
%should be gropued, the end of each group is marked by zeros in following
%row. 
%In case of this fun., the result contains thirds of the night, concerning
%sleep time. Sleep is considered from start, when stage 2 is marked for the first
%time.
%hip - hipnogram,
%epochTime - nb of sec in each classificated epoch.

function [result groupsNames] = findSleepCycles(hip, epochLen) 


minNRem = 15*60;%15 minutes
lh = length(hip);

%find begin of 1st cycle
result = find(hip == 2, 1, 'first');
count = 1;
%look for REM episode:
%(Each REM sleep episode began with the ﬁrst epoch of REM sleep and
%ended after the last epoch of REM sleep that was followed by at least 15 min of NREM sleep )

hipSleep = zeros(lh, 1);
for i = 1:1:lh
   
    if hip(i) > 1
       
        hipSleep(i) = 1;
        
    end
    
end

i = result; % begin of 1st cycle


while i <= lh
    %i
    if hip(i) == 5
       
        result(count, 2) = i;
        %the REM has started, now we look for 15 consecutive nREM period
       
        iBg = i;
        while i <= lh
            
            if hip(i) == 5

                %fprintf('dupa\n')
                consnonrem = 0;                
                result(count, 3)=i; %At the beginning in result table, bg of cycle, bg of REM and end of REM is kept
                iBg = i + 1;
            else
                
                consnonrem = consnonrem+epochLen;
                
            end
            
            if consnonrem > minNRem
                count = count+1;
                result(count, 1) = result(count-1, 3)+1;
                break;
            end
           
            i=i+1;
            
        end
%        i=i-1;
        i = iBg;   
        
        
%     elseif i == 233
%         
%         fprintf('ja jebie i=%d \n', i)
%         sum(hipSleep(result(count, 1):i))
%         hipSleep(i)
        
%     elseif sum(hipSleep(result(count, 1):i)) < 4 && ~hipSleep(i)
%         
%         fprintf('dupa i=%d \n', i)
        
    elseif ((sum(hipSleep(result(count, 1):i)) <= ceil((i-result(count, 1) + 1)/2) && sum(hipSleep(result(count, 1):i)) < 8) ||...
            sum(hipSleep(result(count, 1):i)) < 6) && ~hipSleep(i) %thanks to this 
        %when there is not enough st 2 at bg of cycle, because its
        %disturbed by wake or st 1 or mv, we move a few epoch forward(till nxt st 2 occurs)
        

%        i
%        sum(hipSleep(result(count, 1):i)) <= ceil((i-result(count, 1) + 1)/2)
%        sum(hipSleep(result(count, 1):i))
%        ceil((i-result(count, 1) + 1)/2)
        
        if isempty(find(hipSleep(i:lh), 1, 'first'))
            i = lh;
        else
            result(count, 1) = i+ find(hipSleep(i:lh), 1, 'first')-1;
            i = i+ find(hipSleep(i:lh), 1, 'first') - 2;
            
            %newI = i+ find(hipSleep(i:lh), 1, 'first') - 2
        end
    end
    i=i+1;
    
end






if count ~= 1

    if result(count, 3) == 0

        result = result(1:count-1, :);

    end
    
else


    result(1, 2:3) = find(hip > 1, 1, 'last');


end

%result
for i = 1:1:size(result, 1)
   
    groupsNames{i} = strcat('cycle', num2str(i));
    
end

%groupsNames
%result
%StefsNb = 20;

% resultTRY = zeros(size(result, 1), size(result, 2));
% for i = 1:1:size(result, 1)
% 
%     resultTRY(i, 1) = result(i, 1) + find(hip(result(i, 1):lh)==2, 1, 'first') - 1;%result(i, :);
%     resultTRY(i, 2) = resultTRY(i, 1) + find(hip(resultTRY(i, 1):lh)~=2, 1, 'first') - 2;
%     
%     if abs(resultTRY(i, 2) - resultTRY(i, 1)) + 1 > StefsNb
%        
%         resultTRY(i, 2) = resultTRY(i, 1) + StefsNb - 1;
%         
%     end
%     
%     %resTemp(2*i, :) = [0 0];
%     
% end

%result(:, 1:2) = resultTRY(:, 1:2);
result(:, 1:2) = result(:, [1 3]);


%result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%result(:, 2) = result(:, 1) + epochNb - 1;
for i = 1:1:size(result, 1)-1
   
    result(i, 2) = result(i + 1, 1)-1;
    
end

%result(size(result, 1), 2) = lh-1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result(:, 1) = result(:, 1) - 1;
result = result*epochLen;
result(:, 1) = result(:, 1) + 1;
result = result(:, 1:2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultEnd = zeros(size(result, 1)*2-1, size(result, 2));
for kk = 1:2:length(resultEnd)
    
   
    resultEnd(kk, :) = result(floor(kk/2)+1, :);
    
    
end



result = resultEnd;

% if size(result, 1) > 3
%     result = result(1:3, :);
% end
