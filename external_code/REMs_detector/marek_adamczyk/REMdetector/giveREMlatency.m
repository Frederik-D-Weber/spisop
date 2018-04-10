function [REMlatency funName remEpochs] = giveREMlatency(hip, epochLength, exactREMsVct, sampling)



funName = 'REMlatency';


%hipLen = length(hip);

remEpochs = 0;

%REMlatency = 0;

REMstarts = find(hip == 5, 1, 'first');
sleepStarts = find(hip > 1, 1, 'first');

REMlatency = REMstarts - sleepStarts;


% for i = 1:1:hipLen
%         
%     if hip(i) == 5
%         
%         allREMsNb = allREMsNb + 1;
%         remEpochs = remEpochs + 1;
%         
%     end
%     
% end