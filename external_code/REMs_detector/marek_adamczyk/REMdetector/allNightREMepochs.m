function [allREMsNb funName remEpochs] = allNightREMepochs(hip, epochLength, exactREMsVct, sampling)



funName = 'allNightREMep';


hipLen = length(hip);

allREMsNb = 0;
remEpochs = 0;

allREMsNb = length(find(hip == 5));
remEpochs = length(find(hip == 5));

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