function [numberECG] = getECGChannelNumbers(datachannel)
% datachannel = data.label
ecg_channel = ft_channelselection('*ECG*',datachannel);
numberECG = -1;
if ~isempty(ecg_channel)
    numberECG = find(strcmp(ecg_channel(1),datachannel));
else
    ecg_channel = ft_channelselection('*EKG*',datachannel);
    if ~isempty(ecg_channel)
        numberECG = find(strcmp(ecg_channel(1),datachannel));
    end
end

end