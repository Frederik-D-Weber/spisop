function [numberEEG numberEOG numberEMG] = getScoringChannelNumbers(datachannel)

eeg_channel = ft_channelselection('*C3*',datachannel);
eog_channel = ft_channelselection('*EOG*',datachannel);
emg_channel = ft_channelselection('*EMG*',datachannel);

if ~isempty(eeg_channel)
    numberEEG = find(strcmp(eeg_channel(1),datachannel));
else
    numberEEG = 2;
end

if ~isempty(eog_channel)
    numberEOG = find(strcmp(eog_channel(1),datachannel));
else
    numberEEG = 1;
end

if ~isempty(emg_channel)
    numberEMG = find(strcmp(emg_channel(1),datachannel));
else
    numberEEG = numel(datachannel);
end

end

