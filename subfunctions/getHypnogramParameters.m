function [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = getHypnogramParameters(hypn,epochLengthSamples)

hypnStages = [arrayfun(@sleepStage2str,hypn(:,1),'UniformOutput',0) ...
    arrayfun(@sleepStage2str_alt,hypn(:,1),'UniformOutput',0) ...
    arrayfun(@sleepStage2str_alt2,hypn(:,1),'UniformOutput',0)];

hypnEpochs = (1:size(hypn,1));

hypnEpochsBeginsSamples = (((hypnEpochs - 1) * epochLengthSamples) + 1)';
hypnEpochsEndsSamples = (hypnEpochs * epochLengthSamples)';
end