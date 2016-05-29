function [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(filepath,epochLengthSamples)
hypn = load(filepath);

[hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = getHypnogramParameters(hypn,epochLengthSamples);
end

