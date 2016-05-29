function [begins, ends] = getROIsByHypnogram(hypnFilepath,epochLengthSamples,sleepStagesOfInterst)
[hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnFilepath,epochLengthSamples);

columnOfInterestHypnStages = 1; 
if ismember('SWS',sleepStagesOfInterst)
    columnOfInterestHypnStages = 2; 
elseif ismember('NonREM',sleepStagesOfInterst)
    columnOfInterestHypnStages = 3; 
end

epochsOfInterst = hypnEpochs(ismember(hypnStages(:,columnOfInterestHypnStages),sleepStagesOfInterst) & (hypn(:,2) == 0));

nEpochsOfInterest = length(epochsOfInterst);

[consecBegins, consecEnds] = consecutiveBeginsAndEnds(epochsOfInterst,1);

begins = hypnEpochsBeginsSamples(consecBegins);
ends = hypnEpochsEndsSamples(consecEnds);

end