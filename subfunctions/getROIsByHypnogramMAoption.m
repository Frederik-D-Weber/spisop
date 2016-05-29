function [begins, ends] = getROIsByHypnogramMAoption(hypnFilepath,epochLengthSamples,sleepStagesOfInterst,includeMAepochs)
[hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnFilepath,epochLengthSamples);

columnOfInterestHypnStages = 1; 
if ismember('SWS',sleepStagesOfInterst)
    columnOfInterestHypnStages = 2; 
elseif ismember('NonREM',sleepStagesOfInterst)
    columnOfInterestHypnStages = 3; 
end


epochsOfInterst = hypnEpochs(ismember(hypnStages(:,columnOfInterestHypnStages),sleepStagesOfInterst) & ((hypn(:,2) == 0) | includeMAepochs));

nEpochsOfInterest = length(epochsOfInterst);

[consecBegins, consecEnds] = consecutiveBeginsAndEnds(epochsOfInterst,1);

begins = hypnEpochsBeginsSamples(consecBegins);
ends = hypnEpochsEndsSamples(consecEnds);

end