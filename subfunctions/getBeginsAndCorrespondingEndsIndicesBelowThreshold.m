function [begins, ends] = getBeginsAndCorrespondingEndsIndicesBelowThreshold(data,threshold)
    [begins, ends] = getBeginsAndCorrespondingEndsThresholdCrossing(data,threshold,logical(0));
    if length(begins) >= 1 
        begins = begins + 1;
    end
end