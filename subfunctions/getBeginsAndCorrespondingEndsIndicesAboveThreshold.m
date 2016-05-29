function [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(data,threshold)
    [begins, ends] = getBeginsAndCorrespondingEndsThresholdCrossing(data,threshold,logical(1));
    if length(begins) >= 1 
        begins = begins + 1;
    end
end