function [begins, ends] = getBeginsAndCorrespondingEndsThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow)
    tx = getThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow);
    if length(tx) >= 2
        if mod(length(tx),2) ~= 0
            tx = tx(1:end-1);
        end
        begins = tx(1:2:end);
        ends = tx(2:2:end);
    else
        begins = [];
        ends = [];
    end
end


function tx = getThresholdCrossing(data,threshold,isRisingAboveNotFallingBelow)
% get indices of threshold crossing in data starting with a rising/falling edge,
%i.e. from below to above threshold
    data = data - threshold;
    temp = data(1:end-1) .* data(2:end); % scalar multiplication of two arrays. Second array is shifted to left by 1.
    p = find(temp<0);
    if length(p) >= 2
        if isRisingAboveNotFallingBelow
            if(data(p(1))<0) % transition from below to above.
                tx = p(1:1:end);
            else             % transition from above to below.
                tx = p(2:1:end);
            end
        else
            if(data(p(1))<0) % transition from above to below.
                tx = p(2:1:end);
            else             % transition from below to above.
                tx = p(1:1:end);
            end
        end
        tx=tx';
    else
        tx = [];
    end
    
end