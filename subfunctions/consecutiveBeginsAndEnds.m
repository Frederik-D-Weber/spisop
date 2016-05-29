function [begins, ends] =  consecutiveBeginsAndEnds(x,maxConsecutiveDiff)
    nx = length(x);
    iter = 1;
    begins(iter) = x(1);
    ends = [];
    if nx > 1
        for i = 2:nx
            if x(i) > (x(i-1) + maxConsecutiveDiff)
                ends(iter) = x(i-1);
                iter = iter + 1;
                begins(iter) = x(i);
            end
        end
    end
    if length(ends) < length(begins)
        ends(iter) = x(nx);
    end
end