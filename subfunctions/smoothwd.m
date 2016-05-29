%smoothing,  with sliding window of wd data points, i.e. the center sample +-floor(wd/2), 
%therefore wd is recomended to be an odd positive integer and smaller than points in the data. 
%At the boundaries the mean is smoothed using a reduced window size that shrinks to the maximum samples left over to the left or right of the center sample.
%smooth(mean(tfa.powspctrm,1),round(SmoothFreqWindowSize/(1/SegmentLength)),'lowess');

function r = smoothwd(x,wd)
if (mod(wd,2) == 0) || (wd < 3) || (wd > length(x))
    error('the window size wd is supposed to be an positive odd number greater than 2, but smaller than the number of data points x!');
end

l = length(x);
r = zeros(1,l);
lastSample = l - wd + 1;
tempAddISmpl = wd - 1;
tempWindows = zeros(lastSample,wd);

for iSmpl = 1:lastSample
    tempWindows(iSmpl,:) = x(iSmpl:(iSmpl + tempAddISmpl));
end

firstSmplstimeWndw = (wd+1)/2 ;
lastSmplstimeWndw = l - ((wd-1)/2);
r(firstSmplstimeWndw:lastSmplstimeWndw) = mean(tempWindows,2);

r(1:(firstSmplstimeWndw-1)) = smoothLeftBorder(x(1:((firstSmplstimeWndw*2)-1)),wd);
r((lastSmplstimeWndw+1):end) = smoothRightBorder(x(((lastSmplstimeWndw-firstSmplstimeWndw)+1):end),wd);

end


function r = smoothRightBorder(x,wd)
bu = floor(wd/2);
r = zeros(1,bu);
for i=1:bu
    r(i) = mean(x((2*i-1):(2*bu-1)));
end
end


function r = smoothLeftBorder(x,wd)
    r = smoothRightBorder(x(end:-1:1),wd);
    r = r(end:-1:1);
end