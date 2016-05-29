%rms of data x with sliding window of wd data points, i.e. the center sample +-floor(wd/2), 
%therefore wd is recomended to be an odd positive integer and smaller than points in the data. 
%At the boundaries the RMS is smoothed using a reduced window size that shrinks to the maximum samples left over to the left or right of the center sample.

function r = smoothRMSwd(x,wd)
if (mod(wd,2) == 0) || (wd < 3) || (wd > length(x))
    error('the window size wd is supposed to be an positive odd number greater than 2, but smaller than the number of data points x!');
end
l = length(x);
r = zeros(1,l);
lastRMSSample = l - wd + 1;
tempAddISmpl = wd - 1;
tempRMSWindows = zeros(lastRMSSample,wd);

for iSmpl = 1:lastRMSSample
    tempRMSWindows(iSmpl,:) = x(iSmpl:(iSmpl + tempAddISmpl));
end

firstSmplsRMStimeWndw = (wd+1)/2 ;
lastSmplsRMStimeWndw = l - ((wd-1)/2);
r(firstSmplsRMStimeWndw:lastSmplsRMStimeWndw) = rms(tempRMSWindows,2);

r(1:(firstSmplsRMStimeWndw-1)) = smoothLeftBorderRMS(x(1:((firstSmplsRMStimeWndw*2)-1)),wd);
r((lastSmplsRMStimeWndw+1):end) = smoothRightBorderRMS(x(((lastSmplsRMStimeWndw-firstSmplsRMStimeWndw)+1):end),wd);

end


function r = smoothRightBorderRMS(x,wd)
bu = floor(wd/2);
r = zeros(1,bu);
for i=1:bu
    r(i) = rms(x((2*i-1):(2*bu-1)));
end
end


function r = smoothLeftBorderRMS(x,wd)
    r = smoothRightBorderRMS(x(end:-1:1),wd);
    r = r(end:-1:1);
end




            
            