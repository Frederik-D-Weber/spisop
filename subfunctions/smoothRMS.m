function r = smoothRMS(x,wd)
l = length(x);
bu = floor(wd/2);
r = zeros(1,l);
for i=1:l
    if (i<=bu)
        r(i) = rms(x(1:i));
    elseif i > (l-bu)
        r(i) = rms(x(i:end));
    else
        r(i) = rms(x((i-bu):(i+bu)));
    end
end