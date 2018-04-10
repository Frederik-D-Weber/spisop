% function gets a vector window (probably some signal or derivative of
% that) and looks in the area around jj, to finds a place, where value is
% lower than than thrMin for the first time (explanation for positive 
% values, but works also for negative-> in this case higher than...)
% as a result, function returns length of the area, where values - thrMin 
%have the same sign

function res = checkRange(jj, window, thrMin)


value = window(jj);
wl = length(window);

if value > 0 
    window = window - thrMin;
else
    window = window + thrMin;
end



for i = jj:-1:1
    
    bg = i;
    if abs(window(i) + value) < abs(value)
  
        break;
        
    end
   
end

for i = jj:1:wl
    
    eend = i;
    if abs(window(i) + value) < abs(value)
       
        break;
        
    end
   
end

%value
%bg
%eend

res = eend - bg;