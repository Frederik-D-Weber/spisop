function eogNoise = determineSignalNoise( EOGsignal, orgSampling, eogClass )

sampling = 100;
EOGsignal = resample(EOGsignal, sampling, orgSampling);

HdPass = load('filterFirBandPass100_6_7_to_16_17.mat');
bandPassFltrdEog = filter(HdPass.Hd, EOGsignal);
bandPassFltrdEog = bandPassFltrdEog(ceil(length(HdPass.Hd.numerator)/2):length(bandPassFltrdEog));


stdTable = zeros(ceil(length(bandPassFltrdEog)/sampling), 1);
count = 0;

for i = 1:sampling:length(bandPassFltrdEog) - sampling

    placeInClass = floor(i/sampling) + 1; %in which second are we????
    
    if length(eogClass) < placeInClass || eogClass( placeInClass ) < 1        
        continue;        
    end
    
    piece = bandPassFltrdEog(i:i+sampling - 1);
    count = count + 1;
    stdTable(count) = std(piece);


end


stdTable = stdTable(1:count);



eogNoise = median(stdTable);

