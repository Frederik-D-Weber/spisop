%function estimates where searching for REMs makes sense, looking for
%uncorelated fragments without artifacts, where eog are likely. This
%version should be robust enough to handle people with EEG disturbances

%eogClass - result, where each second of the signal is scored

%example: 
%eogClass = eogClinicClassify(signal(:, 1), signal(:, 2), 250, hip, true);


function eogClass = eogClinicClassifyFlexThr(eog1, eog2, sampling, hip, verbose)


%--------------------------------
%used only for performance monitoring ("verbose" option)
%values which mark possible artifacts, 
highAmplitude = 1000;
highFrequency = 100;
highCorrelation = 10;
%--------------------------------

%sampling I want -> I do not need high freq and this improves the speed
mySampling = 100;

%tic
fprintf('  Ori REMsMaAd code: resampling... \n')
eog1 = resample(eog1, mySampling, sampling);
eog2 = resample(eog2, mySampling, sampling);
%toc

sampling = mySampling;


sigLen = length(eog1);

% here our results will be kept, by default everything is fine. If
% something wrong is found in the signal, given parameter will be
% substracted
eogClass = ones(ceil(sigLen/sampling), 1);



eogNoise1 = determineSignalNoise( eog1, mySampling, ones(length(eogClass), 1) );
eogNoise2 = determineSignalNoise( eog1, mySampling, ones(length(eogClass), 1) );

eogNoise = (eogNoise1 + eogNoise2)/2;




%tic
%1. dont classify fragments with huge amplitude
fprintf('  Ori REMsMaAd code: Looking for high amplitudes... \n')

mxAmpEOG = 200*eogNoise;%850;

count = 6;
for i = 5*sampling+1:sampling:sigLen - 11*sampling

    sig1 = abs(eog1(i:i+sampling-1));
    sig2 = abs(eog2(i:i+sampling-1));
    
    
    if max(sig1) > mxAmpEOG || max(sig2) > mxAmpEOG  
       
        eogClass(count-5:count+10) = eogClass(count-5:count+10) - highAmplitude;
        
    end
    
    count = count + 1;

end


%toc

%2. dont classify fragments affected with high frequency
fprintf('  Ori REMsMaAd code: Looking for high frequencies... \n');

%tic

HdHigh = load('filterFirHighPass100_19_20.mat');
fprintf('  Ori REMsMaAd code: highpass filtering eog1... \n');
eog1Filterd = filter(HdHigh.HdFirHighPass100_19_20, eog1);

fprintf('  Ori REMsMaAd code: highpass filtering eog2... \n');
eog2Filterd = filter(HdHigh.HdFirHighPass100_19_20, eog2);

%get rid of fragment where is nothing because filtering starts
eog1Filterd = eog1Filterd(ceil(length(HdHigh.HdFirHighPass100_19_20.numerator)/2):length(eog1Filterd));
eog2Filterd = eog2Filterd(ceil(length(HdHigh.HdFirHighPass100_19_20.numerator)/2):length(eog2Filterd));

HdHigh = load('filterFirLowPass100_48_45.mat');
fprintf('  Ori REMsMaAd code: lowpass filtering eog1... \n');
eog1Filterd = filter(HdHigh.HdFirLowPass100_48_45, eog1Filterd);

fprintf('  Ori REMsMaAd code: lowpass filtering eog2... \n');
eog2Filterd = filter(HdHigh.HdFirLowPass100_48_45, eog2Filterd);

%get rid of fragment where is nothing because filtering starts
eog1Filterd = eog1Filterd(ceil(length(HdHigh.HdFirLowPass100_48_45.numerator)/2):length(eog1Filterd));
eog2Filterd = eog2Filterd(ceil(length(HdHigh.HdFirLowPass100_48_45.numerator)/2):length(eog2Filterd));


freqThr = 1.5*eogNoise;%6;

count = 3;
range = 2*sampling;

for i = range+1:sampling:length(eog1Filterd) - range - sampling

    % to better estimate freq with fft, remove linear trends 
    eogFragment1 = detrend(eog1Filterd(i:i + sampling-1));
    eogFragment2 = detrend(eog2Filterd(i:i + sampling-1));

    if std(eogFragment1) > freqThr || std(eogFragment2) > freqThr
      
        for jj = count-1:1:count+1
            eogClass(jj) =  eogClass(jj) - highFrequency;
        end
        
    end

    
    %try also 0.5sec later
    eogFragment1 = detrend(eog1Filterd(i + floor(0.5*sampling):i + sampling - 1 + floor(0.5*sampling)));
    eogFragment2 = detrend(eog2Filterd(i + floor(0.5*sampling):i + sampling - 1 + floor(0.5*sampling)));
    
    if std(eogFragment1) > freqThr || std(eogFragment2) > freqThr
      
        for jj = count-1:1:count+1
            eogClass(jj:jj+1) =  eogClass(jj:jj+1) - highFrequency;
        end
        
    end
    
    
    count = count + 1;

end

%toc



%3. dont classify highly correlated fragments
fprintf('  Ori REMsMaAd code: Looking for highly correlated fragments... \n');

% Im interested in correlation of waves, bot I want to alod very low
% freq artifacts which should not affect eye movements
HdHigh = load('filterFirHighPass100_03_05.mat');
%tic
fprintf('  Ori REMsMaAd code: highpass filtering eog1... \n')
eog1Filtered = filter(HdHigh.HdFirHighPass100_03_05, eog1);
%toc

%tic
fprintf('  Ori REMsMaAd code: highpass filtering eog2... \n')
eog2Filtered = filter(HdHigh.HdFirHighPass100_03_05, eog2);
%toc

%get rid of fragment where is nothing because filtering starts
eog1Filtered = eog1Filtered(floor(length(HdHigh.HdFirHighPass100_03_05.numerator)/2) : length(eog1Filtered));
eog2Filtered = eog2Filtered(floor(length(HdHigh.HdFirHighPass100_03_05.numerator)/2) : length(eog2Filtered));
toAdd = floor(length(HdHigh.HdFirHighPass100_03_05.numerator)/2) - 1;
eog1Filtered = cat(1, eog1Filtered, eog1Filtered(length(eog1Filtered))*ones(toAdd, 1));
eog2Filtered = cat(1, eog2Filtered, eog2Filtered(length(eog2Filtered))*ones(toAdd, 1));


secNb = 5; %correlation is checked on 5 sec gap
gap = secNb*sampling;


corrThr = 0.6;
highVarianceCorrThr = 0.15;
count = 1;


% one could also try to filter out very low freq (below 0.2) and try again,
% but everything should be rechecked again (influence in REMs also)
for i = 1:sampling:length(eog1) - gap
   
    temp = myCorr(eog1Filtered(i:i+gap-1), eog2Filtered(i:i+gap-1));

%    if temp >= corrThr || (temp >= highVarianceCorrThr && min(std(eog1(i:i+gap-1)), std(eog2(i:i+gap-1))) > 15)
    if temp >= corrThr || (temp >= highVarianceCorrThr && min(std(eog1(i:i+gap-1)), std(eog2(i:i+gap-1))) > 3.75*eogNoise)
        
        eogClass(count+1:count+3) = eogClass(count+1:count+3) - highCorrelation + temp;
        
    end
    
    count = count + 1;
end




if ~isempty(hip) && verbose % hipnogram is only for talking

    %hip = hypno2vec(fileName);
    hipClas = zeros(ceil(sigLen/sampling), 1);
        
    count = 1;
    eplen = 30; 
    for i = 1:1:length(hip) - 2
    
        if length(hipClas) >=count + eplen-1      
            hipClas(count:count + eplen-1) = hip(i)/4;        
        end
        
        tm = getTime(i, eplen);
        fprintf('  Ori REMsMaAd code: time: %s, epoch: %d   ', tm, hip(i))
        fprintf('  Ori REMsMaAd code: %2.2f ', eogClass(count:count + eplen-1));
        fprintf('\n');
        
        count = count + eplen;
        
    end  
    
    plot(eogClass)
    hold on
    plot(hipClas, 'r')
    
end


%%
%function returns time as string taking the number of the epoch and epoch length
function tm = getTime(num, epLen)

num = num-1;%begin of theh epoch
num = num*epLen;
hr = floor(num/3600);

num = mod(num, 3600);
mn = floor(num/60);
sc = mod(num,60);

tm = strcat(num2str(hr, '%02d'), ':',num2str(mn, '%02d'), ':',num2str(sc, '%02d'));


