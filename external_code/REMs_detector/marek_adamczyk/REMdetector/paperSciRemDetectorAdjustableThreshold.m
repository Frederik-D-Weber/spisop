%function computes REMs
%eog1, eog2 - signals
%sampling - sampling rate
%epLen - length of the epoch
%eogClass - result of eogClassify function
%countAll = false : divide epoch on 3 sec fragments and add all frag. with REMs included
%countAll = true : compute all REMs in the epoch


%example:
%remDens = REMDetectorOnClassificatedAllSig(test1(:, 1), test1(:, 2), 250, 30, eogClass, false, false);

% remDens: each entry is an epoch, each value is number of miniepochs (3
% sec) in a given epoch with REMs (standard MPI rem density definition)

function [remDens, exactlyMarkedREMs, meanREMSpeed, meanREMWidth, EOG1, EOG2, EOG1_lp, EOG2_lp, EOG1_lphp, EOG2_lphp, EOG1_lp_dev, EOG2_lp_dev, EOG1_lphp_dev, EOG2_lphp_dev] ...
    = paperSciRemDetectorAdjustableThreshold(EOG1, EOG2, orgSampling, epLen, eogClass, countAll,epochOfInterestToPlot, RelaxedThresholdFactorOfBasicThreshold, BasicThresholdFactorOfEOGnoise)

epochOfInterest = epochOfInterestToPlot;
%sampling I want -> I do not need high freq and this improves the speed
sampling = 100;
miniEp = epLen / 10;

%tic
fprintf('resampling... \n')
EOG1 = resample(EOG1, sampling, orgSampling);
EOG2 = resample(EOG2, sampling, orgSampling);
%toc

% eog1 = originalEog1;
% eog2 = originalEog2;




HdLow = load('filterFirLowPass100_6_5.mat');
%tic
fprintf('  Ori REMsMaAd code: lowpass filtering eog1... \n')
EOG1_lp = filter(HdLow.HdFirLowPass100_6_5, EOG1);
%toc

%tic
fprintf('  Ori REMsMaAd code: lowpass filtering eog1... \n')
EOG2_lp = filter(HdLow.HdFirLowPass100_6_5, EOG2);
%toc

%get rid of fragment where is nothing because filtering starts

EOG1_lp = EOG1_lp(ceil(length(HdLow.HdFirLowPass100_6_5.numerator)/2):length(EOG1_lp));
EOG2_lp = EOG2_lp(ceil(length(HdLow.HdFirLowPass100_6_5.numerator)/2):length(EOG2_lp));


HdHigh = load('filterFirHighPass100_03_05.mat');
%tic
fprintf('  Ori REMsMaAd code: highpass filtering eog1... \n')
EOG1_lphp = filter(HdHigh.HdFirHighPass100_03_05, EOG1_lp);
%toc

%tic
fprintf('  Ori REMsMaAd code: highpass filtering eog2... \n')
EOG2_lphp = filter(HdHigh.HdFirHighPass100_03_05, EOG2_lp);
%toc

%get rid of fragment where is nothing because filtering starts
EOG1_lphp = EOG1_lphp(floor(length(HdHigh.HdFirHighPass100_03_05.numerator)/2) : length(EOG1_lphp));
EOG2_lphp = EOG2_lphp(floor(length(HdHigh.HdFirHighPass100_03_05.numerator)/2) : length(EOG2_lphp));



%fltr = ones(floor(0.1*sampling), 1)*(1/floor(0.1*sampling)); %->sprawdz jak jest w tym najnowszym papierze, cos mi sie wydaje ze podobnie(przejrzyj notatki)
%fltrLen = length(fltr);
gap = epLen*sampling;
sigLen = length(EOG1_lp);

%extend signals a bit, so we will have always counting up to last epoch
toAdd = gap - mod(length(sigLen), gap) + 1;

EOG1_lp = cat(1, EOG1_lp, EOG1_lp(length(EOG1_lp))*ones(toAdd, 1));
EOG2_lp = cat(1, EOG2_lp, EOG2_lp(length(EOG2_lp))*ones(toAdd, 1));

sigLen = length(EOG1_lp);

toAdd = sigLen - length(EOG1_lphp);
EOG1_lphp = cat(1, EOG1_lphp, EOG1_lphp(length(EOG1_lphp))*ones(toAdd, 1));
EOG2_lphp = cat(1, EOG2_lphp, EOG2_lphp(length(EOG2_lphp))*ones(toAdd, 1));


if sigLen ~= length(EOG2_lp)

    error('Ori REMsMaAd code: Chosen EOG channels do not have the same length')

end

exactlyMarkedREMs = zeros(sigLen, 1);


fprintf('  Ori REMsMaAd code: REM detection... \n');
%tic

EOG1_lp_dev = diff(EOG1_lp);%derivative
EOG2_lp_dev = diff(EOG2_lp);%derivative

EOG1_lphp_dev = diff(EOG1_lphp);%derivative
EOG2_lphp_dev = diff(EOG2_lphp);%derivative



% length(Y1)
% length(Y1Filtered)
% sigLen
% gap
% pause

remDens = zeros(1,floor(sigLen/gap));

meanREMSpeed = zeros(1,floor(sigLen/gap));
meanREMWidth = zeros(1,floor(sigLen/gap));



eogNoise1 = determineSignalNoise( EOG1, sampling, eogClass );
eogNoise2 = determineSignalNoise( EOG2, sampling, eogClass );

eogNoise = (eogNoise1 + eogNoise2)/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%epochOfInterest = 100000000000; %if user would like to see exactly how the chosen epoch is analysed, please specify nunber of this epoch here

% values are "hardly typed" -> because signal is resampled, adjustment is
% unneccessary


%thrBS = 2.61; %basic steepness threshold
%thrRS = 1.65; % relaxed steepness threshold
%RelaxedThresholdFactorOfBasicThreshold = 0.66;
%BasicThresholdFactorOfEOGnoise = 0.7;

thrBS = BasicThresholdFactorOfEOGnoise*eogNoise; %basic steepness threshold
thrRS = RelaxedThresholdFactorOfBasicThreshold*thrBS; % relaxed steepness threshold, 0.462

thrRange = 2; % gap range (corresponding to 2*remSearchGap = 140 ms) 
thrMin2 = 0;%thr minimal to interpret increase/decrease as significant
remSearchGap = floor(0.07*sampling); % 70msec
clearNbrhThr = 1.35; %threshold for neighbourhood criterion

classCount = 1;
count = 1;
for i = 1:gap:sigLen - gap % -> for each epoch

    %Let m = length(u) and n = length(v) . Then w is the vector of length
    %m+n-1 (vector after convoluting)

    %    c1 = cat(1, zeros(ceil(length(fltr)/2), 1), Y1(i:i+gap-1), zeros(ceil(length(fltr)/2), 1));  %; fltr);%convolution
    %    c2 = cat(1, zeros(ceil(length(fltr)/2), 1), Y2(i:i+gap-1), zeros(ceil(length(fltr)/2), 1));  %; fltr);%convolution

    c1 = EOG1_lp_dev(i:i+gap-1);
    c2 = EOG2_lp_dev(i:i+gap-1);

    c1Filtered = EOG1_lphp_dev(i:i+gap-1);
    c2Filtered = EOG2_lphp_dev(i:i+gap-1);

    sigPiece1 = EOG1_lp(i:i+gap-1);
    sigPiece2 = EOG2_lp(i:i+gap-1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sprobuj wyrzucic ones, zamiast tego poprostu zeros

    rem1 = zeros(length(c1), 1); % -> length of convoluted signal from 1 epoch
    rem2 = zeros(length(c1), 1);

    rem1Filtered = zeros(length(c1), 1); % -> length of convoluted signal from 1 epoch
    rem2Filtered = zeros(length(c1), 1);


    %     %sprobuj wyrzucic ones, zamiast tego poprostu zeros
    %
    %     rem1 = ones(length(c1), 1)*0.01; % -> length of convoluted signal from 1 epoch
    %     rem2 = ones(length(c1), 1)*0.01;
    %
    %     rem1Filtered = ones(length(c1), 1)*0.01; % -> length of convoluted signal from 1 epoch
    %     rem2Filtered = ones(length(c1), 1)*0.01;


    %EOG
    % mark places in derivative vector, where relaxed steepnes threshold is exceeded
    rem1(find(c1>thrRS)') = thrRS;
    rem2(find(c2>thrRS)') = thrRS;
    rem1(find(c1<-thrRS)') = -thrRS;
    rem2(find(c2<-thrRS)') = -thrRS;

    % mark places in derivative vector, where basic steepnes threshold is exceeded
    rem1(find(c1>thrBS)') = thrBS;
    rem2(find(c2>thrBS)') = thrBS;
    rem1(find(c1<-thrBS)') = -thrBS;
    rem2(find(c2<-thrBS)') = -thrBS;


    %Filtered EOG
    % mark places in derivative vector, where relaxed steepnes threshold is exceeded
    rem1Filtered(find(c1Filtered>thrRS)') = thrRS;
    rem2Filtered(find(c2Filtered>thrRS)') = thrRS;
    rem1Filtered(find(c1Filtered<-thrRS)') = -thrRS;
    rem2Filtered(find(c2Filtered<-thrRS)') = -thrRS;

    % mark places in derivative vector, where basic steepnes threshold is exceeded
    rem1Filtered(find(c1Filtered>thrBS)') = thrBS;
    rem2Filtered(find(c2Filtered>thrBS)') = thrBS;
    rem1Filtered(find(c1Filtered<-thrBS)') = -thrBS;
    rem2Filtered(find(c2Filtered<-thrBS)') = -thrBS;




    plusminusREM1 = zeros(gap, 1);  % -> length of signal from 1 epoch
    plusminusREM2 = zeros(gap, 1);

    % mark if signal amplitude value is positive or negative
    plusminusREM1(find(sigPiece1>0)') = thrRS;
    plusminusREM1(find(sigPiece1<0)') = -thrRS;
    plusminusREM2(find(sigPiece2 > 0)') = thrRS;
    plusminusREM2(find(sigPiece2 < 0)') = -thrRS;

    remRes = zeros(length(c1), 1);  % -> length of convoluted signal from 1 epoch
    candidateP1 = zeros(length(c1), 1);  % -> length of convoluted signal from 1 epoch
    candidateP2 = zeros(length(c1), 1);  % -> length of convoluted signal from 1 epoch

    lastRemUsingEog2 = thrRange*remSearchGap;

    for jj = thrRange*remSearchGap:1:length(c1) - thrRange*remSearchGap

        %IF
        %1
        %basic threshold is crossed &
        %         or
        %2
        %relaxed threshold was crossed &
        %in 20msec before something was marked all the time (but what??) &
        %signal grows/decreases for at least 0.1 sec &
        %signal is on the same side as increase direction or signal
        %amplitude has to be low enough (below +/-100 uV)

        %3 we want to avoid the situation when channels in chosen place
        %both increase/decrease


        %%%%%%%%%%%%%%%%%%%  just for plots  %%%%%%%%%%%%%%%%%%%%%%
        if  (abs(rem1(jj)) == thrBS  || ...
                ( abs(rem1(jj)) == thrRS && checkRange(jj, c1, thrMin2) > (0.1*sampling) && ...
                ( abs(sigPiece1(jj)) < 100 || plusminusREM1(jj) == rem1(jj)) )) && abs(c2(jj) - rem1(jj)) > thrRS

            candidateP1(jj) = 1;

        end
        if (abs(rem2(jj)) == thrBS || ...
                ( abs(rem2(jj)) == thrRS && checkRange(jj, c2, thrMin2) > (0.1*sampling)  && ...
                ( abs(sigPiece2(jj)) < 100 || plusminusREM2(jj) == rem2(jj)))) && abs(c1(jj) - rem2(jj)) > thrRS

            candidateP2(jj) = 1;

        end
        %%%%%%%%%%%%%%%%%%%  just for plots  %%%%%%%%%%%%%%%%%%%%%%

        if  (abs(rem1(jj)) == thrBS  || ... % 1
                ( abs(rem1(jj)) == thrRS && checkRange(jj, c1, thrMin2) > (0.1*sampling) && ... % 2
                ( abs(sigPiece1(jj)) < 100 || plusminusREM1(jj) == rem1(jj)) )) && ... % 2 signal on the same side as direction of movement or amplitude < 100 uV
                abs(c2(jj) - rem1(jj)) > thrRS % 3


            if lastRemUsingEog2 > jj - remSearchGap

                remStart = lastRemUsingEog2 + 1;

            else

                remStart = jj - remSearchGap + 1;

            end

            if jj + remSearchGap < length(c1) - ceil(thrRange*remSearchGap)

                remEnd = jj + remSearchGap;

            else

                remEnd = length(c1) - ceil(thrRange*remSearchGap);

            end



            if epochOfInterest == count

                fprintf('eog1 cndt pnt jj: %d, remStart: %d\n', jj, remStart)

            end
            

            while remStart < remEnd


                if epochOfInterest == count
                    fprintf('$');
                end


                %first, remStart is not exactly jj (at the beginning of an epoch, it starts a bit earlier)

                %IF
                %1
                %basic threshold is crossed &
                %movement is in the opposite direction &   -> in addition
                %to condition in EOG1
                %                     or
                %2
                %relaxed threshold was crossed &
                %movement is in the opposite direction &   -> in addition
                %to condition in EOG1
                %in 20msec before something was marked all the time (but
                %what??) &
                %signal grows/decreases for at least 0.1 sec &
                %signal is on the same side as increase direction or signal
                %amplitude has to be low enough (below +/-100 uV)

                %3 we want to avoid the situation when channels in chosen
                %place both increase/decrease

                %4 direction must be opposite to eog1


                %look for value above threshold in 2 eog
                if ( abs(rem2(remStart)) == thrBS || ...  % 1
                        ( abs(rem2(remStart)) == thrRS && checkRange(remStart, c2, thrMin2) > (0.1*sampling)  && ... % 2
                        ( abs(sigPiece2(remStart)) < 100 || plusminusREM2(remStart) == rem2(remStart))) ) && ...   % 2
                        abs(c1(remStart) - rem2(remStart)) > thrRS && abs(rem2(remStart)) < abs(rem2(remStart)-rem1(jj)) % 3&4



                    %Basic threshold must be exceeded in at least 1 EOG
                    %channel
                    if abs(rem2(remStart)) + abs(rem1(jj)) < thrBS + thrRS

                        remStart = remStart + 1;
                        continue;

                    end


                    %REM starts in a place, which is identified faster
                    %in 2 EOGs and finishes in a place, where we are
                    %now
                    if jj < remStart
                        bg = jj;
                        eend = remStart;
                    else
                        bg = remStart;
                        eend = jj;
                    end


                    %if there is real REM, during it increase/decrease of the
                    %signal should be more abrupt than decreases/increases in
                    %the neighborhood

                    %whole REM area + 100 msec from each site
                    nbBg = bg-floor(0.1*sampling)+1;
                    nbEnd = eend+floor(0.1*sampling)-1;

                    if nbBg < 1
                        nbBg = 1;
                    end

                    if nbEnd > length(c1)
                        nbEnd = length(c1);
                    end

                    %neighbourhood of whole REM area + 100 msec from each site
                    rn1 = c1(nbBg:nbEnd); % rem neighborhood in DR1
                    rn2 = c2(nbBg:nbEnd); % rem neighborhood in DR2

                    %neighbourhood of places +/- 30 msec, where movement was accepted in each channel
                    lpn1 = c1(jj-floor(0.03*sampling)+1:jj+floor(0.03*sampling)-1); % local point 1 neighborhood (LPN)
                    lpn2 = c2(remStart-floor(0.03*sampling)+1:remStart+floor(0.03*sampling)-1); % local point 2 neighborhood (LPN)



                    if epochOfInterest == count
                        fprintf('Pair check jj: %d, remStart: %d\n', jj, remStart)
                    end



                    % if increase/decrease of signal recognized as REM
                    % is stronger than decrease/increase in the neighbourhood
                    if (rem1(jj) > 0 && clearNbrhThr*abs(min(rn1)) < max(lpn1) ) && ( rem2(remStart) < 0 && clearNbrhThr*max(rn2) < abs(min(lpn2)) ) ||...
                            (rem1(jj) < 0 && clearNbrhThr*max(rn1) < abs(min(lpn1)) ) && ( rem2(remStart) > 0 && clearNbrhThr*abs(min(rn2)) < max(lpn2))

                        %if everything fits till now, the probability
                        %of REM is very high. But one more problem are
                        %slow EOG movements which could make EOG a
                        %little faster than REM itself. To determine
                        %it, check values again in filtered signal

                        if abs(rem1Filtered(jj) - rem2Filtered(remStart)) >= thrRS + thrBS

                            if epochOfInterest == count

                                fprintf('Ngbrhd check jj: %d, remStart: %d, normal values: %f, filtered values: %f \n', jj, remStart, abs(rem1(jj) - rem2(remStart)), abs(rem1Filtered(jj) - rem2Filtered(remStart)))

                            end
                            %mark REM exactly in the middle
                            remRes(floor((remStart + jj)/2)) = 3;

                            %mark new place in EOG2, from which program
                            %will start to look for next REM in EOG1
                            lastRemUsingEog2 = remStart;

                            %stop searching
                            break;

                        elseif epochOfInterest == count

                            fprintf('not marked because of filtered sig: jj: %d = %f, remstart: %d = %f\n', jj, rem1Filtered(jj), remStart, rem2Filtered(remStart));

                        end

                    else

                        %if its not true, that is our candidate for
                        %REM turned out to have smaller steepnes
                        %than something next to it but directed in
                        %opposite direction, mark this place as REM
                        %search failure


                        % In which EOG do we have a problem?
                        %check EOG2
                        if (rem2(remStart) < 0 && max(rn2) >= abs(min(lpn2)) ) || ( rem2(remStart) > 0 && abs(min(rn2)) >= max(lpn2))

                            if count == epochOfInterest
                                fprintf('problem c2 in: %d, rem2(remStart): %f\n', remStart, rem2(remStart));
                            end
                        end

                        %check EOG1
                        if (rem1(jj) > 0 && abs(min(rn1)) >= max(lpn1) ) || (rem1(jj) < 0 && max(rn1) >= abs(min(lpn1)) )

                            if count == epochOfInterest
                                fprintf('problem c1 in: %d, rem1(jj): %f\n', jj, rem1(jj));
                            end

                            % if there is a problem in EOG 1 it makes
                            % no sense to look for REMs in jj, lets move forward
                            break;
                        end

                    end % if increase/decrease of signal recognized as REM
                    % is stronger than decrease/increase in the
                    % neighbourhood

                end % look for value above threshold in 2 eog

                remStart = remStart + 1;

            end %  while remStart < remEnd

        end

    end

    pairs = remRes;

    % plot results before consolidating
    if count == epochOfInterest

        subplot(3, 1, 3)
        plot(remRes, 'c','LineWidth',2);

    end

    %delete movements, which are:
    %1: very short (must exceed 30ms)
    
    %2 If REM is shorter than 70ms, it is considered as short. Short REMs should be clearly visible in
    % both EOG channels in order to minimize false detections. Therefore, we required for short REMs
    % that the median of values in the two derivative functions (left and right) should be above the basic
    % threshold. 

    strt = 0;
    for jj = 1:1:length(remRes)

        if remRes(jj) == 3

            strt = strt + 1;

        elseif strt && remRes(jj) ~= 3

            if strt < floor(0.04*sampling) || ... % 1
                    (strt < floor(0.07*sampling) && (abs(median(c1(jj-strt:jj - 1))) < 0.95*thrBS || abs(median(c2(jj-strt:jj- 1))) < 0.95*thrBS)) % 2

                remRes(jj-strt:jj- 1) = 0;
                if count == epochOfInterest
                    fprintf('deleting fragments that are too short %d:%d strt: %d, med1: %f, med2: %f\n', jj-strt,jj- 1,...
                        strt, median(c1(jj-strt:jj - 1)), median(c2(jj-strt:jj - 1)))
                end
            end

            strt = 0;

        end
    end

    %delete movements, which are too close to previous ones (0.11 sec)
    frag = floor(0.11*sampling);
    temp1 = find(remRes>0); %places where there are some REMs marked

    jj = 2;
    while jj <= length(temp1)

        if temp1(jj) > temp1(jj-1)+1  && temp1(jj) < temp1(jj-1) + frag  %if gap is too short...

            remRes(temp1(jj)) = 0;
            for jjj = jj+1:1:length(temp1)


                if temp1(jjj) == temp1(jjj-1)+1
                    remRes(temp1(jjj)) = 0;
                else

                    temp1 = find(remRes>0);
                    jj = 2;
                    break;

                end

            end

        end

        jj = jj+1;
    end



    %check founded REMs with eogClass

%    length(eogClass)
%    classCount
%    classCount + epLen-1
    
    
    while classCount + epLen-1 > length(eogClass)
        
        eogClass(length(eogClass) + 1) = -888;
        
    end
    
    fragm = eogClass(classCount:classCount + epLen-1);
    
    
    for jj = 1:1:epLen

        if fragm(jj) < 1

            remRes((jj-1)*sampling + 1:jj*sampling) = 0;

        end

    end


    if count == epochOfInterest
        measure = zeros(length(remRes), 1);
    end




    %exactlyMarkedREMs(i:i+gap-1) = remRes;

    REMNumPnt = 1; %remRes has the length of 1 epoch, exactlyMarkedREMs is a whole signal, 
                   %therefore pointer is used here
    for jj = i:1:i+gap-1

        if remRes(REMNumPnt) > 0
            exactlyMarkedREMs(jj) = 1;
        elseif remRes(REMNumPnt) <= 0
            exactlyMarkedREMs(jj) = 0;
        end

        REMNumPnt = REMNumPnt + 1;
    end


    %     exactlyMarkedREMs
    %     c1 = Y1(i:i+gap-1);
    %
    %
    %     REMNumPnt = 1;
    %     for jj = i:1:i+gap-1
    %
    %
    %         exactlyMarkedREMs(jj)
    %
    %     end



    if countAll

        %mark number of rem's for each epoch (Nb of REMs)
        REMNum = 0;
        previous = 0;
        for jj = 1:1:length(remRes)

            if ~previous && remRes(jj) > 0
                previous = previous+1;
                REMNum = REMNum+1;
            elseif previous && remRes(jj) <= 0
                previous = 0;
            end

        end

    else

        %mark number of rem's for each epoch (Nb of 3 sec fragments with REMs)
        REMNum = 0;

        for jj = 1:miniEp*sampling:length(remRes) - sampling

            %            jj
            %            length(remRes)

            if count == epochOfInterest
                measure(jj) = max(abs(EOG1_lp(i:i+gap-1)));
            end

            fragMiniEp = remRes(jj:jj+miniEp*sampling-1);

            if ~isempty(find(fragMiniEp > 0, 1))
                REMNum = REMNum+1;
            end

        end

    end


    %%%%%%%%%%%%%%%%%% statistics purposes %%%%%%%%%%%%%%%%%%%%%%%
    %meanREMSpeed = zeros(1,floor(sigLen/gap));
    %meanREMWidth = zeros(1,floor(sigLen/gap));

    meanREMWidthTemp = [];
    meanREMSpeedTemp = [];
    previous = 0;

    for jj = 1:1:length(remRes)
        if remRes(jj) > 0

            previous = previous+1;

        elseif previous && remRes(jj) <= 0
            %            previous
            meanREMWidthTemp(length(meanREMWidthTemp) + 1) = previous;
            meanREMSpeedTemp(length(meanREMSpeedTemp) + 1) = mean(abs(c1(jj-previous:jj))) + mean(abs(c2(jj-previous:jj)));
            %            mean(abs(c1(jj-previous:jj)))

            previous = 0;

        end
    end


    if(~isempty(meanREMSpeedTemp))

        meanREMSpeed(count) = mean(meanREMSpeedTemp);
        meanREMWidth(count) = mean(meanREMWidthTemp);

    end
    %%%%%%%%%%%%%%%%%% statistics purposes %%%%%%%%%%%%%%%%%%%%%%%


    if count == epochOfInterest

        makeDetailedPlotSaveVectorsData(i, orgSampling, remRes, EOG1, EOG2, ...
            EOG1_lp, EOG2_lp, measure, EOG1_lp_dev, EOG2_lp_dev, REMNum, EOG1_lphp, EOG2_lphp, EOG1_lphp_dev, EOG2_lphp_dev, gap, candidateP1, candidateP2, pairs);

    end


    remDens(count) = REMNum;
    count = count + 1;
    classCount = classCount + epLen;


end

%toc


%%

function makeDetailedPlotSaveVectorsData(i, orgSampling, remRes, originalEog1, originalEog2, LowPassFltrdEog1, LowPassFltrdEog2, measure, Y1, Y2, REMNum, eog1Filtered, eog2Filtered, Y1Filtered, Y2Filtered, gap, candidateP1, candidateP2, pairs)
figure
subplot(5, 1, 1)
plot(remRes*max(LowPassFltrdEog1(i:i+gap-1))/4, 'r');
hold on
plot(LowPassFltrdEog1(i:i+gap-1), 'b')
plot(LowPassFltrdEog2(i:i+gap-1), 'g')
plot(measure, 'm')
title('EOG');
hold off

subplot(5, 1, 2)
plot(Y1(i:i+gap-1), 'b')
hold on
plot(Y2(i:i+gap-1), 'g')
title(strcat('REMs found: ', num2str(REMNum)))
plot(remRes, 'r');
hold off

% subplot(5, 1, 3)
% hold on
% plot(c1, 'b')
% plot(c2, 'g')
%
% hold off


subplot(4, 1, 3)
hold on
plot(eog1Filtered(i:i+gap-1), 'b')
plot(eog2Filtered(i:i+gap-1), 'g')
title('Highpass filtered EOG');
hold off

subplot(4, 1, 4)
hold on
plot(Y1Filtered(i:i+gap-1), 'b')
plot(Y2Filtered(i:i+gap-1), 'g')
hold off


figure
%subplot(3, 1, 1)

hold on
distance = orgSampling;
plot(LowPassFltrdEog1(i:i+gap-1) + distance, 'Color', [0 0 0])
plot(zeros(1, length(i:i+gap-1)), 'Color', [0 0 0])

plot(LowPassFltrdEog2(i:i+gap-1) - distance, 'Color', [0 0 0])
title('EOG');

hold off


reog = originalEog1(i:i+gap-1);
leog = originalEog2(i:i+gap-1);
rDR = Y1(i:i+gap-1);
lDR = Y2(i:i+gap-1);
lpAndHpFiltRDR = Y1Filtered(i:i+gap-1);
lpAndHpFiltLDR = Y2Filtered(i:i+gap-1);

lowpassFiltREOG = LowPassFltrdEog1(i:i+gap-1);
lowpassFiltLEOG = LowPassFltrdEog2(i:i+gap-1);
lowpassAndHighpassFiltREOG = eog1Filtered(i:i+gap-1);
lowpassAndHighpassFiltLEOG = eog2Filtered(i:i+gap-1);

reog = reog';
leog = leog';
rDR = rDR';
lDR = lDR';
remRes = remRes';

fid = fopen('/home/marek/myMatlabScripts/REM density scoring/paperSignalPicture.txt', 'w');

fprintf(fid, 'reog = c(');
fprintf(fid, '%d, ',reog(1:length(reog) - 1));
fprintf(fid, ' %d) \n', reog(length(reog)));

fprintf(fid, 'leog = c(');
fprintf(fid, '%d, ',leog(1:length(leog) - 1));
fprintf(fid, ' %d) \n', leog(length(leog)));

fprintf(fid, 'lpFiltReog = c(');
fprintf(fid, '%d, ',lowpassFiltREOG(1:length(lowpassFiltREOG) - 1));
fprintf(fid, ' %d) \n', lowpassFiltREOG(length(lowpassFiltREOG)));

fprintf(fid, 'lpFiltLeog = c(');
fprintf(fid, '%d, ',lowpassFiltLEOG(1:length(lowpassFiltLEOG) - 1));
fprintf(fid, ' %d) \n', lowpassFiltLEOG(length(lowpassFiltLEOG)));

fprintf(fid, 'lpAndHpFiltReog = c(');
fprintf(fid, '%d, ',lowpassAndHighpassFiltREOG(1:length(lowpassAndHighpassFiltREOG) - 1));
fprintf(fid, ' %d) \n', lowpassAndHighpassFiltREOG(length(lowpassAndHighpassFiltREOG)));

fprintf(fid, 'lpAndHpFiltLeog = c(');
fprintf(fid, '%d, ',lowpassAndHighpassFiltLEOG(1:length(lowpassAndHighpassFiltLEOG) - 1));
fprintf(fid, ' %d) \n', lowpassAndHighpassFiltLEOG(length(lowpassAndHighpassFiltLEOG)));

fprintf(fid, 'rDR = c(');
fprintf(fid, '%d, ',rDR(1:length(rDR) - 1));
fprintf(fid, ' %d) \n', rDR(length(rDR)));

fprintf(fid, 'lDR = c(');
fprintf(fid, '%d, ',lDR(1:length(lDR) - 1));
fprintf(fid, ' %d) \n', lDR(length(lDR)));

fprintf(fid, 'lpAndHpFiltRDR = c(');
fprintf(fid, '%d, ',lpAndHpFiltRDR(1:length(lpAndHpFiltRDR) - 1));
fprintf(fid, ' %d) \n', lpAndHpFiltRDR(length(lpAndHpFiltRDR)));

fprintf(fid, 'lpAndHpFiltLDR = c(');
fprintf(fid, '%d, ',lpAndHpFiltLDR(1:length(lpAndHpFiltLDR) - 1));
fprintf(fid, ' %d) \n', lpAndHpFiltLDR(length(lpAndHpFiltLDR)));

fprintf(fid, 'candidateREOG = c(');
fprintf(fid, '%d, ',candidateP1(1:length(candidateP1) - 1));
fprintf(fid, ' %d) \n', candidateP1(length(candidateP1)));

fprintf(fid, 'candidateLEOG = c(');
fprintf(fid, '%d, ',candidateP2(1:length(candidateP2) - 1));
fprintf(fid, ' %d) \n', candidateP2(length(candidateP2)));

fprintf(fid, 'pairs = c(');
fprintf(fid, '%d, ',pairs(1:length(pairs) - 1));
fprintf(fid, ' %d) \n', pairs(length(pairs)));


fprintf(fid, 'remRes = c(');
fprintf(fid, '%d, ',remRes(1:length(remRes) - 1));
fprintf(fid, ' %d) \n', remRes(length(remRes)));


fclose(fid);

pause




















