function [res_filters, res_band_channels, res_full_channels] = spisop_pow_fluct_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% determine (average) power of specific frequency bands
% Copyright Frederik D. Weber

DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);
ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
AVGoverChannels = getParam('AVGoverChannels',listOfParameters);
FrequencyBandsFileName = getParam('FrequencyBandsFileName',listOfParameters);%delete

CenterFrequenciesFileName = getParam('CenterFrequenciesFileName',listOfParameters);
if exist([pathInputFolder filesep CenterFrequenciesFileName],'file') ~= 2
    error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep CenterFrequenciesFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
listOfCenterFrequencies = dlmread([pathInputFolder filesep CenterFrequenciesFileName],',')';

preCenterFreqFilterTo_FpassLeft = str2num(getParam('preCenterFreqFilterTo_FpassLeft',listOfParameters)); % in Hz
postCenterFreqFilterTo_FpassRight = str2num(getParam('postCenterFreqFilterTo_FpassRight',listOfParameters)); % in Hz

FreqSteps = str2num(getParam('FreqSteps',listOfParameters)); % in Hz
TimeSteps = str2num(getParam('TimeSteps',listOfParameters)); % in s

ECGchannelNameFileName = getParam('ECGchannelNameFileName',listOfParameters);
if exist([pathInputFolder filesep CenterFrequenciesFileName],'file') ~= 2
    error(['ECGchannelNameFileName file ' [pathInputFolder filesep ECGchannelNameFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
listOfECGchannel = read_mixed_csv([pathInputFolder filesep ECGchannelNameFileName],',');

if exist([pathInputFolder filesep DataSetPathsFileName],'file') ~= 2
    error(['DataSetPathsFileName file ' [pathInputFolder filesep DataSetPathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
if ~strcmp(IgnoreDataSetHeader,'yes')
    if exist([pathInputFolder filesep DataSetHeaderPathsFileName],'file') ~= 2
        error(['DataSetHeaderPathsFileName file ' [pathInputFolder filesep DataSetHeaderPathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
end
if exist([pathInputFolder filesep HypnogramsFileName],'file') ~= 2
    error(['HypnogramsFileName file ' [pathInputFolder filesep HypnogramsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
if exist([pathInputFolder filesep ChannelsOfInterestFileName],'file') ~= 2
    error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep ChannelsOfInterestFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
if exist([pathInputFolder filesep FrequencyBandsFileName],'file') ~= 2
    error(['FrequencyBandsFileName file ' [pathInputFolder filesep FrequencyBandsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz
%MaxPreDetectionFilterFreq = str2num(getParam('MaxPreDetectionFilterFreq',listOfParameters));%in Hz

%FreqStepSize = str2num(getParam('FreqStepSize',listOfParameters));

epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
%sleepStagesOfInterst = {'S3','S4'};
%sleepStagesOfInterst = {'SWS','S2'};
sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));

SegmentLength = str2num(getParam('SegmentLength',listOfParameters));% in seconds
SegmentOverlap = str2num(getParam('SegmentOverlap',listOfParameters));% in fraction of Segment length between 0 and 1


FrqOfSmplWished = str2num(getParam('FrqOfSmplWished',listOfParameters));%samples per second / Hz

DataSetsWhich = getParam('DataSetsWhich',listOfParameters);%Datasets to be processed either all or subset if subset then DataSetsNumbers is used for selection default all
DataSetsNumbers = str2num(getParam('DataSetsNumbers',listOfParameters));%The line numbers of the Datasets to be processed if DataSetsWich parameter is set to subset


listOfDatasetsPaths = read_mixed_csv([pathInputFolder filesep DataSetPathsFileName],',');
listOfDatasetHeaderPaths = [];
if strcmp(IgnoreDataSetHeader,'no')
    listOfDatasetHeaderPaths = read_mixed_csv([pathInputFolder filesep DataSetHeaderPathsFileName],',');
    if ~(all(size(listOfDatasetsPaths) == size(listOfDatasetHeaderPaths)))
        error('files or number of Datasetspaths and Headerpaths are invalid or do not aggree')
    end
end
listOfHypnogramPaths = read_mixed_csv([pathInputFolder filesep HypnogramsFileName],',');
listOfChannelsOfInterest = read_mixed_csv([pathInputFolder filesep ChannelsOfInterestFileName],',');
listOfFrequencyBands = read_mixed_csv([pathInputFolder filesep FrequencyBandsFileName],',');


if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
    error('files or number of Datasetspaths Hypnogramsfiles ChannelsOfInterest are invalid or do not aggree')
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end


if SegmentLength > epochLength
    error(['Parameter SegmentLength ' num2str(SegmentLength) 's must not be greater than epoch length ' num2str(epochLength) ' s!'])
end

minBandFreq = min(str2num(strjoin(listOfFrequencyBands(:,2)')));
maxBandFreq = max(str2num(strjoin(listOfFrequencyBands(:,3)')));

if ~(maxBandFreq*3 < FrqOfSmplWished)
    error(['FrqOfSmplWished of ' num2str(FrqOfSmplWished) ' Hz must be MORE THAN three-fold (i.e. 3-fold) the maximal band frequency of ' num2str(maxBandFreq) ' Hz, even lower than requested by Nyquist-Shannon sample theorem, choose higer FrqOfSmplWished or exclude higher frequency bands!'])
end

if (PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff > minBandFreq/2)
    error(['PreDownSampleHighPassFilter_FpassLeft of ' num2str(PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff) ' Hz must be smaller than half of minimal band frequency of ' num2str(minBandFreq) ' Hz!'])
end

% if (MaxPreDetectionFilterFreq < maxBandFreq)
% 	error(['MaxPreDetectionFilterFreq of ' num2str(MaxPreDetectionFilterFreq) ' Hz must be at least the maximal band frequency of ' num2str(maxBandFreq) ' Hz!'])
% end

if SegmentLength < (1/minBandFreq)
    error(['Parameter SegmentLength of ' num2str(SegmentLength) ' s is to short to support minimum band frequency of ' num2str(minBandFreq) ' Hz, choose either mimimum band ferquency of ' num2str((1/SegmentLength)) ' Hz or adapt sequence SegmentLength to ' num2str((1/minBandFreq)) ' s!'])
end

pretestHeaderForPersistentSampleFrequencies(IgnoreDataSetHeader,iDatas,listOfDatasetHeaderPaths,listOfChannelsOfInterest,FrqOfSmplWished);

SignalMultiplicator = getParam('SignalMultiplicator',listOfCoreParameters);%factor that signals should be muliplicated with either a number or mixed. e.g. -1 means inverted. in case of mixed DataSetSignalMultiplicatorFileName is used. default 1 (nothing)
DataSetSignalMultiplicatorFileName = getParam('DataSetSignalMultiplicatorFileName',listOfCoreParameters);%Filename of file containing an muliplicatoion factor (for example -1 for inversion) applied to each signal per line for respective dataset

if (strcmp(SignalMultiplicator,'mixed'))
    if exist([pathInputFolder filesep DataSetSignalMultiplicatorFileName],'file') ~= 2
        error(['DataSetSignalMultiplicatorFileName file ' [pathInputFolder filesep DataSetSignalMultiplicatorFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfDataSetSignalMultiplicator = load([pathInputFolder filesep DataSetSignalMultiplicatorFileName]);
    if ~(all(size(listOfDataSetSignalMultiplicator) == size(listOfDatasetsPaths)))
        error('files or number of Datasetspaths and DataSetSignalMultiplicator are invalid or do not aggree')
    end
else 
    listOfDataSetSignalMultiplicator = repmat(str2num(SignalMultiplicator),length(listOfDatasetsPaths),1);
end

SignalMultiplicator_ECG = getParam('SignalMultiplicator_ECG',listOfParameters);%factor that signals should be muliplicated with either a number or mixed. e.g. -1 means inverted. in case of mixed DataSetSignalMultiplicatorFileName is used. default 1 (nothing)
DataSetSignalMultiplicator_ECG_FileName = getParam('DataSetSignalMultiplicator_ECG_FileName',listOfParameters);%Filename of file containing an muliplicatoion factor (for example -1 for inversion) applied to each signal per line for respective dataset

if (strcmp(SignalMultiplicator_ECG,'mixed'))
    if exist([pathInputFolder filesep DataSetSignalMultiplicator_ECG_FileName],'file') ~= 2
        error(['DataSetSignalMultiplicator_ECG_FileName file ' [pathInputFolder filesep DataSetSignalMultiplicatorFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfDataSetSignalMultiplicator_ECG = load([pathInputFolder filesep DataSetSignalMultiplicator_ECG_FileName]);
    if ~(all(size(listOfDataSetSignalMultiplicator_ECG) == size(listOfDatasetsPaths)))
        error('files or number of Datasetspaths and DataSetSignalMultiplicator_ECG are invalid or do not aggree')
    end
else 
    listOfDataSetSignalMultiplicator_ECG = repmat(str2num(SignalMultiplicator_ECG),length(listOfDatasetsPaths),1);
end

DataSetOffsetSamples = getParam('DataSetOffsetSamples',listOfCoreParameters);%offset in the data in samples of the original sampling frequency of the file either a constant for all datasets or mixed.
DataSetOffsetSamplesFileName = getParam('DataSetOffsetSamplesFileName',listOfCoreParameters);%Filename of file containing an offset factor (for example -1 for inversion) applied to each signal per line for respective dataset

if (strcmp(DataSetOffsetSamples,'mixed'))
    if exist([pathInputFolder filesep DataSetOffsetSamplesFileName],'file') ~= 2
        error(['DataOffsetSamplesFileName file ' [pathInputFolder filesep DataSetOffsetSamplesFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfDataSetOffsetSamples = load([pathInputFolder filesep DataSetOffsetSamplesFileName]);
    if ~(all(size(listOfDataSetOffsetSamples) == size(listOfDatasetsPaths)))
        error('files or number of Datasetspaths and DataSetOffsetSamples (in File) are invalid or do not aggree, also check for empty lines in corresponding files')
    end
else 
    listOfDataSetOffsetSamples = repmat(str2num(DataSetOffsetSamples),length(listOfDatasetsPaths),1);
end

core_cfg = [];
UseFTfiltfilt = getParam('UseFTfiltfilt',listOfCoreParameters);
core_cfg.use_ft_filtfilt = strcmp(UseFTfiltfilt,'yes');

core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);
core_cfg.precision     = getParam('ft_cfg_precision',listOfCoreParameters);

core_cfg.dftfilter     = getParam('ft_cfg_dftfilter',listOfCoreParameters);
core_cfg.dftfreq       = str2num(getParam('ft_cfg_dftfreq',listOfCoreParameters));



% core_cfg.bpfilttype    = getParam('ft_cfg_bpfilttype',listOfCoreParameters);
% core_cfg.bpfiltdir     = getParam('ft_cfg_bpfiltdir',listOfCoreParameters);
% core_cfg.bpinstabilityfix = getParam('ft_cfg_bpinstabilityfix',listOfCoreParameters);

% core_cfg.lpfilttype    = getParam('ft_cfg_lpfilttype',listOfCoreParameters);
% core_cfg.lpfiltdir     = getParam('ft_cfg_lpfiltdir',listOfCoreParameters);
% core_cfg.lpinstabilityfix = getParam('ft_cfg_lpinstabilityfix',listOfCoreParameters);

core_cfg.hpfilttype    = getParam('ft_cfg_hpfilttype',listOfCoreParameters);
core_cfg.hpfiltdir     = getParam('ft_cfg_hpfiltdir',listOfCoreParameters);
core_cfg.hpinstabilityfix = getParam('ft_cfg_hpinstabilityfix',listOfCoreParameters);

% if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
%     error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
% end

% if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
%     error(['low pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
% end

if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

% if ~strcmp(core_cfg.bpfilttype,'FIRdesigned')
%     error(['filter type for band pass not supported, only FIRdesigned allowed'])
% end

% if ~strcmp(core_cfg.lpfilttype,'FIRdesigned')
%     error(['filter type for low pass not supported, only FIRdesigned allowed'])
% end
%
% if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned') )
%     error(['filter type for band pass not supported, only FIRdesigned or IIRdesigned allowed'])
% end

Apass = str2num(getParam('Apass',listOfCoreParameters)); %Attenuation of bandpass ripples in db
AstopLeft = str2num(getParam('AstopLeft',listOfCoreParameters)); %Attenuation of left stop band (<FstopLeft) ripples in db
% AstopRight = str2num(getParam('AstopRight',listOfCoreParameters)); %Attenuation of right stop band (>FstopRight) ripples in db

% Apass_bp = Apass;
% AstopLeft_bp = AstopLeft;
% AstopRight_bp = AstopRight;

% Apass_lp = Apass;
% AstopRight_lp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;



% StopToPassTransitionWidth_bp = str2num(getParam('StopToPassTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
% PassToStopTransitionWidth_bp = str2num(getParam('PassToStopTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

% PassToStopTransitionWidth_lp = str2num(getParam('PassToStopTransitionWidth_lp',listOfCoreParameters)); %for low pass filter frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

% StopToPassTransitionWidth_hp = str2num(getParam('StopToPassTransitionWidth_hp',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2

StopToPassTransitionWidth_hp_predownsample = str2num(getParam('StopToPassTransitionWidth_hp_predownsample',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2



% UseFixedFilterOrder_bp = (getParam('UseFixedFilterOrder_bp',listOfCoreParameters));

% UseFixedFilterOrder_lp = (getParam('UseFixedFilterOrder_lp',listOfCoreParameters));

UseFixedFilterOrder_hp = (getParam('UseFixedFilterOrder_hp',listOfCoreParameters));

% UseTwoPassAttenuationCorrection_bp = (getParam('UseTwoPassAttenuationCorrection_bp',listOfCoreParameters));

% UseTwoPassAttenuationCorrection_lp = (getParam('UseTwoPassAttenuationCorrection_lp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_hp = (getParam('UseTwoPassAttenuationCorrection_hp',listOfCoreParameters));

% FilterOrder_bp = str2num(getParam('FilterOrder_bp',listOfCoreParameters));

% FilterOrder_lp = str2num(getParam('FilterOrder_lp',listOfCoreParameters));

FilterOrder_hp = str2num(getParam('FilterOrder_hp',listOfCoreParameters));

MaximizeFilterOrderIfFixedFilterOrderIsUsed = str2num(getParam('MaximizeFilterOrderIfFixedFilterOrderIsUsed',listOfCoreParameters));

% useTwoPassFiltering_bp = 'no';

% useTwoPassFiltering_lp = 'no';

useTwoPassFiltering_hp = 'no';

% if ~isempty(strfind(core_cfg.bpfiltdir,'two'))
%     useTwoPassFiltering_bp = 'yes';
% end
%
% if ~isempty(strfind(core_cfg.lpfiltdir,'two'))
%     useTwoPassFiltering_lp = 'yes';
% end

if ~isempty(strfind(core_cfg.hpfiltdir,'two'))
    useTwoPassFiltering_hp = 'yes';
end

%filtfilt --  The length of the input x must be more than three times the filter order (N) defined as max(length(b)-1,length(a)-1).
%For best results, make sure the sequence you are filtering has length at least three times the filter order and tapers to zero on both edges.i.e.:
minSignalLengthSamples = epochLength*FrqOfSmplWished;
maxFilterOrder = floor((minSignalLengthSamples) / 3) - 1;

% if strcmp(UseFixedFilterOrder_bp,'yes') && strcmp(useTwoPassFiltering_bp,'yes') && (FilterOrder_bp > maxFilterOrder)
% 	error(['filter order for band pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
% elseif (FilterOrder_bp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
% 	FilterOrder_bp = maxFilterOrder;
% end

% if strcmp(UseFixedFilterOrder_lp,'yes') && strcmp(useTwoPassFiltering_lp,'yes') && (FilterOrder_lp > maxFilterOrder)
% 	error(['filter order for low pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
% elseif (FilterOrder_lp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
% 	FilterOrder_lp = maxFilterOrder;
% end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end


% if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
%     error('band pass filter order must be an even number')
% end

% if strcmp(UseFixedFilterOrder_lp,'yes') && logical(mod(FilterOrder_lp,2))
%     error('low pass order must be an even number')
% end

if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
    error('high pass order must be an even number')
end


% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
% if strcmp(useTwoPassFiltering_bp,'yes') && strcmp(UseTwoPassAttenuationCorrection_bp,'yes')
%     Apass_bp = Apass_bp/2;
%     AstopLeft_bp = AstopLeft_bp/2;
%     AstopRight_bp = AstopRight_bp/2;
% end

% if strcmp(useTwoPassFiltering_lp,'yes') && strcmp(UseTwoPassAttenuationCorrection_lp,'yes')
%     Apass_lp = Apass_lp/2;
%     AstopRight_lp = AstopRight_lp/2;
% end

if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end

ft_power_cfg_taper = getParam('ft_power_cfg_taper',listOfParameters);
ft_power_cfg_tapsmofrq = str2num(getParam('ft_power_cfg_tapsmofrq',listOfParameters));

%if ~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'hamming') || strcmp(ft_power_cfg_taper,'dpss'))
if ~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'dpss'))
    %error(['Parameter ft_power_cfg_taper = ' ft_power_cfg_taper ' is not supported, for now only ft_power_cfg_taper = hanning or ft_power_cfg_taper = hamming or ft_power_cfg_taper = dpss are supported !'])
    error(['Parameter ft_power_cfg_taper = ' ft_power_cfg_taper ' is not supported, for now only ft_power_cfg_taper = hanning or ft_power_cfg_taper = dpss are supported !'])
end

NENBW = [];
NENBW.hanning = 1.5;
%NENBW.hamming = 1.3628;


    flucs_powers = {};
    flucs_powers_freqs = {};
    flucs_powers_norm = {};
    
    datas = {};
    datas_complete = {};
    
    datas_ECG =  {};
    datas_ECG_complete =  {};
    
    datas_ECG_hilbert = {};
    datas_ECG_hilbert_complete = {};
    
    hypns_complete = {};
    
    number_of_bouts = {};
    mean_used_of_bout_length = {};
    
    ECG_signal_concat = {};
    ECG_R_wave_threshold = {};
    
    InstHRsmin = {};
    InstHRsmin_change = {};
    InstHRsmin_SignalPeaksSamples = {};

    
    ECG_signal_concat_compete = {};
    ECG_R_wave_threshold_compete = {};
    
    InstHRsmin_compete = {};
    InstHRsmin_complete_change = {};
    InstHRsmin_SignalPeaksSamples_complete = {};
    
    power_signal_freq2 = {};
    power_signal_freq2_complete = {};
    
    power_signal_freq2_fig = {};
    power_signal_freq2_complete_fig = {};
    
    power_signal_relative_fig = {};
    
    
tic
memtic
fprintf('POW fluctuation function initialized\n');
conseciDatas = 1:length(iDatas);%conseciDatas = 13:length(iDatas);
for conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    FrqOfSmplWishedPar = FrqOfSmplWished;
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    
    channelsOfInterest = listOfChannelsOfInterest(iData,:);
    channelsOfInterest = channelsOfInterest(~(cellfun(@isempty,channelsOfInterest)));
    signalMultiplicator = listOfDataSetSignalMultiplicator(iData);
    signalMultiplicator_ECG = listOfDataSetSignalMultiplicator_ECG(iData);
    signalOffsetSamples = listOfDataSetOffsetSamples(iData);


    centerFreqFilter = listOfCenterFrequencies(iData);
   
    
%     if minFreq < MinDetectionFrequency_FpassLeft
%         minFreq = MinDetectionFrequency_FpassLeft;
%     end
%     
%     if maxFreq > MaxDetectionFrequency_FpassRight
%         maxFreq = MaxDetectionFrequency_FpassRight;
%     end
    ECGchannel = listOfECGchannel(iData);
    ECGchannel = ECGchannel(~(cellfun(@isempty,ECGchannel)));

    
    
    hdr = [];
    preDownsampleFreq = 0;
    if strcmp(IgnoreDataSetHeader,'no')
        headerPath = listOfDatasetHeaderPaths{iData};
        hdr = ft_read_header(headerPath);
        if (FrqOfSmplWishedPar > hdr.Fs)
            warning(['dataset ' num2str(iData) ': designated frequency not supported by data, will use: ' num2str(hdr.Fs) ' Hz instead!']);
            FrqOfSmplWishedPar = hdr.Fs;
        end
        preDownsampleFreq = hdr.Fs;
    elseif strcmp(IgnoreDataSetHeader,'yes')
        fprintf('dataset %i: try to ignore header of data file\n',iData);
        cfg = [];
        cfg.roiBegins = [1];
        cfg.roiEnds = [100];
        cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
        cfg.feedback = core_cfg.feedback;
        cfg = ft_definetrial(cfg);
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        cfg.dataset = datasetsPath;
        cfg.channel = 1;
        cfg.feedback = core_cfg.feedback;
        tempdata = ft_fw_preprocessing(cfg);
        preDownsampleFreq = tempdata.fsample;
        tempdata = [];%clear
    else
        error('wrong parameter for IgnoreDataSetHeader either yes or no');
    end
    
    fprintf('dataset %i: process ROI from hypnogram info\n',iData);
    
    %ROI
    epochLengthSamples = epochLength * preDownsampleFreq;
    [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
    
    if length(roiEnds) < 1
        error(['no ROI in data left for analysis']);
    end
    
    indexLastIncludedROIinData = length(roiBegins);
    nSampleLength = -1;
    if strcmp(IgnoreDataSetHeader,'no')
        nSampleLength = hdr.nSamples*hdr.nTrials + hdr.nSamplesPre;
        if (roiEnds(end) > nSampleLength)
            nMissingSamples  = roiEnds(end) - nSampleLength;
            warning ([ num2str(nMissingSamples) ' scored (hypnogram) samples (' num2str(nMissingSamples/preDownsampleFreq) ' seconds) missing in the dataset ' num2str(iData)]);
            
            indexLastIncludedROIinData = find(roiBegins < nSampleLength,1,'last');
            roiBegins = roiBegins(1:indexLastIncludedROIinData);
            roiEnds = roiEnds(1:indexLastIncludedROIinData);
            
            if (roiEnds(end) > nSampleLength)
                roiEnds(indexLastIncludedROIinData) = nSampleLength;
            end;
            
        end
    end
    
    if length(roiEnds) < 1
        error(['no ROI in data left for analysis']);
    end
    
    if (signalOffsetSamples ~= 0)
        signalOffsetSeconds = signalOffsetSamples/preDownsampleFreq;
        roiBegins = roiBegins + signalOffsetSamples;
        roiEnds = roiEnds + signalOffsetSamples;
    end
    
    cfg = [];
    cfg = core_cfg;
    cfg.roiBegins = roiBegins;
    cfg.roiEnds = roiEnds;
    cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
    cfg.feedback = core_cfg.feedback;
    cfg = ft_definetrial(cfg);
    
    cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
    
    FpassLeft = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff; %left pass frequency in Hz
    FstopLeft = FpassLeft - StopToPassTransitionWidth_hp_predownsample; %left stop frequency in Hz
    usedFilterOrder_hp_preDS = NaN;
    hp_preDS_hdm = NaN;
    if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff ~= 0
        cfg.hpfilter = 'yes';

        if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned')
            hp_preDS_d = [];
            hp_preDS_hd = [];
            if strcmp(UseFixedFilterOrder_hp,'yes')
                hp_preDS_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,preDownsampleFreq);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
            else
                hp_preDS_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,preDownsampleFreq);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
            end
            fprintf('dataset %i: designing high pass filter for pre downsampling filtering \n',iData);
            if strcmp(core_cfg.hpfilttype,'IIRdesigned')
                hp_preDS_hd = design(hp_preDS_d,'butter'); %isstable(hp_preDS_hd)
            elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
                hp_preDS_hd = design(hp_preDS_d,'equiripple','MinOrder', 'even');
            else
                error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
            end
            usedFilterOrder_hp_preDS = hp_preDS_hd.order;
            cfg.hpfilterdesign = hp_preDS_hd;
            hp_preDS_hdm = measure(hp_preDS_hd);
        end
    else
        cfg.hpfilter = 'no';
    end
    if strcmp(UseFixedFilterOrder_hp,'yes')
        cfg.hpfiltord     = FilterOrder_hp;
    end
    cfg.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function
    
    cfg.dataset = datasetsPath;
    if strcmp(IgnoreDataSetHeader,'no')
        cfg.channel = ft_channelselection(channelsOfInterest, hdr.label);
    else
        cfg.channel = cellstr(channelsOfInterest');
    end
    fprintf('dataset %i: preprocess and pre filter data\n',iData);
    cfg.feedback = core_cfg.feedback;
    data = ft_fw_preprocessing(cfg);
    
    
    data_fsample_read_in = data.fsample;
    cfg.trl = [1 data_fsample_read_in*60*120 0];
    
    data_complete = ft_fw_preprocessing(cfg);

    
    
   
    
    
    
    
    if strcmp(AVGoverChannels,'yes')
        for iTr = 1:size(data.trial,2)
            data.trial{iTr} = mean(data.trial{iTr},1);
            
        end;
        data.label = {'meanOverChannels'};
    end
    
    
    %FrqOfSmplWished = 400;
    if (FrqOfSmplWishedPar < data.fsample)
        fprintf('dataset %i:resample data from %i to %i Hz\n',iData,data.fsample,FrqOfSmplWishedPar);
        cfg = [];
        cfg.resamplefs = FrqOfSmplWishedPar;%frequency at which the data will be resampled (default = 256 Hz)
        cfg.detrend = 'no';
        cfg.feedback = core_cfg.feedback;
        data = ft_resampledata(cfg,data);
        data_complete = ft_resampledata(cfg,data_complete);
    end
    
    if (signalMultiplicator ~= 1)
        data = ft_fw_factorMultiplicationOnSignal(data,'trial',signalMultiplicator);
        data_complete = ft_fw_factorMultiplicationOnSignal(data_complete,'trial',signalMultiplicator);
    end
    
    FrqOfSmpl = data.fsample;%data.hdr.Fs;%samples per second / Hz
    
    datas{iData} = data;
    datas_complete{iData} = data_complete;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%% ECG %%%
    
    
    
    
    
    cfg = [];
    cfg = core_cfg;
    cfg.roiBegins = roiBegins;
    cfg.roiEnds = roiEnds;
    cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
    cfg.feedback = core_cfg.feedback;
    cfg = ft_definetrial(cfg);
    
    cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
    
    FpassLeft = 20;%PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff; %left pass frequency in Hz
    FstopLeft = FpassLeft - 5; %left stop frequency in Hz
    
    usedFilterOrder_hp_preDS = NaN;
    hp_preDS_hdm = NaN;
    %if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff ~= 0
        cfg.hpfilter = 'yes';

        if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned')
            hp_preDS_d = [];
            hp_preDS_hd = [];
            if strcmp(UseFixedFilterOrder_hp,'yes')
                hp_preDS_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,preDownsampleFreq);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
            else
                hp_preDS_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,preDownsampleFreq);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
            end
            fprintf('dataset %i: designing high pass filter for pre downsampling filtering \n',iData);
            if strcmp(core_cfg.hpfilttype,'IIRdesigned')
                hp_preDS_hd = design(hp_preDS_d,'butter'); %isstable(hp_preDS_hd)
            elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
                hp_preDS_hd = design(hp_preDS_d,'equiripple','MinOrder', 'even');
            else
                error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
            end
            usedFilterOrder_hp_preDS = hp_preDS_hd.order;
            cfg.hpfilterdesign = hp_preDS_hd;
            hp_preDS_hdm = measure(hp_preDS_hd);
        end
    %else
    %    cfg.hpfilter = 'no';
    %end
    if strcmp(UseFixedFilterOrder_hp,'yes')
        cfg.hpfiltord     = FilterOrder_hp;
    end
    cfg.hpfreq = [FpassLeft];%dummy values are overwritten by low level function
    
    cfg.dataset = datasetsPath;
    if strcmp(IgnoreDataSetHeader,'no')
        cfg.channel = ft_channelselection(ECGchannel, hdr.label);
    else
        cfg.channel = cellstr(ECGchannel');
    end
    fprintf('dataset %i: preprocess and pre filter data\n',iData);
    cfg.feedback = core_cfg.feedback;
    data_ECG = ft_fw_preprocessing(cfg);
    
    cfg_trl = cfg.trl;
    
    cfg.trl = [1 data_fsample_read_in*60*120 0];
    data_ECG_complete = ft_fw_preprocessing(cfg);
    
    
   cfg.trl = cfg_trl;
    
    cfg.hilbert = 'abs';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 45;
    cfg.lpfiltord = 4;
    cfg.lpfilttype = 'but';
    %cfg.lpfilterdesign = lp_preDS_hd;
    data_ECG_hilbert = ft_fw_preprocessing(cfg);
    
    cfg.trl = [1 data_fsample_read_in*60*120 0];
    data_ECG_hilbert_complete = ft_fw_preprocessing(cfg);
    
    

    AVGoverChannels_ECG = 'no';
    if strcmp(AVGoverChannels_ECG,'yes')
        for iTr = 1:size(data_ECG.trial,2)
            data_ECG.trial{iTr} = mean(data_ECG.trial{iTr},1);
            
        end;
        data_ECG.label = {'meanOverECGChannels'};
    end
    
    
    %FrqOfSmplWishedPar_ECG = 400;
    FrqOfSmplWishedPar_ECG = 100;

    if (FrqOfSmplWishedPar_ECG < data_ECG.fsample)
        fprintf('dataset %i:resample data from %i to %i Hz\n',iData,data.fsample,FrqOfSmplWishedPar);
        cfg = [];
        cfg.resamplefs = FrqOfSmplWishedPar_ECG;%frequency at which the data will be resampled (default = 256 Hz)
        cfg.detrend = 'no';
        cfg.feedback = core_cfg.feedback;
        data_ECG = ft_resampledata(cfg,data_ECG);
        data_ECG_complete = ft_resampledata(cfg,data_ECG_complete);
        data_ECG_hilbert = ft_resampledata(cfg,data_ECG_hilbert); 
        data_ECG_hilbert_complete = ft_resampledata(cfg,data_ECG_hilbert_complete);
    end
    
    if (signalMultiplicator ~= 1)
        data_ECG = ft_fw_factorMultiplicationOnSignal(data_ECG,'trial',signalMultiplicator);
    end
    
    
    datas_ECG{iData} = data_ECG;
    datas_ECG_complete{iData} = data_ECG_complete;
    
    datas_ECG_hilbert{iData} = data_ECG_hilbert;
    datas_ECG_hilbert_complete{iData} = data_ECG_hilbert_complete;
    
    
%  
%     figure
%      xxx = abs(hilbert(data_ECG.trial{1}));
%      plot(xxx)
%     figure
%      xxx = data_ECG.trial{1};
%      plot(xxx)

% plot(xxx)
% plot(interpol_x_pow_HRmin)
%     [c lags] = xcorr(interpol_x_pow_HRmin(200:200+200*800),interpol_x_HRmin(200:200+200*800),200*200);
%     plot(lags/200,c)
%     plot((1:length(interpol_x_pow_HRmin))./200,interpol_x_pow_HRmin,'Color',[1 0 0])
    
    FrqOfSmpl_ECG = data_ECG.fsample;%data.hdr.Fs;%samples per second / Hz
    
    
    tempAllDataValues = [];
    tempAllDataValues_ECG = [];
    for iTr = 1:size(data_ECG_hilbert.trial,2)
        if iTr == 1
            tempAllDataValues = data_ECG_hilbert.trial{iTr};
            tempAllDataValues_ECG = data_ECG.trial{iTr};
        else
            tempAllDataValues = cat(2,tempAllDataValues,data_ECG_hilbert.trial{iTr});
            tempAllDataValues_ECG = cat(2,tempAllDataValues_ECG,data_ECG.trial{iTr});
        end
    end;
    
    %candSignal = data_ECG_hilbert.trial{1};

    candSignal = tempAllDataValues(:);
    candSignal_ECG = tempAllDataValues_ECG(:);
    threshold_ECG = 2*std(tempAllDataValues(:));
    
    ECG_signal_concat{iData} = candSignal;
    ECG_R_wave_threshold{iData} = threshold_ECG;
    
    
    minPeakDistanceSamples = FrqOfSmpl_ECG * 0.2;%200ms the refractory time of a heart beat
    [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(candSignal,'MINPEAKHEIGHT',threshold_ECG,'MINPEAKDISTANCE',minPeakDistanceSamples);
    candSignalPeaks = candSignal(tempCandSignalPeaksSamples);
    candSignalPeaksSamples = tempCandSignalPeaksSamples;%    candSignalPeaksSamples = currentRawDataSampleOffset + tempCandSignalPeaksSamples;

    
%     bin_length_seconds = 2.5;
%     bin_length_samples = FrqOfSmpl_ECG*bin_length_seconds;
%     [bincounts,ind] = histc(candSignalPeaksSamples,1:bin_length_samples:length(candSignal))
%     
%     bincounts
    samples_diff = diff(candSignalPeaksSamples);
    tempinstHRmin_pre = cat(1,0,60./(samples_diff(1:end-1)/FrqOfSmpl_ECG));
    tempinstHRmin_post = cat(1,0,60./(samples_diff(2:end)/FrqOfSmpl_ECG));
    tempinstHRmin_change = [0 ; tempinstHRmin_post-tempinstHRmin_pre];
%     plot([0 ; tempinstHRmin_post-tempinstHRmin_pre])
%     hold on
%     plot (tempinstHRmin(1:end))
    tempinstHRmin = cat(1,0,60./(samples_diff/FrqOfSmpl_ECG));
    nCandSignalPeaks = length(candSignalPeaksSamples);
    
%     figure;
%     plot(candSignal);
%     %plot(tempAllDataValues(:));
%     hold on
%     hline(threshold_ECG,'LineStyle','--','Color',[0 0 0]);
%     
%     for iENP = 1:length(candSignalPeaksSamples)
%         plot(candSignalPeaksSamples(iENP),candSignalPeaks(iENP),'o','MarkerFaceColor',[1 0 0])
%     end
%     
%     plot(candSignalPeaksSamples,tempinstHRmin,'Color',[0 1 0])
%     
%     interpol_x_HRmin = interp1(candSignalPeaksSamples,tempinstHRmin,1:length(candSignal));
%     %plot(candSignalPeaksSamples,tempinstHRmin,'o',1:length(candSignal),interpol_x_HRmin,':.');
%     
%     plot(candSignalPeaksSamples,tempinstHRmin,'Color',[0 1 0])
%     plot(candSignal_ECG,'Color',[0 1 1])
% 
% 
%     hold off

    InstHRsmin{iData} = tempinstHRmin;
    InstHRsmin_change{iData} = tempinstHRmin_change;
    InstHRsmin_SignalPeaksSamples{iData} = candSignalPeaksSamples;
    
tempAllDataValues_complete = [];
    tempAllDataValues_complete_ECG = [];
    
    for iTr = 1:size(data_ECG_hilbert_complete.trial,2)
        if iTr == 1
            tempAllDataValues_complete = data_ECG_hilbert_complete.trial{iTr};
            tempAllDataValues_complete_ECG = data_ECG_complete.trial{iTr};
        else
            tempAllDataValues_complete = cat(2,tempAllDataValues_complete,data_ECG_hilbert_complete.trial{iTr});
            tempAllDataValues_complete_ECG = cat(2,tempAllDataValues_complete_ECG,data_ECG_complete.trial{iTr});
        end
    end;
    
    %candSignal = data_ECG_hilbert.trial{1};

    candSignal = tempAllDataValues_complete(:);
    candSignal_ECG = tempAllDataValues_complete_ECG(:);
    threshold_ECG = 2*std(tempAllDataValues_complete(:));
    
    ECG_signal_concat_complete{iData} = candSignal;
    ECG_R_wave_threshold_complete{iData} = threshold_ECG;
    
    
    minPeakDistanceSamples = FrqOfSmpl_ECG * 0.2;%200ms the refractory time of a heart beat
    [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(candSignal,'MINPEAKHEIGHT',threshold_ECG,'MINPEAKDISTANCE',minPeakDistanceSamples);
    candSignalPeaks = candSignal(tempCandSignalPeaksSamples);
    candSignalPeaksSamples = tempCandSignalPeaksSamples;%    candSignalPeaksSamples = currentRawDataSampleOffset + tempCandSignalPeaksSamples;

    samples_diff = diff(candSignalPeaksSamples);
    tempinstHRmin_pre = cat(1,0,60./(samples_diff(1:end-1)/FrqOfSmpl_ECG));
    tempinstHRmin_post = cat(1,0,60./(samples_diff(2:end)/FrqOfSmpl_ECG));
    tempinstHRmin_change = [0 ; tempinstHRmin_post-tempinstHRmin_pre];
    
    tempinstHRmin_complete = cat(1,0,60./(diff(candSignalPeaksSamples)/FrqOfSmpl_ECG));
    nCandSignalPeaks = length(candSignalPeaksSamples);
    
    
%     figure;
%     plot(candSignal);
%     %plot(tempAllDataValues(:));
%     hold on
%     hline(threshold_ECG,'LineStyle','--','Color',[0 0 0]);
%     
%     for iENP = 1:length(candSignalPeaksSamples)
%         plot(candSignalPeaksSamples(iENP),candSignalPeaks(iENP),'o','MarkerFaceColor',[1 0 0])
%     end
%     
%     plot(candSignalPeaksSamples,tempinstHRmin_complete,'Color',[0 1 0])
%     
%     interpol_x_HRmin = interp1(candSignalPeaksSamples,tempinstHRmin_complete,1:length(candSignal));
%     %plot(candSignalPeaksSamples,tempinstHRmin,'o',1:length(candSignal),interpol_x_HRmin,':.');
%     
%     plot(candSignalPeaksSamples,tempinstHRmin_complete,'Color',[0 1 0])
%     plot(candSignal_ECG,'Color',[0 1 1])
% 
% 
%     hold off
    
    InstHRsmin_complete{iData} = tempinstHRmin_complete;
    InstHRsmin_complete_change{iData} = tempinstHRmin_change;
    InstHRsmin_SignalPeaksSamples_complete{iData} = candSignalPeaksSamples;
    
    
    
    %%%%%%%%% end ECG %%%%%%%%%%
    
    %%%%%%%%% Hypnogram plot data %%%%%%%%%

    
    
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    
    
    hypn_plot = hypn;
    hypn_plot(hypn_plot(:,1) == 8,1)
    hypn_plot(hypn_plot(:,1) == 4,1) = 3;
    hypn_plot(hypn_plot(:,1) == 5,1) = 4;
    hypn_plot(hypn_plot(:,1) == 8,1) = 0;
    hypn_plot(:,1) = hypn_plot(:,1)*-1;
    hypn_plot = hypn_plot(1:(120*2),1);
    hypn_plot_interpol = [];
    for iEp = 1:length(hypn_plot)
        hypn_plot_interpol = [hypn_plot_interpol; repmat(hypn_plot(iEp),30*FrqOfSmpl,1)];
    end
%     %plot(1:length(hypn_plot(:,1)),hypn_plot(:,1))
%      x = 1:size(hypn_plot,1);
%      
%      xq = 1:1/(30*(FrqOfSmpl)):size(hypn_plot,1);
%      v = hypn_plot(:,1);
%      
%     vq1 = interp1(x,v,xq,'cubic');
%     %plot(x,v,'o',xq,vq1)
%     %plot(vq1)
    
    hypns_complete{iData} = hypn_plot_interpol;
    
    %%%%%%%%% end Hypnogram plot data %%%%%%%%%
    
    %ROI2
    epochLengthSamples = epochLength * FrqOfSmpl;
    [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
    
    trlSampleBeginsAndEnds = [roiBegins roiEnds];
    
    trlSampleLengths = roiEnds - roiBegins + 1;
    
    sampleLengthsAcrossROIs = sum(trlSampleLengths);
    lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/FrqOfSmpl; % in seconds
    
    NconsecutiveROIs = size(data.trial,2);
    
    guaranteedROIsegmentCoverage = lengthsAcrossROIsSeconds - ((SegmentLength * (1 - SegmentOverlap)) * NconsecutiveROIs);
    
    
minFreq = 0.5;
maxFreq = 25;

    cfg = [];
    cfg.method    = 'wavelet';%single number (in unit of time, typically seconds) of the required snippets
    cfg.output   = 'pow';%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    cfg.foi = [minFreq:FreqSteps:maxFreq];%
    %cfg.foi = [10:0.25:15];

    cfg.width = 4;%7
    cfg.pad = 'maxperlen';
    cfg.feedback = core_cfg.feedback;
    cfg.keeptrials = 'yes';
    cfg.toi = [0:TimeSteps:max(cellfun(@max,data.time))];
    freq1 = ft_freqanalysis(cfg,data);
    cfg.toi = [0:TimeSteps:max(cellfun(@max,data_complete.time))];
    freq1_complete = ft_freqanalysis(cfg,data_complete);
    
    freqInc = 2;
    % 27 is 2Hz band 3 Hz below fast spindle peak
    % 28 is 2Hz band 3 Hz above fast spindle peak
    % 29 is fast spindle peak +- 1Hz (2Hz band)
    % 30 is SWA 0.5 to 4 Hz
    % 31 is sigma 10 to 15 Hz
    % 32 is theta 6 to 10 Hz
    % 33 is above sigma 16 to 20 Hz
    % 34 is beta 20 to 24 Hz
    % 35 is beta 4 to 8 Hz
    for iFreqMin = [27 28 29 30 31 32 33 34 45]
    %for iFreqMin = [1:1:24 27 28 29 30 31 32 33 34]

        if iFreqMin == 29
            minFreq_temp = centerFreqFilter - preCenterFreqFilterTo_FpassLeft;
            maxFreq_temp = centerFreqFilter + postCenterFreqFilterTo_FpassRight;
        elseif iFreqMin == 27
            minFreq_temp = centerFreqFilter - 5;
            maxFreq_temp = centerFreqFilter - 3;
        elseif iFreqMin == 28
            minFreq_temp = centerFreqFilter + 3;
            maxFreq_temp = centerFreqFilter + 5;
        elseif iFreqMin == 30
            minFreq_temp = 0.49;
            maxFreq_temp = 4.01;
        elseif iFreqMin == 31
            minFreq_temp = 9.99;
            maxFreq_temp = 15.01;
        elseif iFreqMin == 32
            minFreq_temp = 5.99;
            maxFreq_temp = 10.01;
        elseif iFreqMin == 33
            minFreq_temp = 15.99;
            maxFreq_temp = 20.01; 
        elseif iFreqMin == 34
            minFreq_temp = 19.99;
            maxFreq_temp = 24.01;
        elseif iFreqMin == 35
            minFreq_temp = 3.99;
            maxFreq_temp = 8.01;
        else
            minFreq_temp = iFreqMin;
            maxFreq_temp = iFreqMin+freqInc;
        end
    
    cfg_temp = [];
    
    
    cfg_temp.frequency  = [minFreq_temp maxFreq_temp];
    freq1_temp = ft_selectdata(cfg_temp,freq1);
    %size(freq1_temp.powspctrm)
    %size(freq1.powspctrm)
  
    

    freq2 = freq1_temp;
    freq2.powspctrm = squeeze(nanmean(freq1_temp.powspctrm,3));
    freq2.trial = {};
    freq2.time2 = {};
    freq2_fig = freq2;
    
    freq1_complete_temp = ft_selectdata(cfg_temp,freq1_complete);
    
    freq2_complete = freq1_complete_temp;
    freq2_complete.powspctrm = squeeze(nanmean(freq1_complete_temp.powspctrm,3));
    freq2_complete.trial = {};
    freq2_complete.time2 = {};
    freq2_complete_fig = freq2_complete;

    for iChn = 1:size(freq2.powspctrm,2)
        %iChn = 1
        %tempValuesToMean = freq2.powspctrm(:,iChn,:);%all parts
        %mice: normalize pwr to NonREM of the first 100 min after lights off
        %mice: normalize pwr to NonREM of the first 40 min after lights off
        tempValuesToMean = freq2.powspctrm(:,iChn,:);%all parts
        tempValuesToMean = tempValuesToMean(:);
        tempValuesToMean = tempValuesToMean(1:min(floor((1/TimeSteps)*60*100),length(tempValuesToMean(:))));
        temp_mean = nanmean(tempValuesToMean(:));
        for iRpt = 1:size(freq2.powspctrm,1)
            %iRpt = 1
            freq2_fig.powspctrm(iRpt,iChn,:) = freq2_fig.powspctrm(iRpt,iChn,:)./temp_mean;
            
            freq2.powspctrm(iRpt,iChn,:) = smooth(freq2.powspctrm(iRpt,iChn,:),4/TimeSteps,'moving');%smooth the signal to 4 seconds window
            freq2.powspctrm(iRpt,iChn,:) = ( freq2.powspctrm(iRpt,iChn,:) -  temp_mean) ./temp_mean;
        end
        
        freq2_complete_fig.powspctrm(iChn,:) = freq2_complete_fig.powspctrm(iChn,:)./temp_mean;
        
        freq2_complete.powspctrm(iChn,:) = smooth(freq2_complete.powspctrm(iChn,:),4/TimeSteps,'moving');%smooth the signal to 4 seconds window
        freq2_complete.powspctrm(iChn,:) = ( freq2_complete.powspctrm(iChn,:) -  temp_mean) ./temp_mean;
        
    end
    for iRpt = 1:size(freq2.powspctrm,1)
        %iRpt = 1;
        
        %size(freq2.powspctrm)
        freq2.trial{iRpt} = squeeze(freq2.powspctrm(iRpt,:,:));
        tempIndNan = isnan(freq2.trial{iRpt}(1,:));
        freq2.trial{iRpt}(:,tempIndNan) = [];
        freq2.time2{iRpt} = freq2.time(~tempIndNan);
        
        freq2_fig.trial{iRpt} = squeeze(freq2_fig.powspctrm(iRpt,:,:));
        tempIndNan = isnan(freq2_fig.trial{iRpt}(1,:));
        freq2_fig.trial{iRpt}(:,tempIndNan) = [];
        freq2_fig.time2{iRpt} = freq2_fig.time(~tempIndNan);
        
    end
    iRpt = 1;
        freq2_complete.trial{iRpt} = squeeze(freq2_complete.powspctrm(:,:));
        tempIndNan = isnan(freq2_complete.trial{iRpt}(1,:));
        freq2_complete.trial{iRpt}(:,tempIndNan) = [];
        freq2_complete.time2{iRpt} = freq2_complete.time(~tempIndNan);
        
        freq2_complete_fig.trial{iRpt} = squeeze(freq2_complete_fig.powspctrm(:,:));
        tempIndNan = isnan(freq2_complete_fig.trial{iRpt}(1,:));
        freq2_complete_fig.trial{iRpt}(:,tempIndNan) = [];
        freq2_complete_fig.time2{iRpt} = freq2_complete_fig.time(~tempIndNan);
        
        power_signal_freq2{iData,iFreqMin} = freq2;
        power_signal_freq2_complete{iData,iFreqMin} = freq2_complete;
        
        power_signal_freq2_fig{iData,iFreqMin} = freq2_fig;
        power_signal_freq2_complete_fig{iData,iFreqMin} = freq2_complete_fig;
        
     
    
        for iChanNum = 1:length(freq2.label)
                        
            temp_freq_plot = freq2_complete_fig;
            xxx = cat(2,temp_freq_plot.trial{:});
            xxx = xxx(iChanNum,:);
            %      figure
            %      plot(xxx)
            %      figure
            %      xxx = smooth(xxx,4/TimeSteps,'moving');
            %      plot(xxx)
            %
            %
            %
            x = 1:length(xxx);
            xq = 1:(1/(TimeSteps*FrqOfSmpl_ECG)):length(xxx);
            v = xxx;
            
            vq1 = interp1(x,v,xq);
            power_signal_relative_fig{iData,iFreqMin,iChanNum} = vq1;
            %     figure
            %     plot(vq1)
            %     hold on
            %     plot(InstHRsmin_SignalPeaksSamples_complete{iData},InstHRsmin_complete{iData},'Color',[0 1 0])
            %     plot(hypns_complete{iData}*10,'Color',[0 0 0])
            %
            %
          
            
        end
    
    %plot(x,v,'o',xq,vq1)
    %plot(vq1)
     
     
%      x = (1:length(xxx)).*(length(interpol_x_HRmin)/length(xxx));
%     interpol_x_pow_HRmin = interp1(x,xxx,1:length(interpol_x_HRmin));
% % plot(xxx)
% % plot(interpol_x_pow_HRmin)
% %     [c lags] = xcorr(interpol_x_pow_HRmin(200:200+200*800),interpol_x_HRmin(200:200+200*800),200*200);
% %     plot(lags/200,c)
%      plot((1:length(interpol_x_pow_HRmin))./200,interpol_x_pow_HRmin,'Color',[1 0 0])
%      hold on
%      plot((1:length(interpol_x_pow_HRmin))./200,interpol_x_HRmin,'Color',[0 1 0])

    freq2.time = freq2.time2;
    
    freq2 = rmfield(freq2,'dimord');
    freq2 = rmfield(freq2,'cumtapcnt');
    freq2 = rmfield(freq2,'freq');
    freq2 = rmfield(freq2,'powspctrm');
    freq2 = rmfield(freq2,'time2');
    
    freq2.fsample = 1/FreqSteps;
    
    
    
    
    freq2_fig.time = freq2_fig.time2;
    
    freq2_fig = rmfield(freq2_fig,'dimord');
    freq2_fig = rmfield(freq2_fig,'cumtapcnt');
    freq2_fig = rmfield(freq2_fig,'freq');
    freq2_fig = rmfield(freq2_fig,'powspctrm');
    freq2_fig = rmfield(freq2_fig,'time2');
    
    freq2_fig.fsample = 1/FreqSteps;

    
    
 
    
    
    
    
    tempStepSize = 0.5;
    cfg = [];
    cfg.method    = 'wavelet';%single number (in unit of time, typically seconds) of the required snippets
    cfg.output   = 'pow';%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    cfg.foi = [0.001:0.001:0.12];
    cfg.width = 4;
    cfg.pad = 'maxperlen';
    cfg.feedback = core_cfg.feedback;
    cfg.keeptrials = 'yes';
    cfg.toi = [0:tempStepSize:max(cellfun(@max,data.time))];
    freq3 = ft_freqanalysis(cfg,freq2);
    
    if length(freq3.label) > 1
        powspec = squeeze(nanmean(freq3.powspctrm,4));
    else
        powspec = (nanmean(freq3.powspctrm,4));
    end
    
    %     plot(freq3.freq,squeeze(powspec(1,1,:)))
    %     plot(freq3.freq,squeeze(powspec(2,1,:)))
    %     plot(freq3.freq,squeeze(powspec(3,1,:)))
    %     plot(freq3.freq,squeeze(powspec(4,1,:)))
    %     plot(freq3.freq,squeeze(powspec(5,1,:)))
    %     plot(freq3.freq,squeeze(powspec(6,1,:)))
    
    
    minSecBouts = 120;%4 30s epochs
    parfor iChanNum = 1:length(freq3.label)
        temp_exampleChannel =         char(strrep(freq3.label(iChanNum),'_','_'));
        titleName =  [ouputFilesPrefixString num2str(minFreq_temp) ' to ' num2str(maxFreq_temp) ' Hz Power fluctuations dataset ' num2str(iData) ' freq iter' num2str(iFreqMin) ' in ' temp_exampleChannel];
        
        
        %size(powspec)
        
        w = cellfun(@max,freq2.time);
        w2 = w(w>=minSecBouts);
        
        
        number_of_bouts{iData,iFreqMin,iChanNum} = length(w2);
        mean_used_of_bout_length{iData,iFreqMin,iChanNum} =  mean(w2);
        if length(w2) <2
            weighted_mean = squeeze(powspec(w>=minSecBouts,iChanNum,:))';
        else
            W = w2./sum(w2);
            W = W';
            A = squeeze(powspec(w>=minSecBouts,iChanNum,:));
            %A = squeeze(powspec(1:2,iChanNum,:));
            
            
            weighted_mean = zeros(1,size(A,2));
            for iFreq = 1:size(A,2)
                weighted_mean(iFreq) = nansum(W.*A(:,iFreq))/sum(W(~isnan(A(:,iFreq))));
            end
        end
        
        fh = figure;
        plot(freq3.freq,weighted_mean);
        
        title(titleName);
        xlabel('Frequncy (Hz)');
        ylabel('Power fluctuation relative to NonREM bouts');
        
        figure_width = 6;     % Width in inches
        figure_height = 5;    % Height in inches
        pos = get(fh, 'Position');
        set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
        % Here we preserve the size of the image when we save it.
        set(fh,'InvertHardcopy','on');
        set(fh,'PaperUnits', 'inches');
        papersize = get(fh, 'PaperSize');
        left = (papersize(1)- figure_width)/2;
        bottom = (papersize(2)- figure_height)/2;
        myfiguresize = [left, bottom, figure_width, figure_height];
        set(fh,'PaperPosition', myfiguresize);
        set(fh,'PaperOrientation', 'portrait');
        %print([titleName '.png'],'-dpng','-r600');
        saveas(fh, [titleName '.fig']);
        saveas(fh, [titleName '.png']);
        
        close(fh);
        
        
        flucs_powers{iData,iFreqMin,iChanNum} = weighted_mean;
        flucs_powers_norm{iData,iFreqMin,iChanNum} = weighted_mean./nanmean(weighted_mean);
        flucs_powers_freqs{iData,iFreqMin,iChanNum} = freq3.freq;
        %plot(flucs_powers_freqs{iData,iFreqMin},flucs_powers_norm{iData,iFreqMin})
        
    end
    
    end
    
    
%     if mod(conseciData, 20) == 0
%         iTemp = ceil(conseciData/2);
%         freq1 = [];
%     save([ouputFilesPrefixString 'process_fluct_d_16_all_temp' num2str(iTemp) '_freq.mat'],'-v7.3')
%     end
    %iData = 1;
    %figure;
    %plot(flucs_powers_freqs{iData},flucs_powers{iData})
    
    
    
    %     cfg = [];
    %     cfg.length    = SegmentLength;%single number (in unit of time, typically seconds) of the required snippets
    %     cfg.overlap   = SegmentOverlap;%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    %     cfg.feedback = core_cfg.feedback;
    %     data = ft_redefinetrial(cfg,data);
    %
    %     NSegments = length(data.time);
    %     NSamplesPerSegment = length(data.time{1,1});
    %     freqResolutionCalculation = (FrqOfSmpl/NSamplesPerSegment);
    %     ENBW = NaN;
    %     if strcmp(ft_power_cfg_taper,'hanning')
    %         ENBW = NENBW.hanning * freqResolutionCalculation;
    %         %     elseif strcmp(ft_power_cfg_taper,'hamming')
    %         %        ENBW = NENBW.hamming * freqResolutionCalculation;
    %
    %
%         
% %         NSamplesPerSegment = 1000;
% %         FrqOfSmpl = 100;
% %         freqResolutionCalculation = (FrqOfSmpl/NSamplesPerSegment);
% %         windowFunction = 'hanning';
% %         windowFunctionValues =  window(windowFunction, NSamplesPerSegment);
% %         plot(windowFunctionValues);
% %         
% %         %windowFunctionValues = windowFunctionValues ./ norm(windowFunctionValues);
% %         
% %         windowFunctionFactor = 0.1
% %         temp_windowFunctionValues = window(windowFunction, floor(2*windowFunctionFactor*NSamplesPerSegment));
% %         
% %         windowFunctionValuesLeft = temp_windowFunctionValues(1:floor(end/2))
% %         windowFunctionValuesRight = temp_windowFunctionValues(floor(end/2)+1:end)
% %         windowFunctionValues = ones(NSamplesPerSegment,1);
% %         windowFunctionValues(1:length(windowFunctionValuesLeft)) = windowFunctionValuesLeft;
% %         windowFunctionValues(end-length(windowFunctionValuesRight)+1:end) = windowFunctionValuesRight;
% %         
% %         plot(windowFunctionValues);
% %         S1 = sum(windowFunctionValues);
% %         S2 = sum(windowFunctionValues.^2);
% %         ENBW = NSamplesPerSegment*(S2/(S1^2))* freqResolutionCalculation;
%         
%     end
%     
%     
%     
%     cfg = [];
%     cfg.method = 'mtmfft';
%     cfg.output = 'pow';
%     %cfg.pad = 'maxperlen';
%     cfg.foilim  = [minBandFreq maxBandFreq];
%     %cfg.foi  = [minFreq:FreqStepSize:maxFreq];
%     cfg.taper = ft_power_cfg_taper;%'hanning';
%     
%     %if (~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'hamming')) && strcmp(ft_power_cfg_taper,'dpss'))
%     if (~strcmp(ft_power_cfg_taper,'hanning')) && strcmp(ft_power_cfg_taper,'dpss')
%         cfg.tapsmofrq = ft_power_cfg_tapsmofrq;%0.1;
%     end
%     
%     
%     if strcmp(IgnoreDataSetHeader,'no')
%         cfg.channel = ft_channelselection(channelsOfInterest, hdr.label);
%     else
%         cfg.channel = cellstr(channelsOfInterest');
%     end
%     cfg.keeptrials = 'yes';
%     cfg.feedback = core_cfg.feedback;
%     fprintf('dataset %i: filter %i to %i Hz\n',iData,minBandFreq,maxBandFreq);
%     tfa = ft_freqanalysis(cfg,data);
%     
%     pFreq = tfa.freq;
%     pPower = tfa.powspctrm;
%     
%     tfa = [];%clear
%     
%     
%     nChannels = length(data.label);
%     %W = (trlSampleLengths./sampleLengthsAcrossROIs); %Nx1
%     
%     band_ch_meanPowerSumOverSegments = [];
%     band_ch_meanPowerMeanOverSegments = [];
%     
%     
%     for iBand = 1:(size(listOfFrequencyBands,1))
%         %iBand = 1;
%         band_meanPower = [];
%         
%         bandName = listOfFrequencyBands{iBand,1};
%         bandMinFreq = str2num(listOfFrequencyBands{iBand,2});
%         bandMaxFreq = str2num(listOfFrequencyBands{iBand,3});
%         bandfoiIndex = find((pFreq >= bandMinFreq) & (pFreq <= bandMaxFreq));
%         
%         
%         for iChan = 1:nChannels
%             %iChan = 1;
%             fprintf('dataset %i: process band %i to %i Hz in channel %s\n',iData,bandMinFreq,bandMaxFreq,data.label{iChan});
%             
%             trl_meanPower = [];
%             for iTr = 1:size(pPower,1)
%                 trl_meanPower(iTr,:) = pPower(iTr,iChan,bandfoiIndex);
%             end
%             %band_ch_meanPower{iBand,iChan} = sum(W .* mean(trl_meanPower,2))/sum(W);%weighted mean
%             band_ch_meanPowerSumOverSegments{iBand,iChan} = sum(mean(trl_meanPower,2));%sum over meaned power in band of channel
%             band_ch_meanPowerMeanOverSegments{iBand,iChan} = mean(mean(trl_meanPower,2));%mean over meaned power in band of channel
%         end
%         
%     end
%     
%     fprintf('dataset %i: write results\n',iData);
%     
%     if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff == 0
%         usedFilterOrder_hp_preDS = 0;
%         hp_preDS_hdm.Fs = NaN;
%         hp_preDS_hdm.Astop = NaN;
%         hp_preDS_hdm.Fstop = NaN;
%         hp_preDS_hdm.F6dB = NaN;
%         hp_preDS_hdm.F3dB = NaN;
%         hp_preDS_hdm.TransitionWidth = NaN;
%         hp_preDS_hdm.Fpass = NaN;
%         hp_preDS_hdm.Apass = NaN;
%     end
%     
%     if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned'))
%         
%         usedFilterOrder_hp_preDS = NaN;
%         hp_preDS_hdm.Fs = preDownsampleFreq;
%         hp_preDS_hdm.Astop = NaN;
%         hp_preDS_hdm.Fstop = NaN;
%         hp_preDS_hdm.F6dB = NaN;
%         hp_preDS_hdm.F3dB = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff;
%         hp_preDS_hdm.TransitionWidth = NaN;
%         hp_preDS_hdm.Fpass = NaN;
%         hp_preDS_hdm.Apass = NaN;
%     end
%     
%     hp_f_type_detail = '';
%     switch core_cfg.hpfilttype
%         case 'but'
%             hp_f_type_detail = 'IIR_Butterworth_ml_butter';
%         case 'fir'
%             hp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
%         case 'FIRdesigned'
%             hp_f_type_detail = 'FIR_equiripple_signal_toolbox';
%         case 'IIRdesigned'
%             hp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
%     end
%     
%     fidff = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
%     %write header
%     fprintf(fidff,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
%         'datasetnum','dataset',...
%         'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB');
%     %write content
%     fprintf(fidff,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' '\n'],...
%         iData,datasetsPath,...
%         core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass);
%     
%     
%     fclose(fidff);
%     
%     
%     
%     %open output files
%     fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
%     
%     %write header of ouptufiles
%     fprintf(fidc,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
%         'datasetnum','band','channel','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
%         'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',...
%         'mean_band_of_sum_power_over_segments','mean_band_of_mean_power_over_segments','mean_band_of_sum_power_over_segments_per_ROI_seconds','mean_band_of_mean_power_over_segments_times_ROI_seconds',...
%         'mean_band_of_sum_powerDensity_over_segments','mean_band_of_mean_powerDensity_over_segments','mean_band_of_sum_powerDensity_over_segments_per_ROI_seconds','mean_band_of_mean_powerDensity_over_segments_times_ROI_seconds',...
%         'min_freq_Hz','max_freq_Hz','lengths_ROI_seconds');
%     
%     
%     for iBand = 1:(size(listOfFrequencyBands,1))
%         %iBand = 1;
%         band_meanPower = [];
%         
%         bandName = listOfFrequencyBands{iBand,1};
%         bandMinFreq = str2num(listOfFrequencyBands{iBand,2});
%         bandMaxFreq = str2num(listOfFrequencyBands{iBand,3});
%         
%         for iChan = 1:nChannels
%             
%             ch = data.label{iChan};
%             
%             
%             fprintf(fidc,'%i,',iData);
%             fprintf(fidc,'%s,',bandName);
%             fprintf(fidc,'%s,',ch);
%             fprintf(fidc,'%f,',epochLength);
%             fprintf(fidc,'%f,',SegmentLength);
%             fprintf(fidc,'%f,',SegmentOverlap);
%             fprintf(fidc,'%i,',NSegments);
%             fprintf(fidc,'%i,',NconsecutiveROIs);
%             fprintf(fidc,'%f,',guaranteedROIsegmentCoverage);
%             fprintf(fidc,'%e,',freqResolutionCalculation);
%             %fprintf(fidc,'%e,',FreqStepSize);
%             fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan});
%             fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan});
%             fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan}/lengthsAcrossROIsSeconds);
%             fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan}*lengthsAcrossROIsSeconds);
%             fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan}/ENBW);
%             fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan}/ENBW);
%             fprintf(fidc,'%e,',(band_ch_meanPowerSumOverSegments{iBand,iChan}/ENBW)/lengthsAcrossROIsSeconds);
%             fprintf(fidc,'%e,',(band_ch_meanPowerMeanOverSegments{iBand,iChan}/ENBW)*lengthsAcrossROIsSeconds);
%             fprintf(fidc,'%f,',bandMinFreq);
%             fprintf(fidc,'%f,',bandMaxFreq);
%             fprintf(fidc,'%f\n',lengthsAcrossROIsSeconds);
%         end
%         
%     end
%     fclose(fidc);
%     
%     %open output files
%     fidf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
%     
%     %write header of ouptufiles
%     fprintf(fidf,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
%         'datasetnum','channel','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
%         'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',
%         'freq',...
%         'sum_power_over_segments','mean_power_over_segments','sum_power_over_segments_per_ROI_seconds','mean_power_over_segments_times_ROI_seconds',...
%         'sum_powerDensity_over_segments','mean_powerDensity_over_segments','sum_powerDensity_over_segments_per_ROI_seconds','mean_powerDensity_over_segments_times_ROI_seconds',...
%         'lengths_ROI_Seconds'...
%         );
%     
%     for iChan = 1:nChannels
%         ch = data.label{iChan};
%         tPowSumed = squeeze(sum(pPower(:,iChan,:),1));
%         tPowMeaned = squeeze(mean(pPower(:,iChan,:),1));
%         for iFreq = 1:(length(pFreq))
%             fprintf(fidf,'%i,',iData);
%             fprintf(fidf,'%s,',ch);
%             fprintf(fidf,'%f,',epochLength);
%             fprintf(fidf,'%f,',SegmentLength);
%             fprintf(fidf,'%f,',SegmentOverlap);
%             fprintf(fidf,'%i,',NSegments);
%             fprintf(fidf,'%i,',NconsecutiveROIs);
%             fprintf(fidf,'%f,',guaranteedROIsegmentCoverage);
%             fprintf(fidf,'%e,',freqResolutionCalculation);
%             %fprintf(fidf,'%e,',FreqStepSize);
%             fprintf(fidf,'%f,',pFreq(iFreq));
%             fprintf(fidf,'%e,',tPowSumed(iFreq));
%             fprintf(fidf,'%e,',tPowMeaned(iFreq));
%             fprintf(fidf,'%e,',tPowSumed(iFreq)/lengthsAcrossROIsSeconds);
%             fprintf(fidf,'%e,',tPowMeaned(iFreq)*lengthsAcrossROIsSeconds);
%             fprintf(fidf,'%e,',tPowSumed(iFreq)/ENBW);
%             fprintf(fidf,'%e,',tPowMeaned(iFreq)/ENBW);
%             fprintf(fidf,'%e,',(tPowSumed(iFreq)/ENBW)/lengthsAcrossROIsSeconds);
%             fprintf(fidf,'%e,',(tPowMeaned(iFreq)/ENBW)*lengthsAcrossROIsSeconds);
%             fprintf(fidf,'%f\n',lengthsAcrossROIsSeconds);
%         end
%     end
%     
%     fclose(fidf);
%     data = [];%clear
%     pPower = [];%clear
%     pFreq = [];%clear
end

%iData = 3;
%plot(flucs_powers_freqs{iData},flucs_powers{iData})

%plot(flucs_powers_freqs{iData},flucs_powers_norm{iData})
save([ouputFilesPrefixString 'process_fluct7_all_freq.mat'],'-v7.3')
% 
% for conseciData = 1:27
%     %conseciData = 2
%     iData = iDatas(conseciData);
%     freqInc = 2;
%     % 29 is fast spindle peak
%     % 30 is SWA 0.5 to 4 Hz
%     for iFreqMin = [1:24 27 28 29 30 31]
%         %iFreqMin = 29
%         for iChanNum = 1:size(flucs_powers_freqs,3)
%             %iChanNum = 1
%             
%             %temp_exampleChannel =         char(strrep(freq2.label(iChanNum),'_','_'));
%             %titleName =  [num2str(minFreq_temp) ' to ' num2str(maxFreq_temp) ' in dataset ' num2str(iData) ' freq iter' num2str(iFreqMin) ' in channel number ' temp_exampleChannel];
%             titleName =  ['pow hr stage dataset ' num2str(iData) ' freq iter' num2str(iFreqMin) ' in channel number ' num2str(iChanNum)];
% 
%             fh = figure
%             
%             x_time = 1:length(power_signal_relative_fig{iData,iFreqMin,iChanNum});
%             x_time = x_time./(FrqOfSmplWished*60);
%             nPlots = 5;
%             hAxes1 = subplot(nPlots,1,1);
%             plot(x_time,smooth(power_signal_relative_fig{iData,iFreqMin,iChanNum}*100,FrqOfSmplWished*4,'moving'),'Color',[0 0 1])
%             ylabel('% of mean power');
%             xlabel('time in min');
%             title(titleName)
%             
%             x_time = InstHRsmin_SignalPeaksSamples_complete{iData}./(FrqOfSmplWished*60);
%             %y_smooth = smooth(abs(interp1(InstHRsmin_SignalPeaksSamples_complete{iData},InstHRsmin_complete_change{iData},1:(FrqOfSmplWished*60*120))),FrqOfSmplWished*4,'moving')
%             %x_time_smooth = 1:length(y_smooth);
%             %x_time_smooth = x_time_smooth./(FrqOfSmplWished*60);
%             hAxes2 = subplot(nPlots,1,2);
%             plot(x_time,InstHRsmin_complete{iData},'Color',[1 0 0])
%             ylabel('HR in beats per minute');
%             title('instantaneous HR')
%             hAxes3 = subplot(nPlots,1,3);
%             plot(x_time,abs(InstHRsmin_complete_change{iData}),'Color',[1 0 0])
%             %plot(x_time_smooth,y_smooth,'Color',[1 0 0])
%             ylabel('HR abs change in beats per minute');
%             title('instantaneous HR change')
%             
%             x_time = 1:length(datas_ECG_complete{iData}.trial{1});
%             x_time = x_time./(FrqOfSmplWished*60);
%             hAxes4 = subplot(nPlots,1,4);
%             plot(x_time,datas_ECG_complete{iData}.trial{1},'Color',[1 0 0])
%             ylabel('signal in �V');
%             title('ECG HR')
%             hAxes5 = subplot(nPlots,1,5);
%             plot(x_time,hypns_complete{iData},'Color',[0 0 0])
%             ylabel('stage [WAKE=0/S2=-2/SWS=-3/REM=-4]');
%             title('sleep stage')
%             
%             linkaxes([hAxes1,hAxes2,hAxes3,hAxes4,hAxes5], 'x');
%             
%             
%             figure_width = 6;     % Width in inches
%             figure_height = 8;    % Height in inches
%             pos = get(fh, 'Position');
%             set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
%             % Here we preserve the size of the image when we save it.
%             set(fh,'InvertHardcopy','on');
%             set(fh,'PaperUnits', 'inches');
%             papersize = get(fh, 'PaperSize');
%             left = (papersize(1)- figure_width)/2;
%             bottom = (papersize(2)- figure_height)/2;
%             myfiguresize = [left, bottom, figure_width, figure_height];
%             set(fh,'PaperPosition', myfiguresize);
%             set(fh,'PaperOrientation', 'portrait');
%             %print([titleName '.png'],'-dpng','-r600');
%             saveas(fh, [ titleName '.fig']);
%             %saveas(fh, ['new.fig']);
% 
%             %saveas(fh, [titleName '.png']);
%             
%             close(fh);
%         end
%     end
% end
% 
% fluct_peak_freqs = {};
% fluct_peaks = {};
% 
% freqInc = 2;
% % 29 is fast spindle peak
% % 30 is SWA 0.5 to 4 Hz
% for iFreqMin = [1:24 27 28 29 30 31]
%     %iFreqMin = 29;
%     for iChanNum = 1:size(flucs_powers_freqs,3)
%         %iChanNum = 1
%         %iData = 21
%         
%         
%         fh = figure;
%         val = [];
%         conseciDatas = 1:12
%         tempcolors = jet(length(conseciDatas));
%         for conseciData = conseciDatas
%             iData = iDatas(conseciData);
%             %iData = 2;
%             
%             temp_x_offset = 4;
%             curr_fluct_freqs =  flucs_powers_freqs{iData,iFreqMin,iChanNum}((1+temp_x_offset):end);
%             curr_fluct_signal = flucs_powers_norm{iData,iFreqMin,iChanNum}((1+temp_x_offset):end);
%             curr_fluct_freqs_plot =  flucs_powers_freqs{iData,iFreqMin,iChanNum}(1:end);
%             curr_fluct_signal_plot = flucs_powers_norm{iData,iFreqMin,iChanNum}(1:end);
%             
%             
%             f = fit(curr_fluct_freqs.',curr_fluct_signal.','gauss3');
%             
%             [fluct_peak, fluct_peak_sample] = findpeaks(f(curr_fluct_freqs),'MINPEAKHEIGHT',0,'NPEAKS',1,'SORTSTR','descend');
%             
%             if length(fluct_peak_sample) == 0
%                 fluct_peak_sample = 1;
%             end
%             fluct_peak_freqs{iData,iFreqMin,iChanNum} = curr_fluct_freqs(fluct_peak_sample);
%             fluct_peaks{iData,iFreqMin,iChanNum} = curr_fluct_signal(fluct_peak_sample);
% 
%             
%             plot(curr_fluct_freqs_plot,curr_fluct_signal_plot,'Color',tempcolors(conseciData,:));
%             hold on;
%             plot(curr_fluct_freqs_plot((1:length(curr_fluct_freqs))+temp_x_offset),f(curr_fluct_freqs),'LineStyle','.','Color',tempcolors(conseciData,:))
%             %plot(f,(1:length(curr_fluct_freqs))+temp_x_offset-1,curr_fluct_signal,'LineStyle','.','Color',tempcolors(conseciData,:))
%             %if length(fluct_peak_sample) > 0
%             plot(curr_fluct_freqs_plot(fluct_peak_sample+temp_x_offset),curr_fluct_signal(fluct_peak_sample),'o','MarkerFaceColor',tempcolors(conseciData,:))
%             vline(curr_fluct_freqs_plot(fluct_peak_sample+temp_x_offset),'LineStyle','--','Color',tempcolors(conseciData,:))
%             %end
%             val(iData,:,:) = curr_fluct_signal_plot;
%         end
%         
%         temp_peaks_freq = [fluct_peak_freqs{:,iFreqMin,iChanNum}];
%         temp_peaks_freq_mean = mean(temp_peaks_freq);
%         temp_peaks_freq_std = std(temp_peaks_freq);
%         
%         temp_peaks_freq_left = temp_peaks_freq_mean-temp_peaks_freq_std*0.5;
%         temp_peaks_freq_right = temp_peaks_freq_mean+temp_peaks_freq_std*0.5;
%         
%         temp_peaks = [fluct_peaks{:,iFreqMin,iChanNum}];
%         temp_peaks_mean = mean(temp_peaks);
%         temp_peaks_std = std(temp_peaks);
%         
%         
%         temp_peaks_mean2 = nanmean(nanmean(val(:,:,(curr_fluct_freqs_plot >= temp_peaks_freq_left) & (curr_fluct_freqs_plot <= temp_peaks_freq_right))))
%         
%         size(val)
%         size(squeeze(nanmean(val,1)))
%         %%figure;
%         y = squeeze(nanmean(val,1));
%         plot(curr_fluct_freqs_plot,y,'LineWidth',3,'Color','black');
%        vline(temp_peaks_freq_mean,'LineStyle','--','LineWidth',3,'Color','black')
%        hline(temp_peaks_mean,'LineStyle','--','LineWidth',3,'Color','black')
%               hline(temp_peaks_mean2,'LineStyle','--','LineWidth',1,'Color','black')
% 
%               
%         titleName =  ['Normalized Power fluctuations all dataset for freq iter' num2str(iFreqMin) ' in channel number ' num2str(iChanNum)];
%         
%         title(titleName);
%         xlabel('Frequncy (Hz)');
%         ylabel('normlized Power fluctuation relative to NonREM bouts');
%         
%         figure_width = 6;     % Width in inches
%         figure_height = 5;    % Height in inches
%         pos = get(fh, 'Position');
%         set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
%         % Here we preserve the size of the image when we save it.
%         set(fh,'InvertHardcopy','on');
%         set(fh,'PaperUnits', 'inches');
%         papersize = get(fh, 'PaperSize');
%         left = (papersize(1)- figure_width)/2;
%         bottom = (papersize(2)- figure_height)/2;
%         myfiguresize = [left, bottom, figure_width, figure_height];
%         set(fh,'PaperPosition', myfiguresize);
%         set(fh,'PaperOrientation', 'portrait');
%         %print([titleName '.png'],'-dpng','-r600');
%         saveas(fh, [titleName '.fig']);
%         saveas(fh, [titleName '.png']);
%         
%         close(fh);
%         
%         
% %         for conseciData = conseciDatas
% %             iData = iDatas(conseciData);
% %             curr_fluct_freqs =  flucs_powers_freqs{iData,iFreqMin,iChanNum}(3:end);
% %             curr_fluct_signal = flucs_powers_norm{iData,iFreqMin,iChanNum}(3:end);
% %             
% %             titleName =  [num2str(minFreq) ' to ' num2str(maxFreq) ' Hz Power fluctuations dataset ' num2str(iData) ' in channel number ' num2str(iChanNum)];
% %             
% %             [fluct_peak, fluct_peak_sample] = findpeaks(curr_fluct_signal,'MINPEAKHEIGHT',0,'NPEAKS',1,'SORTSTR','descend');
% %             
% %             plot(curr_fluct_freqs,curr_fluct_signal,'Color',tempcolors(iData,:));
% %             
% %             
% %             fh = figure;
% %             plot(curr_fluct_freqs,curr_fluct_signal);
% %             hold on;
% %             plot(curr_fluct_freqs(fluct_peak_sample),fluct_peak,'o','MarkerFaceColor',tempcolors(iData,:))
% %             vline(curr_fluct_freqs(fluct_peak_sample),'LineStyle','--','Color',tempcolors(iData,:))
% %             
% %             title(titleName);
% %             xlabel('Frequncy (Hz)');
% %             ylabel('Power fluctuation relative to NonREM bouts');
% %             
% %             figure_width = 6;     % Width in inches
% %             figure_height = 5;    % Height in inches
% %             pos = get(fh, 'Position');
% %             set(fh, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
% %             % Here we preserve the size of the image when we save it.
% %             set(fh,'InvertHardcopy','on');
% %             set(fh,'PaperUnits', 'inches');
% %             papersize = get(fh, 'PaperSize');
% %             left = (papersize(1)- figure_width)/2;
% %             bottom = (papersize(2)- figure_height)/2;
% %             myfiguresize = [left, bottom, figure_width, figure_height];
% %             set(fh,'PaperPosition', myfiguresize);
% %             set(fh,'PaperOrientation', 'portrait');
% %             %print([titleName '.png'],'-dpng','-r600');
% %             saveas(fh, [titleName '.fig']);
% %             saveas(fh, [titleName '.png']);
% %             
% %             close(fh);
% %             
% %         end
%     end
% end
% 
% % fprintf('Aggregate results of all datasets\n');
% % %aggregate all results from datasets
% % fidff_all = [];
% % fidc_all = [];
% % fidf_all = [];
% % delimiter = ',';
% % for iData = iDatas
% %     
% %     
% %     fidff = dataset('File', [pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
% %     fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
% %     fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
% %     
% %     if iData == iDatas(1)
% %         fidff_all = fidff;
% %         fidc_all = fidc;
% %         fidf_all = fidf;
% %     else
% %         fidff_all = cat(1,fidff_all,fidff);
% %         fidc_all = cat(1,fidc_all,fidc);
% %         fidf_all = cat(1,fidf_all,fidf);
% %     end
% %     
% % end
% % export(fidff_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
% % export(fidc_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
% % export(fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
% % 
% % res_filters = fidff_all;
% % res_band_channels = fidc_all;
% % res_full_channels = fidf_all;
% % 
% % fprintf('POW fluctuation function finished\n');
% % toc
% % memtoc