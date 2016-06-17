function [res_filters, res_band_channels, res_full_channels] = spisop_pow_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% determine (average) power of specific frequency bands
% Copyright Frederik D. Weber

DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);
ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
AVGoverChannels = getParam('AVGoverChannels',listOfParameters);
FrequencyBandsFileName = getParam('FrequencyBandsFileName',listOfParameters);

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

DoReReference = getParam('DoReReference',listOfCoreParameters);%either yes or no
RerefDefinitionsFileName = getParam('RerefDefinitionsFileName',listOfCoreParameters);
listOfRerefDefinitionFiles = {};
if strcmp(DoReReference,'yes')
    
    if exist([pathInputFolder filesep RerefDefinitionsFileName],'file') ~= 2
        error(['RerefChannelsFileName file ' [pathInputFolder filesep RerefDefinitionsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfRerefDefinitionFiles = read_mixed_csv([pathInputFolder filesep RerefDefinitionsFileName],',');
    if ~(all(size(listOfDatasetsPaths) == size(listOfRerefDefinitionFiles)))
        error('files or number of Datasetspaths RerefChannels are invalid or do not aggree')
    end
    
    for iDefFiles = 1:size(listOfRerefDefinitionFiles,1)
            
        if exist([pathInputFolder filesep listOfRerefDefinitionFiles{iDefFiles}],'file') ~= 2
            error(['The rereference definitions file listed in RerefDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' does not exist'])
        end
        
        try
           temp_table = readtable([pathInputFolder filesep listOfRerefDefinitionFiles{iDefFiles}],'FileType','text','Delimiter',',');
        catch err
           error(['The rereference definitions file listed in RerefDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' is not readable'])
        end
        if (size(temp_table,2) ~= 3) || (size(temp_table,1) <1)
            error(['The rereference definitions file listed in RerefDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' is not readable'])
        end
        temp_table = [];

    end
            
end

ApplyLinearDeviationMontage = getParam('ApplyLinearDeviationMontage',listOfCoreParameters);%either yes or no
DelimiterLinearDeviationMontage = ',';
LinearDeviationMontageDefinitionsFileName = getParam('LinearDeviationMontageDefinitionsFileName',listOfCoreParameters);
listOfLinearDeviationMontageFiles = {};
if strcmp(ApplyLinearDeviationMontage,'yes')
    if exist([pathInputFolder filesep LinearDeviationMontageDefinitionsFileName],'file') ~= 2
        error(['LinearDeviationMontagePathsFileName file ' [pathInputFolder filesep LinearDeviationMontageDefinitionsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfLinearDeviationMontageFiles = read_mixed_csv([pathInputFolder filesep LinearDeviationMontageDefinitionsFileName],',');
    if ~(all(size(listOfDatasetsPaths) == size(listOfLinearDeviationMontageFiles)))
        error('files or number of Datasetspaths LinearDeviationMontagePaths are invalid or do not aggree')
    end
    
     for iDefFiles = 1:size(listOfLinearDeviationMontageFiles,1)
        if exist([pathInputFolder filesep listOfLinearDeviationMontageFiles{iDefFiles}],'file') ~= 2
            error(['The linear deviations montage definitions file listed in LinearDeviationMontageDefinitions for dataset number ' num2str(iDefFiles)  ' does not exist'])
        end
        
         try
          readtable([pathInputFolder filesep listOfLinearDeviationMontageFiles{iDefFiles}],'FileType','text','Delimiter',',');
        catch err
           error(['The linear deviations montage definitions file listed in LinearDeviationMontageDefinitions for dataset number ' num2str(iDefFiles)  ' is not readable'])
        end
    end
    
end

AggregationOfDatasetOutputsOfDetections = getParam('AggregationOfDatasetOutputsOfDetections',listOfCoreParameters);%If the aggregation of datasetOutputfiles should be skipped either full or fast or no default full



WindowProportion = 1; %default
try
    WindowProportion = str2num(getParam('WindowProportion',listOfParameters));%The proportion of hanning window that is applied. e.g. 1 means 100% of hanning window applied and 0.5 means 50% of hanning window with symmetrically 25% of each segment tail (left and right) given a hanning shape 
catch err
    
end

core_cfg = [];
UseFTfiltfilt = getParam('UseFTfiltfilt',listOfCoreParameters);
core_cfg.use_ft_filtfilt = strcmp(UseFTfiltfilt,'yes');

core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);
core_cfg.precision = getParam('ft_cfg_precision',listOfCoreParameters);

core_cfg.dftfilter = getParam('ft_cfg_dftfilter',listOfCoreParameters);
core_cfg.dftfreq = str2num(getParam('ft_cfg_dftfreq',listOfCoreParameters));



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
%NENBW.hanning = 1.5;
%NENBW.hamming = 1.3628;

tic
memtic
fprintf('POW function initialized\n');
conseciDatas = 1:length(iDatas);
parfor conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    FrqOfSmplWishedPar = FrqOfSmplWished;
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    
    channelsOfInterest = listOfChannelsOfInterest(iData,:);
    channelsOfInterest = channelsOfInterest(~(cellfun(@isempty,channelsOfInterest)));
    signalMultiplicator = listOfDataSetSignalMultiplicator(iData);
    signalOffsetSamples = listOfDataSetOffsetSamples(iData);

    
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
    
    if (signalOffsetSamples ~= 0)
        signalOffsetSeconds = signalOffsetSamples/preDownsampleFreq;
        roiBegins = roiBegins + signalOffsetSamples;
        roiEnds = roiEnds + signalOffsetSamples;
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
    
    if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
        
        cfg = [];
        %cfg = core_cfg;
        cfg.feedback = core_cfg.feedback;
        cfg.precision = core_cfg.precision;
        cfg.roiBegins = roiBegins;
        cfg.roiEnds = roiEnds;
        cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
        cfg.feedback = core_cfg.feedback;
        cfg = ft_definetrial(cfg);
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        cfg.dataset = datasetsPath;
        cfg.channel = 'all';
        fprintf('dataset %i: preprocess and pre filter data\n',iData);
        data = ft_fw_preprocessing(cfg);
        
        
        if strcmp(DoReReference,'yes')
            
            fileRerefSettings = listOfRerefDefinitionFiles{iData};
            table_reref = readtable([pathInputFolder filesep fileRerefSettings],'FileType','text','Delimiter',',');
            
            for iReref = size(table_reref,1)
                
                
                %referenceChannels = {'A1', 'A2'};
                %toBeReferencedChannels = {'C*', 'F*'};
                %newImplicitRefChannelLabel = 'Cz';
                
                referenceChannels = strsplit(table_reref.referenceChannels{iReref});
                toBeReferencedChannels = strsplit(table_reref.toBeReferencedChannels{iReref});
                
                referenceChannels = referenceChannels(~(cellfun(@isempty,referenceChannels)));
                toBeReferencedChannels = toBeReferencedChannels(~(cellfun(@isempty,toBeReferencedChannels)));
                
                newImplicitRefChannelLabel = table_reref.newImplicitRefChannelLabel{iReref};
                %             referenceChannels = listOfRerefDefinitionFiles(iData,:);
                %             referenceChannels = referenceChannels(~(cellfun(@isempty,referenceChannels)));
                %             toBeReferencedChannels = listOfToBeReferencedChannels(iData,:);
                %             toBeReferencedChannels = toBeReferencedChannels(~(cellfun(@isempty,toBeReferencedChannels)));
                %             newImplicitRefChannelLabel = listOfNewImplicitRefChannelLabels{iData};
                
                fprintf('dataset %i: reref data\n',iData);
                cfg = [];
                cfg.reref       = DoReReference;
                
                if strcmp(IgnoreDataSetHeader,'no')
                    cfg.channel = ft_channelselection([referenceChannels toBeReferencedChannels], hdr.label);
                else
                    cfg.channel = cellstr([referenceChannels toBeReferencedChannels]');
                end
                cfg.implicitref = newImplicitRefChannelLabel;% the implicit (non-recorded) reference channel is added to the data representation
                cfg.refchannel     = referenceChannels;
                data_reref        = ft_fw_preprocessing(cfg,data);
                
                cfg = [];
                notChan = strcat('-',data_reref.label);
                cfg.channel = ['all'; notChan(:)];
                data = ft_selectdata(cfg, data);
                
                cfg = [];
                data = ft_appenddata(cfg, data, data_reref);
                
                data_reref = [];
            end
        end
        
        
        if strcmp(ApplyLinearDeviationMontage,'yes')
            fprintf('dataset %i: apply linear deviation montage to data\n',iData);
            linearDeviationMontageFile = listOfLinearDeviationMontageFiles{iData};
            %linearDeviationMontagePath = 'T:\\Freddy\\temp_bva\\EEG_data\\save\\LinearDeviation_for_Schlafaus_SpiSOP.txt';

            
            montageTable = dataset('File',[pathInputFolder filesep linearDeviationMontageFile],'Delimiter',DelimiterLinearDeviationMontage,'ReadVarNames',true,'ReadObsNames',true);
            
            montage = [];
            montage.labelorg = get(montageTable,'VarNames');
            montage.labelnew  = get(montageTable,'ObsNames')';
            montage.tra = double(montageTable);
            
            data = ft_apply_montage(data,montage,'keepunused','yes','inverse','no');
        end
        
    end
    
    cfg = [];
    cfg = core_cfg;

    
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
    
    if strcmp(IgnoreDataSetHeader,'no')
        if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
            cfg.channel = ft_channelselection(channelsOfInterest, data.label);
        else
            cfg.channel = ft_channelselection(channelsOfInterest, hdr.label);
        end
    else
        cfg.channel = cellstr(channelsOfInterest');
    end
    cfg.feedback = core_cfg.feedback;
    fprintf('dataset %i: preprocess and pre filter data\n',iData);
    
    if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
        data = ft_fw_preprocessing(cfg,data);
    else
        cfg.roiBegins = roiBegins;
        cfg.roiEnds = roiEnds;
        cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
        cfg.feedback = core_cfg.feedback;
        cfg = ft_definetrial(cfg);
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        cfg.dataset = datasetsPath;
        data = ft_fw_preprocessing(cfg);
    end
    
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
    end
    
    if (signalMultiplicator ~= 1)
        data = ft_fw_factorMultiplicationOnSignal(data,'trial',signalMultiplicator);
    end
    
    FrqOfSmpl = data.fsample;%data.hdr.Fs;%samples per second / Hz
    
    %ROI2
    epochLengthSamples = epochLength * FrqOfSmpl;
    [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
    
    if (signalOffsetSamples ~= 0)
        signalOffsetSamples_downsampled = floor(signalOffsetSeconds*FrqOfSmpl);
        roiBegins = roiBegins + signalOffsetSamples_downsampled;
        roiEnds = roiEnds + signalOffsetSamples_downsampled;
    end
    
    trlSampleBeginsAndEnds = [roiBegins roiEnds];
    
    trlSampleLengths = roiEnds - roiBegins + 1;
    
    sampleLengthsAcrossROIs = sum(trlSampleLengths);
    lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/FrqOfSmpl; % in seconds
    
    NconsecutiveROIs = size(data.trial,2);
    
    guaranteedROIsegmentCoverage = lengthsAcrossROIsSeconds - ((SegmentLength * (1 - SegmentOverlap)) * NconsecutiveROIs);
    
    
    cfg = [];
    cfg.length    = SegmentLength;%single number (in unit of time, typically seconds) of the required snippets
    cfg.overlap   = SegmentOverlap;%single number (between 0 and 1 (exclusive)) specifying the fraction of overlap between snippets (0 = no overlap)
    cfg.feedback = core_cfg.feedback;
    data = ft_redefinetrial(cfg,data);
    
    NSegments = length(data.time);
    NSamplesPerSegment = length(data.time{1,1});
    freqResolutionCalculation = (FrqOfSmpl/NSamplesPerSegment);
    ENBW = NaN;
    if strcmp(ft_power_cfg_taper,'hanning')
        %ENBW = NENBW.hanning * freqResolutionCalculation;
        %     elseif strcmp(ft_power_cfg_taper,'hamming')
        %        ENBW = NENBW.hamming * freqResolutionCalculation;
        %ft_power_cfg_taper = 'hamming'
        
        windowFunction = ft_power_cfg_taper;
        %windowFunctionValues =  window(windowFunction, NSamplesPerSegment);
        %plot(windowFunctionValues);
        %windowFunctionValues = windowFunctionValues ./ norm(windowFunctionValues);
        
        
        temp_windowFunctionValues = window(windowFunction, floor(WindowProportion*NSamplesPerSegment));
        %plot(temp_windowFunctionValues);

        windowFunctionValuesLeft = temp_windowFunctionValues(1:floor(end/2));
        windowFunctionValuesRight = temp_windowFunctionValues(floor(end/2)+1:end);
        windowFunctionValues = ones(NSamplesPerSegment,1);
        windowFunctionValues(1:length(windowFunctionValuesLeft)) = windowFunctionValuesLeft;
        windowFunctionValues(end-length(windowFunctionValuesRight)+1:end) = windowFunctionValuesRight;
        
        %plot(windowFunctionValues);
        S1 = sum(windowFunctionValues);
        S2 = sum(windowFunctionValues.^2);
        
        NENBW = NSamplesPerSegment*(S2/(S1^2));
        ENBW = NENBW * freqResolutionCalculation;
        
    end
    
    
    
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    %cfg.pad = 'maxperlen';
    cfg.foilim  = [minBandFreq maxBandFreq];
    %cfg.foi  = [minFreq:FreqStepSize:maxFreq];
    cfg.taper = ft_power_cfg_taper;%'hanning'
    if strcmp(cfg.taper,'hanning') %%&& (WindowProportion ~= 1)
        cfg.taper = 'hanning_proportion';
        cfg.tapervalues = windowFunctionValues;
    end
    %if (~(strcmp(ft_power_cfg_taper,'hanning') || strcmp(ft_power_cfg_taper,'hamming')) && strcmp(ft_power_cfg_taper,'dpss'))
    if (~strcmp(ft_power_cfg_taper,'hanning')) && strcmp(ft_power_cfg_taper,'dpss')
        cfg.tapsmofrq = ft_power_cfg_tapsmofrq;%0.1;
    end
    
    
    if strcmp(IgnoreDataSetHeader,'no')
        cfg.channel = ft_channelselection(channelsOfInterest, hdr.label);
    else
        cfg.channel = cellstr(channelsOfInterest');
    end
    cfg.keeptrials = 'yes';
    cfg.feedback = core_cfg.feedback;
    fprintf('dataset %i: filter %i to %i Hz\n',iData,minBandFreq,maxBandFreq);
    tfa = ft_freqanalysis(cfg,data);
    
    pFreq = tfa.freq;
    pPower = tfa.powspctrm;
    
    tfa = [];%clear
    
    
    nChannels = length(data.label);
    %W = (trlSampleLengths./sampleLengthsAcrossROIs); %Nx1
    
    band_ch_meanPowerSumOverSegments = [];
    band_ch_meanPowerMeanOverSegments = [];
    
    
    for iBand = 1:(size(listOfFrequencyBands,1))
        %iBand = 1;
        band_meanPower = [];
        
        bandName = listOfFrequencyBands{iBand,1};
        bandMinFreq = str2num(listOfFrequencyBands{iBand,2});
        bandMaxFreq = str2num(listOfFrequencyBands{iBand,3});
        bandfoiIndex = find((pFreq >= bandMinFreq) & (pFreq <= bandMaxFreq));
        
        
        for iChan = 1:nChannels
            %iChan = 1;
            fprintf('dataset %i: process band %i to %i Hz in channel %s\n',iData,bandMinFreq,bandMaxFreq,data.label{iChan});
            
            trl_meanPower = [];
            for iTr = 1:size(pPower,1)
                trl_meanPower(iTr,:) = pPower(iTr,iChan,bandfoiIndex);
            end
            %band_ch_meanPower{iBand,iChan} = sum(W .* mean(trl_meanPower,2))/sum(W);%weighted mean
            band_ch_meanPowerSumOverSegments{iBand,iChan} = sum(mean(trl_meanPower,2));%sum over meaned power in band of channel
            band_ch_meanPowerMeanOverSegments{iBand,iChan} = mean(mean(trl_meanPower,2));%mean over meaned power in band of channel
        end
        
    end
    
    fprintf('dataset %i: write results\n',iData);
    
    if PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff == 0
        usedFilterOrder_hp_preDS = 0;
        hp_preDS_hdm.Fs = NaN;
        hp_preDS_hdm.Astop = NaN;
        hp_preDS_hdm.Fstop = NaN;
        hp_preDS_hdm.F6dB = NaN;
        hp_preDS_hdm.F3dB = NaN;
        hp_preDS_hdm.TransitionWidth = NaN;
        hp_preDS_hdm.Fpass = NaN;
        hp_preDS_hdm.Apass = NaN;
    end
    
    if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned'))
        
        usedFilterOrder_hp_preDS = NaN;
        hp_preDS_hdm.Fs = preDownsampleFreq;
        hp_preDS_hdm.Astop = NaN;
        hp_preDS_hdm.Fstop = NaN;
        hp_preDS_hdm.F6dB = NaN;
        hp_preDS_hdm.F3dB = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff;
        hp_preDS_hdm.TransitionWidth = NaN;
        hp_preDS_hdm.Fpass = NaN;
        hp_preDS_hdm.Apass = NaN;
    end
    
    hp_f_type_detail = '';
    switch core_cfg.hpfilttype
        case 'but'
            hp_f_type_detail = 'IIR_Butterworth_ml_butter';
        case 'fir'
            hp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
        case 'FIRdesigned'
            hp_f_type_detail = 'FIR_equiripple_signal_toolbox';
        case 'IIRdesigned'
            hp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
    end
    
    fidff = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
    %write header
    fprintf(fidff,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'datasetnum','dataset',...
        'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB');
    %write content
    fprintf(fidff,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' '\n'],...
        iData,datasetsPath,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass);
    
    
    fclose(fidff);
    
    
    
    %open output files
    fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    %write header of ouptufiles
    fprintf(fidc,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        'datasetnum','band','channel','mean_band_of_mean_power_over_segments','mean_band_of_mean_powerDensity_over_segments','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
        'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',...
        'mean_band_of_arb_energy_over_segments','mean_band_of_arb_energy_over_segments_per_ROI_seconds','mean_band_of_mean_power_over_segments_times_ROI_seconds',...
        'mean_band_of_arb_energyDensity_over_segments','mean_band_of_arb_energyDensity_over_segments_per_ROI_seconds','mean_band_of_mean_powerDensity_over_segments_times_ROI_seconds',...
        'min_freq_Hz','max_freq_Hz','lengths_ROI_seconds');
    
    
    for iBand = 1:(size(listOfFrequencyBands,1))
        %iBand = 1;
        band_meanPower = [];
        
        bandName = listOfFrequencyBands{iBand,1};
        bandMinFreq = str2num(listOfFrequencyBands{iBand,2});
        bandMaxFreq = str2num(listOfFrequencyBands{iBand,3});
        
        for iChan = 1:nChannels
            
            ch = data.label{iChan};
            
            
            fprintf(fidc,'%i,',iData);
            fprintf(fidc,'%s,',bandName);
            fprintf(fidc,'%s,',ch);
            fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan});
            fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan}/ENBW);
            fprintf(fidc,'%f,',epochLength);
            fprintf(fidc,'%f,',SegmentLength);
            fprintf(fidc,'%f,',SegmentOverlap);
            fprintf(fidc,'%i,',NSegments);
            fprintf(fidc,'%i,',NconsecutiveROIs);
            fprintf(fidc,'%f,',guaranteedROIsegmentCoverage);
            fprintf(fidc,'%e,',freqResolutionCalculation);
            %fprintf(fidc,'%e,',FreqStepSize);
            fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan});
            fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan}/lengthsAcrossROIsSeconds);
            fprintf(fidc,'%e,',band_ch_meanPowerMeanOverSegments{iBand,iChan}*lengthsAcrossROIsSeconds);
            fprintf(fidc,'%e,',band_ch_meanPowerSumOverSegments{iBand,iChan}/ENBW);
            fprintf(fidc,'%e,',(band_ch_meanPowerSumOverSegments{iBand,iChan}/ENBW)/lengthsAcrossROIsSeconds);
            fprintf(fidc,'%e,',(band_ch_meanPowerMeanOverSegments{iBand,iChan}/ENBW)*lengthsAcrossROIsSeconds);
            fprintf(fidc,'%f,',bandMinFreq);
            fprintf(fidc,'%f,',bandMaxFreq);
            fprintf(fidc,'%f\n',lengthsAcrossROIsSeconds);
        end
        
    end
    fclose(fidc);
    
    %open output files
    fidf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    %write header of ouptufiles
    fprintf(fidf,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        'datasetnum','channel','freq','mean_power_over_segments','mean_powerDensity_over_segments','epoch_length_seconds','segment_length_seconds','segments_overlap_proportion',...
        'segment_count','consecutive_ROI_count','guaranteed_ROI_segment_coverage_seconds','frequency_true_resolution_calculation',...%'frequency_stepsize_output',
        'arb_energy_over_segments','arb_energy_over_segments_per_ROI_seconds','mean_power_over_segments_times_ROI_seconds',...
        'arb_energyDensity_over_segments','arb_energyDensity_over_segments_per_ROI_seconds','mean_powerDensity_over_segments_times_ROI_seconds',...
        'lengths_ROI_Seconds'...
        );
    
    for iChan = 1:nChannels
        ch = data.label{iChan};
        tPowSumed = squeeze(sum(pPower(:,iChan,:),1));
        tPowMeaned = squeeze(mean(pPower(:,iChan,:),1));
        for iFreq = 1:(length(pFreq))
            fprintf(fidf,'%i,',iData);
            fprintf(fidf,'%s,',ch);
            fprintf(fidf,'%f,',pFreq(iFreq));
            fprintf(fidf,'%e,',tPowMeaned(iFreq));
            fprintf(fidf,'%e,',tPowMeaned(iFreq)/ENBW);
            fprintf(fidf,'%f,',epochLength);
            fprintf(fidf,'%f,',SegmentLength);
            fprintf(fidf,'%f,',SegmentOverlap);
            fprintf(fidf,'%i,',NSegments);
            fprintf(fidf,'%i,',NconsecutiveROIs);
            fprintf(fidf,'%f,',guaranteedROIsegmentCoverage);
            fprintf(fidf,'%e,',freqResolutionCalculation);
            %fprintf(fidf,'%e,',FreqStepSize);
            fprintf(fidf,'%e,',tPowSumed(iFreq));
            fprintf(fidf,'%e,',tPowSumed(iFreq)/lengthsAcrossROIsSeconds);
            fprintf(fidf,'%e,',tPowMeaned(iFreq)*lengthsAcrossROIsSeconds);
            fprintf(fidf,'%e,',tPowSumed(iFreq)/ENBW);
            fprintf(fidf,'%e,',(tPowSumed(iFreq)/ENBW)/lengthsAcrossROIsSeconds);
            fprintf(fidf,'%e,',(tPowMeaned(iFreq)/ENBW)*lengthsAcrossROIsSeconds);
            fprintf(fidf,'%f\n',lengthsAcrossROIsSeconds);
        end
    end
    
    fclose(fidf);
    data = [];%clear
    pPower = [];%clear
    pFreq = [];%clear
end

%aggregate all results from datasets
temp_fidff_all = [];
temp_fidc_all = [];
temp_fidf_all = [];
delimiter = ',';

if ~strcmp(AggregationOfDatasetOutputsOfDetections,'no')
fprintf('Aggregate results of all datasets\n');

for iData = iDatas
    
    
    temp_fidff = dataset('File', [pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    temp_fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    temp_fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    
    if iData == iDatas(1)
        temp_fidff_all = temp_fidff;
        temp_fidc_all = temp_fidc;
        temp_fidf_all = temp_fidf;
    else
        temp_fidff_all = cat(1,temp_fidff_all,temp_fidff);
        temp_fidc_all = cat(1,temp_fidc_all,temp_fidc);
        temp_fidf_all = cat(1,temp_fidf_all,temp_fidf);
    end
    
end
export(temp_fidff_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
export(temp_fidc_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_band_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
export(temp_fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'pow_full_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);

end
res_filters = temp_fidff_all;
res_band_channels = temp_fidc_all;
res_full_channels = temp_fidf_all;

fprintf('POW function finished\n');
toc
memtoc
end