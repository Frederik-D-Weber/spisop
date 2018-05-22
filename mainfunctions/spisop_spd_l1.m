function [res_filters, res_channels, res_events, res_peaks, res_troughs, res_peaks_freq, res_troughs_freq, res_peaks_filtered_potential, res_troughs_filtered_potential] = spisop_spd_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% spindle and related phenomena detection
% Copyright Frederik D. Weber

functionName = 'spd';
ouputFilesPrefixString_folder = strtrim(ouputFilesPrefixString);
if strcmp(ouputFilesPrefixString_folder,'')
    ouputFilesPrefixString_folder = 'run0';
end
if ~isdir([pathOutputFolder filesep ouputFilesPrefixString_folder])
    mkdir([pathOutputFolder filesep ouputFilesPrefixString_folder]);
end
if ~isdir([pathOutputFolder filesep ouputFilesPrefixString_folder filesep functionName])
    mkdir([pathOutputFolder filesep ouputFilesPrefixString_folder filesep functionName]);
end
pathOutputFolder = [pathOutputFolder filesep ouputFilesPrefixString_folder filesep functionName];




DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);
ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
AVGoverChannels = getParam('AVGoverChannels',listOfParameters);
CenterFrequenciesFileName = getParam('CenterFrequenciesFileName',listOfParameters);

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

if exist([pathInputFolder filesep CenterFrequenciesFileName],'file') ~= 2
    error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep CenterFrequenciesFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end


PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz

preCenterFreqFilterTo_FpassLeft = str2num(getParam('preCenterFreqFilterTo_FpassLeft',listOfParameters)); % in Hz
postCenterFreqFilterTo_FpassRight = str2num(getParam('postCenterFreqFilterTo_FpassRight',listOfParameters)); % in Hz

MinDetectionFrequency_FpassLeft = str2num(getParam('MinDetectionFrequency_FpassLeft',listOfParameters));%in Hz
MaxDetectionFrequency_FpassRight = str2num(getParam('MaxDetectionFrequency_FpassRight',listOfParameters));%in Hz

MinDetectionLength = str2num(getParam('MinDetectionLength',listOfParameters));%in seconds
MaxDetectionLength = str2num(getParam('MaxDetectionLength',listOfParameters));%in seconds

MergeEventsInProximityWithinDetectionMargins = str2num(getParam('MergeEventsInProximityWithinDetectionMargins',listOfParameters));%in seconds

%consider RMS time window might be different from moving average time
%window!
RMSTimeWndw = str2num(getParam('RMSTimeWndw',listOfParameters)); %in seconds
MovAvgTimeWndw = str2num(getParam('MovAvgTimeWndw',listOfParameters)); %in seconds

EnvelopeMethod = getParam('EnvelopeMethod',listOfParameters);% either hilbertEnv or smoothedRMSwd
if ~(strcmp(EnvelopeMethod,'hilbertEnv') || strcmp(EnvelopeMethod,'smoothedRMSwd'))
    error(['EnvelopeMethod ' EnvelopeMethod ' in parameters is unknown, use either hilbertEnv or smoothedRMSwd'])
end


ThresholdFormationBasis = 'std'; %default
try
ThresholdFormationBasis = getParam('ThresholdFormationBasis',listOfParameters);% either mean or std default std
catch e
end
if ~(strcmp(ThresholdFormationBasis,'mean') || strcmp(ThresholdFormationBasis,'std'))
    error(['ThresholdFormationBasis ' ThresholdFormationBasis ' in parameters is unknown, use either mean or std, default std'])
end


ThresholdSignal = 'filtered_signal'; %default
try
ThresholdSignal = getParam('ThresholdSignal',listOfParameters);% either filtered_signal or envelope default signal
catch e
end

if ~(strcmp(ThresholdSignal,'filtered_signal') || strcmp(ThresholdSignal,'envelope'))
    error(['ThresholdSignal ' ThresholdSignal ' in parameters is unknown, use either filtered_signal or envelope, default filtered_signal'])
end




try
factorThresholdBeginEnd = str2num(getParam('factorThresholdBeginEnd',listOfParameters)); %factor in standard deviations for threshold of the signal for smoothed RMS signal
catch e
factorThresholdBeginEnd = str2num(getParam('factorSDbeginEnd',listOfParameters)); %factor in standard deviations for threshold of the signal for smoothed RMS signal
end

try
factorThresholdCriterion = str2num(getParam('factorThresholdCriterion',listOfParameters)); %factor in standard deviations for threshold of the signal for smoothed RMS signal
catch e
    factorThresholdCriterion = str2num(getParam('factorSDcriterion',listOfParameters)); %factor in standard deviations for threshold of the signal for smoothed RMS signal

end

if factorThresholdBeginEnd > factorThresholdCriterion
    error(['factorThresholdBeginEnd ' num2str(factorThresholdBeginEnd) ' in parameters must be less or equal than factorThresholdBeginEnd ' num2str(factorThresholdCriterion)])
end
    
MinAbsoluteDownToUpPeakPotential = str2num(getParam('MinAbsoluteDownToUpPeakPotential',listOfParameters));
MaxAbsoluteDownToUpPeakPotential = str2num(getParam('MaxAbsoluteDownToUpPeakPotential',listOfParameters));


epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
%sleepStagesOfInterst = {'S3','S4'};
%sleepStagesOfInterst = {'SWS','S2'};
sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));

try
    ThresholdAggregationMethod = getParam('ThresholdAggregationMethod',listOfParameters);% 'meanoverchan' or 'valuesoverchan' or 'respectivechan'
catch
    ThresholdAggregationMethod = getParam('SDmethod',listOfParameters);% 'meanoverchan' or 'valuesoverchan' or 'respectivechan'
end

UseAbsoluteEnvelopeThreshold = getParam('UseAbsoluteEnvelopeThreshold',listOfParameters);%If abosolute positive RMS potential threshold for all channels should be used either yes or no default no
AbsoluteEnvelopeThresholdBeginEnd = str2num(getParam('AbsoluteEnvelopeThresholdBeginEnd',listOfParameters));%Abosolute positive RMS potential threshold for all channels default 4 (for microVolts potential)
AbsoluteEnvelopeThresholdCriterion = str2num(getParam('AbsoluteEnvelopeThresholdCriterion',listOfParameters));%Abosolute positive RMS potential threshold for all channels default 4 (for microVolts potential)

if AbsoluteEnvelopeThresholdBeginEnd > AbsoluteEnvelopeThresholdCriterion
    error(['AbsoluteEnvelopeThresholdBeginEnd ' num2str(AbsoluteEnvelopeThresholdBeginEnd) ' in parameters must be less or equal than AbsoluteEnvelopeThresholdCriterion ' num2str(AbsoluteEnvelopeThresholdCriterion)])
end
    

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
listOfCenterFrequencies = dlmread([pathInputFolder filesep CenterFrequenciesFileName],',')';

if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && all(size(listOfDatasetsPaths) == size(listOfCenterFrequencies')) && (size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
    error('files or number of Datasetspaths Hypnogramsfiles ChannelsOfInterest and CenterFrequnecies are invalid or do not aggree')
end

for iCenterFreq = 1:length(listOfCenterFrequencies)
    centerFreq = listOfCenterFrequencies(iCenterFreq);
    if (centerFreq < MinDetectionFrequency_FpassLeft) || (centerFreq >= MaxDetectionFrequency_FpassRight)
        error('at least one of the CenterFrequnecies is not in the boundaries of MinDetectionFrequency_FpassLeft and MaxDetectionFrequency_FpassRight.')
    end
end


iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end

if strcmp(UseAbsoluteEnvelopeThreshold,'yes')
    factorThresholdBeginEnd = 1;
    AbsoluteEnvelopeThresholdCriterionRatio = AbsoluteEnvelopeThresholdCriterion/AbsoluteEnvelopeThresholdBeginEnd;
    factorThresholdCriterion = AbsoluteEnvelopeThresholdCriterionRatio;
end


if epochLength < (1/MinDetectionFrequency_FpassLeft)
    error(['Parameter epochLength ' num2str(epochLength) 's must not be greater in order to support the maximum of MinDetectionFrequency_FpassLeft of ' num2str(MinDetectionFrequency_FpassLeft) ' Hz!'])
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
        error('files or number of Datasetspaths and LinearDeviationMontagePaths are invalid or do not aggree')
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

OverwiteGlobalThresholdsAndUseIndividualDatasetEnvelopeTh = getParam('OverwiteGlobalThresholdsAndUseIndividualDatasetEnvelopeThreholds',listOfParameters);%determine if the Threshold Criterions for all datasets should be overwritten by the Thresholds stated in the file given in Parameter IndividualDatasetEnvelopeThreholdsFileName. either yes or no default no
IndividualDatasetEnvelopeThresholdsFileName = getParam('IndividualDatasetEnvelopeThresholdsFileName',listOfParameters);%Filename of file containing Threshold given in positive number pairs of ThresholdBeginEnd and ThresholdCriterion separated by one comma per line for respective dataset
DelimiterIndividualDatasetEnvelopeThrehold = ',';
listOfIndividualDatasetEnvelopeThrehold = {};
if strcmp(OverwiteGlobalThresholdsAndUseIndividualDatasetEnvelopeTh,'yes')
    if exist([pathInputFolder filesep IndividualDatasetEnvelopeThresholdsFileName],'file') ~= 2
        error(['IndividualDatasetEnvelopeThresholdsFileName file ' [pathInputFolder filesep IndividualDatasetEnvelopeThresholdsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    
    listOfIndividualDatasetEnvelopeThrehold = read_mixed_csv([pathInputFolder filesep IndividualDatasetEnvelopeThresholdsFileName],DelimiterIndividualDatasetEnvelopeThrehold);
    if ~(all(size(listOfDatasetsPaths,1) == size(listOfIndividualDatasetEnvelopeThrehold,1)))
        error('files or number of Datasetspaths and listOfIndividualDatasetEnvelopeThrehold are invalid or do not aggree')
    end
    
    for iThreholdParams = 1:size(listOfIndividualDatasetEnvelopeThrehold,1)
        if length(listOfIndividualDatasetEnvelopeThrehold(iThreholdParams,:)) ~= 2
            error(['The individual threshold definitions in the file given in IndividualDatasetEnvelopeThresholdsFileName are invalid for dataset number ' num2str(iThreholdParams)  ', the arguments must be 2 numbers separated by ' DelimiterIndividualDatasetEnvelopeThrehold])
        end
        if str2num(listOfIndividualDatasetEnvelopeThrehold{iThreholdParams,1}) > str2num(listOfIndividualDatasetEnvelopeThrehold{iThreholdParams,2})
            error(['The individual threshold definitions in the file given in IndividualDatasetEnvelopeThresholdsFileName are invalid for dataset number ' num2str(iThreholdParams)  ', the first number for begin-end criterion is bigger than the second for min reach criterion and this is not allowed.'])
         end
    end
end







core_cfg = [];
UseFTfiltfilt = getParam('UseFTfiltfilt',listOfCoreParameters);
core_cfg.use_ft_filtfilt = strcmp(UseFTfiltfilt,'yes');

core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);
core_cfg.precision     = getParam('ft_cfg_precision',listOfCoreParameters);

core_cfg.dftfilter     = getParam('ft_cfg_dftfilter',listOfCoreParameters);
core_cfg.dftfreq       = str2num(getParam('ft_cfg_dftfreq',listOfCoreParameters));



core_cfg.bpfilttype    = getParam('ft_cfg_bpfilttype',listOfCoreParameters);
core_cfg.bpfiltdir     = getParam('ft_cfg_bpfiltdir',listOfCoreParameters);
core_cfg.bpinstabilityfix = getParam('ft_cfg_bpinstabilityfix',listOfCoreParameters);

% core_cfg.lpfilttype    = getParam('ft_cfg_lpfilttype',listOfCoreParameters);
% core_cfg.lpfiltdir     = getParam('ft_cfg_lpfiltdir',listOfCoreParameters);
% core_cfg.lpinstabilityfix = getParam('ft_cfg_lpinstabilityfix',listOfCoreParameters);

core_cfg.hpfilttype    = getParam('ft_cfg_hpfilttype',listOfCoreParameters);
core_cfg.hpfiltdir     = getParam('ft_cfg_hpfiltdir',listOfCoreParameters);
core_cfg.hpinstabilityfix = getParam('ft_cfg_hpinstabilityfix',listOfCoreParameters);


if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
    error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

% if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
%     error(['low pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
% end

if (strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.hpinstabilityfix,'no'))
    error(['high pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end


% if ~strcmp(core_cfg.bpfilttype,'FIRdesigned')
%     error(['filter type for band pass not supported, only FIRdesigned allowed'])
% end
%
% if ~strcmp(core_cfg.lpfilttype,'FIRdesigned')
%     error(['filter type for low pass not supported, only FIRdesigned allowed'])
% end
%
% if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned') )
%     error(['filter type for band pass not supported, only FIRdesigned or IIRdesigned allowed'])
% end

Apass = str2num(getParam('Apass',listOfCoreParameters)); %Attenuation of bandpass ripples in db
AstopLeft = str2num(getParam('AstopLeft',listOfCoreParameters)); %Attenuation of left stop band (<FstopLeft) ripples in db
AstopRight = str2num(getParam('AstopRight',listOfCoreParameters)); %Attenuation of right stop band (>FstopRight) ripples in db

Apass_bp = Apass;
AstopLeft_bp = AstopLeft;
AstopRight_bp = AstopRight;

% Apass_lp = Apass;
% AstopRight_lp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;



StopToPassTransitionWidth_bp = str2num(getParam('StopToPassTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
PassToStopTransitionWidth_bp = str2num(getParam('PassToStopTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

% PassToStopTransitionWidth_lp = str2num(getParam('PassToStopTransitionWidth_lp',listOfCoreParameters)); %for low pass filter frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

StopToPassTransitionWidth_hp = str2num(getParam('StopToPassTransitionWidth_hp',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2


StopToPassTransitionWidth_hp_predownsample = str2num(getParam('StopToPassTransitionWidth_hp_predownsample',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2



UseFixedFilterOrder_bp = (getParam('UseFixedFilterOrder_bp',listOfCoreParameters));

% UseFixedFilterOrder_lp = (getParam('UseFixedFilterOrder_lp',listOfCoreParameters));

UseFixedFilterOrder_hp = (getParam('UseFixedFilterOrder_hp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_bp = (getParam('UseTwoPassAttenuationCorrection_bp',listOfCoreParameters));

% UseTwoPassAttenuationCorrection_lp = (getParam('UseTwoPassAttenuationCorrection_lp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_hp = (getParam('UseTwoPassAttenuationCorrection_hp',listOfCoreParameters));

FilterOrder_bp = str2num(getParam('FilterOrder_bp',listOfCoreParameters));

% FilterOrder_lp = str2num(getParam('FilterOrder_lp',listOfCoreParameters));

FilterOrder_hp = str2num(getParam('FilterOrder_hp',listOfCoreParameters));

MaximizeFilterOrderIfFixedFilterOrderIsUsed = str2num(getParam('MaximizeFilterOrderIfFixedFilterOrderIsUsed',listOfCoreParameters));

useTwoPassFiltering_bp = 'no';

% useTwoPassFiltering_lp = 'no';

useTwoPassFiltering_hp = 'no';

if ~isempty(strfind(core_cfg.bpfiltdir,'two'))
    useTwoPassFiltering_bp = 'yes';
end

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

if strcmp(UseFixedFilterOrder_bp,'yes') && strcmp(useTwoPassFiltering_bp,'yes') && (FilterOrder_bp > maxFilterOrder)
    error(['filter order for band pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_bp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_bp = maxFilterOrder;
end

% if strcmp(UseFixedFilterOrder_lp,'yes') && strcmp(useTwoPassFiltering_lp,'yes') && (FilterOrder_lp > maxFilterOrder)
% 	error(['filter order for low pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
% elseif (FilterOrder_lp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
% 	FilterOrder_lp = maxFilterOrder;
% end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end


if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
    error('band pass filter order must be an even number')
end

% if strcmp(UseFixedFilterOrder_lp,'yes') && logical(mod(FilterOrder_lp,2))
%     error('low pass order must be an even number')
% end

if strcmp(UseFixedFilterOrder_hp,'yes') && logical(mod(FilterOrder_hp,2))
    error('high pass order must be an even number')
end


% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
if strcmp(useTwoPassFiltering_bp,'yes') && strcmp(UseTwoPassAttenuationCorrection_bp,'yes')
    Apass_bp = Apass_bp/2;
    AstopLeft_bp = AstopLeft_bp/2;
    AstopRight_bp = AstopRight_bp/2;
end

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

tic
memtic
fprintf('SPD function initialized\n');
conseciDatas = 1:length(iDatas);
parfor conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    
    %iData = 1
    FrqOfSmplWishedPar = FrqOfSmplWished;
    
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    centerFreqFilter = listOfCenterFrequencies(iData);
    channelsOfInterest = listOfChannelsOfInterest(iData,:);
    channelsOfInterest = channelsOfInterest(~(cellfun(@isempty,channelsOfInterest)));
    signalMultiplicator = listOfDataSetSignalMultiplicator(iData);
    signalOffsetSamples = listOfDataSetOffsetSamples(iData);
    
    minFreq = centerFreqFilter - preCenterFreqFilterTo_FpassLeft;
    maxFreq = centerFreqFilter + postCenterFreqFilterTo_FpassRight;
    
    if minFreq < MinDetectionFrequency_FpassLeft
        minFreq = MinDetectionFrequency_FpassLeft;
    end
    
    if maxFreq > MaxDetectionFrequency_FpassRight
        maxFreq = MaxDetectionFrequency_FpassRight;
    end
    
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
        fprintf('dataset %i: resample data from %i to %i Hz\n',iData,data.fsample,FrqOfSmplWishedPar);
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
    
    
    
    
    cfg = [];
    cfg = core_cfg;
    cfg.bpfilter = 'yes';
    FpassLeft = minFreq; %left pass frequency in Hz
    FpassRight = maxFreq; %right pass frequency in Hz
    FstopLeft = FpassLeft - StopToPassTransitionWidth_bp; %left stop frequency in Hz
    FstopRight = FpassRight + PassToStopTransitionWidth_bp; %left stop frequency in Hz
    
    usedFilterOrder_bp = NaN;
    bp_hdm = NaN;
    if strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned')
        bp_d = [];
        bp_hd = [];
        fprintf('dataset %i: designing band pass filter for spd band\n',iData);
        if strcmp(UseFixedFilterOrder_bp,'yes')
            bp_d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2',FilterOrder_bp,FstopLeft,FpassLeft,FpassRight,FstopRight,FrqOfSmpl);
            bp_hd = design(bp_d,'equiripple');
        else
            bp_d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',FstopLeft,FpassLeft,FpassRight,FstopRight,AstopLeft_bp,Apass_bp,AstopRight_bp,FrqOfSmpl);
            bp_hd = design(bp_d,'equiripple','MinOrder', 'even');
        end
        usedFilterOrder_bp = bp_hd.order;
        cfg.bpfilterdesign = bp_hd;
        bp_hdm = measure(bp_hd);
    end
    if strcmp(UseFixedFilterOrder_bp,'yes')
        cfg.bpfiltord     = FilterOrder_bp;
    end
    cfg.bpfreq        = [FpassLeft FpassRight];%dummy values are overwritten by low level function
    cfg.feedback = core_cfg.feedback;
    fprintf('dataset %i: reprocess and apply band filter\n',iData);
    data = ft_fw_preprocessing(cfg,data);
 
    
    
    
    
    
    
    
    
    
    %ROI2
    epochLengthSamples = epochLength * FrqOfSmpl;
    [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
    
    if (signalOffsetSamples ~= 0)
        signalOffsetSamples_downsampled = floor(signalOffsetSeconds*FrqOfSmpl);
        roiBegins = roiBegins + signalOffsetSamples_downsampled;
        roiEnds = roiEnds + signalOffsetSamples_downsampled;
    end
    
    if strcmp(IgnoreDataSetHeader,'no')
        roiBegins = roiBegins(1:indexLastIncludedROIinData);
        roiEnds = roiEnds(1:indexLastIncludedROIinData);
        nSampleLength = floor(nSampleLength*(FrqOfSmplWishedPar/preDownsampleFreq));
        if (roiEnds(end) > nSampleLength)
            roiEnds(indexLastIncludedROIinData) = nSampleLength;
        end;
    end
    
    
    trlSampleBeginsAndEnds = [roiBegins roiEnds];
    
    trlSampleLengths = roiEnds - roiBegins + 1;
    
    sampleLengthsAcrossROIs = sum(trlSampleLengths);
    lengthsAcrossROIsSeconds = sampleLengthsAcrossROIs/FrqOfSmpl; % in seconds
    
    
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    
    
    smplsMinDetectionLength  = round(MinDetectionLength*FrqOfSmpl);
    smplsMaxDetectionLength  = round(MaxDetectionLength*FrqOfSmpl);
    
    minFreqPostFreqBorderBufferLength = round(FrqOfSmpl/minFreq);
    
    %smplsRMSTimeWndw = round(RMSTimeWndw*FrqOfSmpl);
    
    smplsRMSTimeWndw = 0;
    %assure odd samplesize for moving average time window
    if (mod(floor(RMSTimeWndw*FrqOfSmpl),2) == 0)
        smplsRMSTimeWndw = floor(RMSTimeWndw*FrqOfSmpl)+1;
    else
        smplsRMSTimeWndw = floor(RMSTimeWndw*FrqOfSmpl);
    end
    
    
    smplsMovAvgTimeWndw = 0;
    %assure odd samplesize for moving average time window
    if (mod(floor(MovAvgTimeWndw*FrqOfSmpl),2) == 0)
        smplsMovAvgTimeWndw = floor(MovAvgTimeWndw*FrqOfSmpl)+1;
    else
        smplsMovAvgTimeWndw = floor(MovAvgTimeWndw*FrqOfSmpl);
    end
    
    
    SDfrqBndPssSignal = [];
    if ~strcmp(UseAbsoluteEnvelopeThreshold,'yes')
        
        tempAllDataValues = [];
        
        for iTr = 1:size(data.trial,2)
            if iTr == 1
                tempAllDataValues = data.trial{iTr};
            else
                tempAllDataValues = cat(2,tempAllDataValues,data.trial{iTr});
            end
        end;
        
        if strcmp(ThresholdSignal,'filtered_signal') 
        elseif strcmp(ThresholdSignal,'envelope')
            for iChan_temp = 1:size(tempAllDataValues,1)
                tempAllDataValues(iChan_temp,:) = abs(hilbert(tempAllDataValues(iChan_temp,:)));
            end
        end

        
        if strcmp('meanoverchan',ThresholdAggregationMethod)
            if strcmp(ThresholdFormationBasis,'mean')
                SDfrqBndPssSignal = mean(mean(tempAllDataValues,2));
            elseif strcmp(ThresholdFormationBasis,'std')
                SDfrqBndPssSignal = mean(std(tempAllDataValues,0,2));
            end
        elseif strcmp('valuesoverchan',ThresholdAggregationMethod)
            if strcmp(ThresholdFormationBasis,'mean')
                SDfrqBndPssSignal = mean(tempAllDataValues(:));
            elseif strcmp(ThresholdFormationBasis,'std')
                SDfrqBndPssSignal = std(tempAllDataValues(:));
            end
        elseif strcmp('respectivechan',ThresholdAggregationMethod)
            if strcmp(ThresholdFormationBasis,'mean')
                SDfrqBndPssSignal = mean(tempAllDataValues,2);
            elseif strcmp(ThresholdFormationBasis,'std')
                SDfrqBndPssSignal = std(tempAllDataValues,0,2);
            end
        else
            error('ThresholdAggregationMethod for calculating mean or standard deviation over channels is unknown');
        end
        tempAllDataValues = [];%clear
    end
    
    if strcmp(OverwiteGlobalThresholdsAndUseIndividualDatasetEnvelopeTh,'yes')
        tempThBeginEnd = str2num(listOfIndividualDatasetEnvelopeThrehold{iData,1});
        tempThCriterion = str2num(listOfIndividualDatasetEnvelopeThrehold{iData,2});
        
        if ~strcmp(UseAbsoluteEnvelopeThreshold,'yes')
            dataset_factorThresholdBeginEnd = tempThBeginEnd;
            dataset_factorThresholdCriterion = tempThCriterion;
        end
    else
        dataset_factorThresholdBeginEnd = factorThresholdBeginEnd;
        dataset_factorThresholdCriterion = factorThresholdCriterion;
    end

    
    nChannels = length(data.label);
    
    maxPeaksOrTroughsPerSpindle = ceil(maxFreq*MaxDetectionLength+1);
    
    ch_detectedLengthSamples = [];
    ch_detectedBeginSample = [];
    ch_detectedEndSample = [];
    ch_detectedPeak2Peaks = [];
    ch_detectedPeaksSamples  = [];
    ch_detectedTroughsSamples  = [];
    ch_detectedSignalTroughsSamples = [];
    ch_detectedSignalPeaksSamples = [];
    ch_detectedSDofFilteredSignal = [];
    ch_detectedMergeCount = [];
    ch_nDetected = [];
    ch_nMerged = [];
    ch_densityPerEpoch = [];
    ch_SDfrqBndPssSignal = [];
    ch_detectedEnvelopeMaxs = [];
    ch_detectedEnvelopeMaxSamples = [];
    
    ch_detected_linear_regression_freq_slope = [];
    ch_detected_linear_regression_freq_offset = [];
    ch_detected_linear_regression_freq_R_squared = [];
    
    ch_detected_inst_freq_troughs = [];
    ch_detected_inst_freq_peaks = [];
    ch_detectedTroughsPotential = [];
    ch_detectedPeaksPotential = [];
    
    ch_contigSegment = [];
    
    for iChan = 1:nChannels
        %iChan = 1;
        
        fprintf('dataset %i: process channel %s ...\n',iData,data.label{iChan});
        if strcmp(UseAbsoluteEnvelopeThreshold,'yes')
                if strcmp(OverwiteGlobalThresholdsAndUseIndividualDatasetEnvelopeTh,'yes')
                    tempThBeginEnd = str2num(listOfIndividualDatasetEnvelopeThrehold{iData,1});
                    tempThCriterion = str2num(listOfIndividualDatasetEnvelopeThrehold{iData,2});
                    ch_SDfrqBndPssSignal{iChan} = tempThBeginEnd;
                    dataset_factorThresholdCriterion = tempThCriterion/tempThBeginEnd;
                else
                    ch_SDfrqBndPssSignal{iChan} = AbsoluteEnvelopeThresholdBeginEnd;
                end
        else
            if strcmp('respectivechan',ThresholdAggregationMethod)
                ch_SDfrqBndPssSignal{iChan} = SDfrqBndPssSignal(iChan);
            else
                ch_SDfrqBndPssSignal{iChan} = SDfrqBndPssSignal;
            end
        end
        
        %iChan = 1;
        cfg = [];
        cfg.feedback = core_cfg.feedback;
        cfg.channel = ft_channelselection(data.label{iChan}, data.label);
        chData = ft_selectdata(cfg,data);
        
        
        trl_detectedLengthSamples = [];
        trl_detectedBeginSample = [];
        trl_detectedEndSample = [];
        trl_detectedPeak2Peaks = [];
        trl_detectedPeaksSamples  = [];
        trl_detectedTroughsSamples  = [];
        trl_detectedSignalTroughsSamples = [];
        trl_detectedSignalPeaksSamples = [];
        trl_detectedSDofFilteredSignal = [];
        trl_detectedMergeCount = [];
        trl_nDetected = 0;
        trl_nMerged = 0;
        trl_detectedEnvelopeMaxs = [];
        trl_detectedEnvelopeMaxSamples = [];
        
        trl_detected_linear_regression_freq_slope = [];
        trl_detected_linear_regression_freq_offset = [];
        trl_detected_linear_regression_freq_R_squared = [];
        
        trl_detected_inst_freq_troughs = [];
        trl_detected_inst_freq_peaks = [];
        trl_detectedTroughsPotential = [];
        trl_detectedPeaksPotential = [];
        
        trl_contigSegment = [];
        
        for iTr = 1:size(chData.trial,2)
            %iTr = 4;
            fprintf('dataset %i: channel %s, subpart %i, preselect events in envelope\n',iData,data.label{iChan},iTr);
            rawDataSampleOffset = trlSampleBeginsAndEnds(iTr,1) - 1;
            
            frqBndPssSignal = chData.trial{iTr};
            frqBndPssSignal_hilbert = hilbert(frqBndPssSignal);
            
                   
            
            
            thresholdForDetectionBeginEnd = ch_SDfrqBndPssSignal{iChan}*dataset_factorThresholdBeginEnd;
            thresholdForDetectionCriterion = ch_SDfrqBndPssSignal{iChan}*dataset_factorThresholdCriterion;
            
                        
            
            lengthSignal = length(frqBndPssSignal);
            
            envelope = [];
            if strcmp(EnvelopeMethod,'hilbertEnv')
                envelope = abs(frqBndPssSignal_hilbert)';
            elseif strcmp(EnvelopeMethod,'smoothedRMSwd')
                envelope = smoothRMSwd(frqBndPssSignal,smplsRMSTimeWndw);
                if exist('smooth','file') == 2
                    envelope = smooth(envelope,smplsMovAvgTimeWndw);
                else
                    envelope = smoothwd(envelope,smplsMovAvgTimeWndw)';
                end
            end
            
%             plot(frqBndPssSignal,'b')
%             hold on
%             plot(smooth(envelope,smplsMovAvgTimeWndw),'g--')
%             plot(smooth(envelope2,smplsMovAvgTimeWndw),'r--')
%             plot(envelope2,'r')
%             hline(thresholdForDetection)

            
            %             RMSWindows = zeros(1,lengthSignal);
            %
            %             lastRMSSample = lengthSignal - smplsRMSTimeWndw + 1;
            %             tempAddISmpl = smplsRMSTimeWndw - 1;
            %             tempRMSWindows = zeros(lastRMSSample,smplsRMSTimeWndw);
            %             for iSmpl = 1:lastRMSSample
            %                 tempRMSWindows(iSmpl,:) = frqBndPssSignal(iSmpl:(iSmpl + tempAddISmpl));
            %             end
            %             firstSmplsRMStimeWndw = (smplsRMSTimeWndw+1)/2 ;
            %             lastSmplsRMStimeWndw = lengthSignal - ((smplsRMSTimeWndw-1)/2);
            %             RMSWindows(firstSmplsRMStimeWndw:lastSmplsRMStimeWndw) = rms(tempRMSWindows,2);
            %
            %             RMSWindows(1:(firstSmplsRMStimeWndw-1)) = smoothRMS(frqBndPssSignal(1:(firstSmplsRMStimeWndw-1)),smplsRMSTimeWndw);
            %             RMSWindows((lastSmplsRMStimeWndw+1):end) = smoothRMS(frqBndPssSignal((lastSmplsRMStimeWndw+1):end),smplsRMSTimeWndw);
            
            
            
            
            
            [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(envelope,thresholdForDetectionBeginEnd);
            
            
            %indicesValidSamples = find((begins >= firstSmplsRMStimeWndw) & (ends <= lastSmplsRMStimeWndw)); %consider border effects of RMS
            firstSmplsMinFreqPostFreqBorderBufferLength = minFreqPostFreqBorderBufferLength;
            lastSmplsMinFreqPostFreqBorderBufferLength = lengthSignal - minFreqPostFreqBorderBufferLength;
            
            indicesValidSamples = [];
            if strcmp(EnvelopeMethod,'hilbertEnv')
                indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects of filter
            elseif strcmp(EnvelopeMethod,'smoothedRMSwd')
                firstSmplsRMStimeWndw = (smplsRMSTimeWndw+1)/2 ;
                lastSmplsRMStimeWndw = lengthSignal - ((smplsRMSTimeWndw-1)/2);
                indicesValidSamples = find((begins >= firstSmplsRMStimeWndw) & (ends <= lastSmplsRMStimeWndw) & (begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects of RMS and filter
            end
            
            begins = begins(indicesValidSamples);
            ends = ends(indicesValidSamples);
            
            tempCandidatesLengths = ends - begins + 1;
            
            indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
            
            smplsMergeEventsInProximityWithinDetectionMargins = round(MergeEventsInProximityWithinDetectionMargins*FrqOfSmpl);
            
            nDetected = length(indicesCandiates);
            tempNmergedOuter = 0;
            detectedMergeCount = zeros(1,nDetected);
            
            if (nDetected > 0)
                %detectedLengthSamples = tempCandidatesLengths(indicesCandiates);
                detectedBeginSample = begins(indicesCandiates);
                detectedEndSample = ends(indicesCandiates);
                
                if smplsMergeEventsInProximityWithinDetectionMargins > 0
                    fprintf('dataset %i: channel %s, subpart %i, merge events\n',iData,data.label{iChan},iTr);
                    %%%%%%%%
                    tempNmergedInner = 1;
                    
                    while (nDetected > 1) && (tempNmergedInner > 0)
                        tempNmergedInner = 0;
                        
                        candidateDistanceSamples = detectedBeginSample(2:end) - detectedEndSample(1:(end-1));
                        
                        [sortCandidateDistanceSamples,preMergeDetectionIndices] = sort(candidateDistanceSamples);
                        
                        %preMergeDetectionIndices = uint64(preMergeDetectionIndices);
                        %preMergeDetectionIndices = 1:nDetected;
                        
                        tempMergeDetectionIndices = [];
                        
                        iIterCand = 1;
                        
                        
                        while (iIterCand <= length(preMergeDetectionIndices))
                            %iIterCand = 1
                            %if iIterCand < length(preMergeDetectionIndices)
                            if (  ((detectedEndSample(preMergeDetectionIndices(iIterCand)) + smplsMergeEventsInProximityWithinDetectionMargins) >= detectedBeginSample(preMergeDetectionIndices(iIterCand)+1)) && ...
                                    (((detectedEndSample(preMergeDetectionIndices(iIterCand)+1)-detectedBeginSample(preMergeDetectionIndices(iIterCand)) + 1)) <= smplsMaxDetectionLength) )
                                tempNmergedInner = tempNmergedInner + 1;
                                tempNmergedOuter = tempNmergedOuter + 1;
                                tempMergeDetectionIndices(tempNmergedInner) = preMergeDetectionIndices(iIterCand);
                                preMergeDetectionIndices(preMergeDetectionIndices == (preMergeDetectionIndices(iIterCand)+1)) = [];
                                preMergeDetectionIndices(preMergeDetectionIndices == (preMergeDetectionIndices(iIterCand)-1)) = [];
                            end
                            iIterCand = iIterCand + 1;
                        end
                        
                        detectedEndSample(tempMergeDetectionIndices) = detectedEndSample(tempMergeDetectionIndices+1);
                        detectedMergeCount(tempMergeDetectionIndices) = detectedMergeCount(tempMergeDetectionIndices) + detectedMergeCount(tempMergeDetectionIndices+1) + 1;
                        detectedMergeCount(tempMergeDetectionIndices+1) = [];
                        detectedBeginSample(tempMergeDetectionIndices+1) = [];
                        detectedEndSample(tempMergeDetectionIndices+1) = [];
                        
                        nDetected = nDetected - tempNmergedInner;
                        
                    end
                end
                detectedLengthSamples = detectedEndSample - detectedBeginSample + 1 ;
                
                
                %%%%%%%%
                
                
                % %                 tempPostMergeDetectedLengthSamples = zeros(1,nDetected);
                % %                 tempPostMergeDetectedBeginSample = zeros(1,nDetected);
                % %                 tempPostMergeDetectedEndSample = zeros(1,nDetected);
                
                minPeakDistanceSamples = ceil(((1/(maxFreq)) * FrqOfSmpl)/2); % half of max freq of interest in samples
                
                detectedPeak2Peaks = zeros(1,nDetected);
                detectedPeaksSamples  = zeros(1,nDetected);
                detectedTroughsSamples  = zeros(1,nDetected);
                detectedSignalTroughsSamples = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detectedSignalPeaksSamples = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detectedSignalTroughsSamples(:,:) = NaN;
                detectedSignalPeaksSamples(:,:) = NaN;
                detectedSDofFilteredSignal = zeros(1,nDetected);
                
                detectedEnvelopeMaxs = zeros(1,nDetected);
                detectedEnvelopeMaxSamples = zeros(1,nDetected);
                
                detected_linear_regression_freq_slope = zeros(1,nDetected);
                detected_linear_regression_freq_offset = zeros(1,nDetected);
                detected_linear_regression_freq_R_squared = zeros(1,nDetected);
                
                detected_inst_freq_troughs = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detected_inst_freq_peaks = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detected_inst_freq_troughs(:,:) = NaN;
                detected_inst_freq_peaks(:,:) = NaN;
                
                detectedTroughsPotential = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detectedPeaksPotential = zeros(nDetected,maxPeaksOrTroughsPerSpindle);
                detectedTroughsPotential(:,:) = NaN;
                detectedPeaksPotential(:,:) = NaN;

                % %                 detectedMergeCount = zeros(1,nDetected);
                
                fprintf('dataset %i: channel %s, subpart %i, annotate events\n',iData,data.label{iChan},iTr);
                for iIterCand = 1:nDetected
                    % %                 preMergeDetectionIndices = 1:nDetected;
                    % %                 postMergeDetectionIndicesBegins = [];
                    % %                 postMergeDetectionIndicesEnds = [];
                    % %                 iIterCand = 1;
                    % %                 tempNmerged = 0;
                    % %                 tempMergeCount = 0;
                    % %                 %tempMergeDistances = [];%xxxx
                    % %                 while (iIterCand <= length(preMergeDetectionIndices))
                    %iIterCand = 1
                    
                    
                    % %                     if iIterCand < length(preMergeDetectionIndices)
                    % %                        if (  ((detectedEndSample(preMergeDetectionIndices(iIterCand)) + smplsMergeEventsInProximityWithinDetectionMargins) >= detectedBeginSample(preMergeDetectionIndices(iIterCand+1))) && ...
                    % %                                (((detectedEndSample(preMergeDetectionIndices(iIterCand+1))-detectedBeginSample(preMergeDetectionIndices(iIterCand)) + 1)) <= smplsMaxDetectionLength) )
                    % %                            postMergeDetectionIndicesBegins(iIterCand) = preMergeDetectionIndices(iIterCand);
                    % %                            postMergeDetectionIndicesEnds(iIterCand) = preMergeDetectionIndices(iIterCand+1);
                    % %                            tempNmerged = tempNmerged + 1;
                    % %                            tempMergeCount = tempMergeCount + 1;
                    % %                            %tempMergeDistances(tempMergeCount) = detectedBeginSample(preMergeDetectionIndices(iIterCand+1)) - detectedEndSample(preMergeDetectionIndices(iIterCand));%xxxx
                    % %                            preMergeDetectionIndices(iIterCand+1) = [];
                    % %                            continue;
                    % %                        end
                    % %                     end
                    % %
                    % %                     if tempMergeCount < 1
                    % %                         postMergeDetectionIndicesBegins(iIterCand) = preMergeDetectionIndices(iIterCand);
                    % %                         postMergeDetectionIndicesEnds(iIterCand) = preMergeDetectionIndices(iIterCand);
                    % %                     end
                    %postMergeDetectionIndices(iIterCand) = preMergeDetectionIndices(iIterCand);
                    
                    currentRawDataSampleOffset = rawDataSampleOffset + detectedBeginSample(iIterCand) - 1;
                    candSignal = frqBndPssSignal(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                    candSignal_hilbert = frqBndPssSignal_hilbert(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                    
                    
                    tempCandSignalminSample = find(min(candSignal) == candSignal);
                    tempCandSignalmaxSample = find(max(candSignal) == candSignal);
                    
                    
                    candSignalminSample = currentRawDataSampleOffset + tempCandSignalminSample;
                    candSignalmaxSample = currentRawDataSampleOffset + tempCandSignalmaxSample;
                    
                    candSignalmin = candSignal(tempCandSignalminSample);
                    candSignalmax = candSignal(tempCandSignalmaxSample);
                    candPeak2Peak = candSignalmax - candSignalmin;
                    
                    %candSignal = candSignal + 5e-6;
                    
                    [tempCandSignalPeaks, tempCandSignalPeaksSamples] = findpeaks(candSignal,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',minPeakDistanceSamples);
                    candSignalPeaks = candSignal(tempCandSignalPeaksSamples);
                    candSignalPeaksSamples = currentRawDataSampleOffset + tempCandSignalPeaksSamples;
                    
                    tempAlignFactor = 1.0001*max(candSignal);
                    candSignalInv = tempAlignFactor - candSignal;
                    [tempCandSignalTroughs, tempCandSignalTroughsSamples] = findpeaks(candSignalInv,'MINPEAKHEIGHT',tempAlignFactor,'MINPEAKDISTANCE',minPeakDistanceSamples);
                    candSignalTroughs = candSignal(tempCandSignalTroughsSamples);
                    %validBelowZero = find(candSignalTroughs < 0);
                    candSignalTroughsSamples = currentRawDataSampleOffset + tempCandSignalTroughsSamples;
                    
                    nCandSignalTroughs = length(candSignalTroughsSamples);
                    nCandSignalPeaks = length(candSignalPeaksSamples);
                    
                    
                    candSignalEnvelope = envelope(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                    tempCandSignalEnvelopemaxSample = find(max(candSignalEnvelope) == candSignalEnvelope);
                    candSignalEnvelopemax = candSignalEnvelope(tempCandSignalEnvelopemaxSample);
                    candSignalEnvelopemaxSample = currentRawDataSampleOffset + tempCandSignalEnvelopemaxSample;
                    
                    %chirp
%                     tempConsecTroughsPeaks = sort([tempCandSignalPeaksSamples tempCandSignalTroughsSamples]/FrqOfSmpl);
%                     temp_inst_freq2 = 1./(diff(tempConsecTroughsPeaks)*2);
%                     
%                     tempConsecTroughs = sort([tempCandSignalTroughsSamples]/FrqOfSmpl);
%                     temp_inst_freq2_troughs = 1./(diff(tempConsecTroughs));
%                     
%                     tempConsecPeaks = sort([tempCandSignalPeaksSamples]/FrqOfSmpl);
%                     temp_inst_freq2_peaks = 1./(diff(tempConsecPeaks));
                    
                    
                    %inst_freq = (diff(unwrap(angle(hilbert(candSignal))))/(2*pi))*FrqOfSmpl;
                    temp_inst_freq = (diff(unwrap(angle(candSignal_hilbert)))/(2*pi))*FrqOfSmpl;

                    temp_time = (1:(length(candSignal)))./FrqOfSmpl;
                    
                    %[r_time_freq,p_time_freq] = corrcoef([temp_time(2:end)' temp_inst_freq']);
                    
                    
                    temp_time_freq_regression = fitlm(temp_time(2:end),temp_inst_freq);
%                     temp_time_freq_regression.Coefficients.Estimate(1)
%                     temp_time_freq_regression.Coefficients.Estimate(2)
%                     temp_time_freq_regression.Rsquared.Ordinary
                    
                    
                    temp_linear_regression_freq_slope = temp_time_freq_regression.Coefficients.Estimate(2);
                    temp_linear_regression_freq_offset = temp_time_freq_regression.Coefficients.Estimate(1);
                    temp_linear_regression_freq_R_squared = temp_time_freq_regression.Rsquared.Ordinary;
                    
                    temp_inst_freq_troughs = temp_inst_freq(tempCandSignalTroughsSamples);
                    temp_inst_freq_peaks = temp_inst_freq(tempCandSignalPeaksSamples);
                    
                    
%                     figure
%                     subplot(2,1,1);
%                     plot(temp_time(2:end),temp_inst_freq); hold on;
%                     plot(temp_time(2:end),temp_time_freq_regression.Coefficients.Estimate(2) * temp_time(2:end) + temp_time_freq_regression.Coefficients.Estimate(1))
%                     hold off;
%                     
%                     subplot(2,1,2);
%                     plot(temp_time(2:end),candSignal(2:end))
               


 
                    
                    %                     figure
                    %                     plot((currentRawDataSampleOffset+1):(currentRawDataSampleOffset+(detectedEndSample(iIterCand)-detectedBeginSample(iIterCand))+1),candSignal)
                    %                     hold on
                    %                     %plot(candSignalPeaksSamples,candSignalPeaks,'g*')
                    %                     %plot(begins(iCand):minPeakDistanceSamples:ends(iCand),0,'r*')
                    %                     plot(candSignalminSample,candSignalmin,'b*',candSignalmaxSample,candSignalmax,'b*')
                    %                     plot(candSignalTroughsSamples,candSignalTroughs,'r*',candSignalPeaksSamples,candSignalPeaks,'g*')
                    
                    detectedPeak2Peaks(iIterCand) = candPeak2Peak;
                    detectedPeaksSamples(iIterCand) = candSignalmaxSample;
                    detectedTroughsSamples(iIterCand) = candSignalminSample;
                    detectedSignalTroughsSamples(iIterCand,1:nCandSignalTroughs) = candSignalTroughsSamples;
                    detectedSignalPeaksSamples(iIterCand,1:nCandSignalPeaks) = candSignalPeaksSamples;
                    detectedSDofFilteredSignal(iIterCand) = std(candSignal);
                    
                    detectedEnvelopeMaxs(iIterCand) = candSignalEnvelopemax;
                    detectedEnvelopeMaxSamples(iIterCand) = candSignalEnvelopemaxSample;
                    
                    
                    
                    detected_linear_regression_freq_slope(iIterCand) = temp_linear_regression_freq_slope;
                    detected_linear_regression_freq_offset(iIterCand) = temp_linear_regression_freq_offset;
                    detected_linear_regression_freq_R_squared(iIterCand) = temp_linear_regression_freq_R_squared;
                    
                    detected_inst_freq_troughs(iIterCand,1:nCandSignalTroughs) = temp_inst_freq_troughs;
                    detected_inst_freq_peaks(iIterCand,1:nCandSignalPeaks) = temp_inst_freq_peaks;
                    
                    detectedTroughsPotential(iIterCand,1:nCandSignalTroughs) = candSignalTroughs;
                    detectedPeaksPotential(iIterCand,1:nCandSignalPeaks) = candSignalPeaks;                    

                    
                    
                    % %                     detectedMergeCount(iIterCand) = tempMergeCount;
                    
                    % %                     iIterCand = iIterCand + 1;
                    % %                     tempMergeCount = 0;
                    %tempMergeDistances = [];%xxxx
                end
                
                
                
                
                trl_detectedBeginSample = cat(2,trl_detectedBeginSample,detectedBeginSample + rawDataSampleOffset - 1);
                trl_detectedEndSample = cat(2,trl_detectedEndSample,detectedEndSample + rawDataSampleOffset - 1);
                trl_detectedLengthSamples = cat(2,trl_detectedLengthSamples,detectedLengthSamples);
                
                
                
                trl_detectedPeak2Peaks = cat(2,trl_detectedPeak2Peaks,detectedPeak2Peaks);
                trl_detectedPeaksSamples = cat(2,trl_detectedPeaksSamples,detectedPeaksSamples);
                trl_detectedTroughsSamples  = cat(2,trl_detectedTroughsSamples ,detectedTroughsSamples);
                trl_detectedSignalTroughsSamples = cat(1,trl_detectedSignalTroughsSamples,detectedSignalTroughsSamples);
                trl_detectedSignalPeaksSamples = cat(1,trl_detectedSignalPeaksSamples,detectedSignalPeaksSamples);
                trl_detectedSDofFilteredSignal = cat(2,trl_detectedSDofFilteredSignal,detectedSDofFilteredSignal);
                trl_detectedMergeCount = cat(2,trl_detectedMergeCount,detectedMergeCount);
                trl_nDetected = trl_nDetected + nDetected ;
                trl_nMerged = trl_nMerged + tempNmergedOuter; % unused
                trl_detectedEnvelopeMaxs = cat(2,trl_detectedEnvelopeMaxs,detectedEnvelopeMaxs);
                trl_detectedEnvelopeMaxSamples = cat(2,trl_detectedEnvelopeMaxSamples,detectedEnvelopeMaxSamples);
                
                trl_detected_linear_regression_freq_slope = cat(2,trl_detected_linear_regression_freq_slope,detected_linear_regression_freq_slope);
                trl_detected_linear_regression_freq_offset = cat(2,trl_detected_linear_regression_freq_offset,detected_linear_regression_freq_offset);
                trl_detected_linear_regression_freq_R_squared = cat(2,trl_detected_linear_regression_freq_R_squared,detected_linear_regression_freq_R_squared);
                
                trl_detected_inst_freq_troughs = cat(1,trl_detected_inst_freq_troughs,detected_inst_freq_troughs);
                trl_detected_inst_freq_peaks = cat(1,trl_detected_inst_freq_peaks,detected_inst_freq_peaks);
                
                trl_detectedTroughsPotential = cat(1,trl_detectedTroughsPotential,detectedTroughsPotential);
                trl_detectedPeaksPotential = cat(1,trl_detectedPeaksPotential,detectedPeaksPotential);  
                
                trl_contigSegment = cat(2,trl_contigSegment,repmat(iTr,1,nDetected));

                
            end
        end
        
        fprintf('dataset %i: channel %s, select events\n',iData,data.label{iChan});
        tempIndexWithinThresholds = find((trl_detectedPeak2Peaks >= MinAbsoluteDownToUpPeakPotential) & (trl_detectedPeak2Peaks <= MaxAbsoluteDownToUpPeakPotential) & (trl_detectedEnvelopeMaxs >= thresholdForDetectionCriterion));
        
        
        %         %MergeEventsInProximityWithinDetectionMargins
        %
        %         smplsMergeEventsInProximityWithinDetectionMargins = round(MergeEventsInProximityWithinDetectionMargins*FrqOfSmpl);
        %         iAboveTh = 1;
        %         while (iAboveTh < length(tempIndexAboveMeanThreshold))
        %              if ( ((trl_detectedEndSample(tempIndexAboveMeanThreshold(iAboveTh)) + smplsMergeEventsInProximityWithinDetectionMargins) >= trl_detectedBeginSample(tempIndexAboveMeanThreshold(iAboveTh+1))) ...
        %                 && (((trl_detectedEndSample(tempIndexAboveMeanThreshold(iAboveTh+1))-trl_detectedBeginSample(tempIndexAboveMeanThreshold(iAboveTh)) + 1)) <= smplsMaxDetectionLength) )
        %
        %                 trl_detectedEndSample(tempIndexAboveMeanThreshold(iAboveTh)) = trl_detectedEndSample(tempIndexAboveMeanThreshold(iAboveTh+1));
        %
        %                 [tempMax, tempIndexMax]= max(trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold(iAboveTh)),trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold(iAboveTh+1)));
        %                 trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold(iAboveTh)) = tempMax;
        %
        %                 trl_detectedPeaksSamples(tempIndexAboveMeanThreshold(iAboveTh)) =  trl_detectedPeaksSamples(tempIndexAboveMeanThreshold(iAboveTh+(tempIndexMax-1)));
        %                 trl_detectedTroughsSamples(tempIndexAboveMeanThreshold(iAboveTh)) =  trl_detectedTroughsSamples(tempIndexAboveMeanThreshold(iAboveTh+(tempIndexMax-1)));
        %                 trl_detectedSignalTroughsSamples(tempIndexAboveMeanThreshold(iAboveTh),:)  = [trl_detectedSignalTroughsSamples(tempIndexAboveMeanThreshold(iAboveTh),:) , trl_detectedSignalTroughsSamples(tempIndexAboveMeanThreshold(iAboveTh+1),:)]
        %                 trl_detectedSignalPeaksSamples(tempIndexAboveMeanThreshold(iAboveTh),:)  = [trl_detectedSignalPeaksSamples(tempIndexAboveMeanThreshold(iAboveTh),:) , trl_detectedSignalPeaksSamples(tempIndexAboveMeanThreshold(iAboveTh+1),:)]
        %
        %                 tempW = [trl_detectedLengthSamples(tempIndexAboveMeanThreshold(iAboveTh)) trl_detectedLengthSamples(tempIndexAboveMeanThreshold(iAboveTh+1))];
        %                 trl_detectedSDofFilteredSignal(tempIndexAboveMeanThreshold(iAboveTh))  = mean(tempW.*[trl_detectedSDofFilteredSignal(tempIndexAboveMeanThreshold(iAboveTh)),trl_detectedSDofFilteredSignal(tempIndexAboveMeanThreshold(iAboveTh+1))])/mean(tempW);
        %
        %                 trl_detectedLengthSamples(tempIndexAboveMeanThreshold(iAboveTh)) = trl_detectedEndSample(tempIndexAboveMeanThreshold(iAboveTh)) - trl_detectedBeginSample(tempIndexAboveMeanThreshold(iAboveTh)) + 1;
        %
        %
        %                 tempIndexAboveMeanThreshold(iAboveTh+1) = [];
        %             else
        %                 iAboveTh = iAboveTh + 1;
        %             end
        %
        %         end
        
        
        
        
        ch_detectedLengthSamples{iChan} = trl_detectedLengthSamples(tempIndexWithinThresholds);
        ch_detectedBeginSample{iChan} = trl_detectedBeginSample(tempIndexWithinThresholds);
        ch_detectedEndSample{iChan} = trl_detectedEndSample(tempIndexWithinThresholds);
        ch_detectedPeak2Peaks{iChan} = trl_detectedPeak2Peaks(tempIndexWithinThresholds);
        ch_detectedPeaksSamples{iChan} = trl_detectedPeaksSamples(tempIndexWithinThresholds);
        ch_detectedTroughsSamples{iChan}  = trl_detectedTroughsSamples(tempIndexWithinThresholds);
        ch_detectedSignalTroughsSamples{iChan} = trl_detectedSignalTroughsSamples(tempIndexWithinThresholds,:);
        ch_detectedSignalPeaksSamples{iChan} = trl_detectedSignalPeaksSamples(tempIndexWithinThresholds,:);
        ch_detectedSDofFilteredSignal{iChan} = trl_detectedSDofFilteredSignal(tempIndexWithinThresholds);
        ch_detectedMergeCount{iChan} = trl_detectedMergeCount(tempIndexWithinThresholds);
        ch_nDetected{iChan} = length(tempIndexWithinThresholds);
        ch_nMerged{iChan} = sum(trl_detectedMergeCount(tempIndexWithinThresholds));
        ch_densityPerEpoch{iChan} = length(tempIndexWithinThresholds)/(lengthsAcrossROIsSeconds/epochLength);
        ch_detectedEnvelopeMaxs{iChan} = trl_detectedEnvelopeMaxs(tempIndexWithinThresholds);
        ch_detectedEnvelopeMaxSamples{iChan} = trl_detectedEnvelopeMaxSamples(tempIndexWithinThresholds);
        
        ch_detected_linear_regression_freq_slope{iChan} = trl_detected_linear_regression_freq_slope(tempIndexWithinThresholds);
        ch_detected_linear_regression_freq_offset{iChan} = trl_detected_linear_regression_freq_offset(tempIndexWithinThresholds);
        ch_detected_linear_regression_freq_R_squared{iChan} = trl_detected_linear_regression_freq_R_squared(tempIndexWithinThresholds);
        
        ch_detected_inst_freq_troughs{iChan} = trl_detected_inst_freq_troughs(tempIndexWithinThresholds,:);
        ch_detected_inst_freq_peaks{iChan} = trl_detected_inst_freq_peaks(tempIndexWithinThresholds,:);
        
        ch_detectedTroughsPotential{iChan} = trl_detectedTroughsPotential(tempIndexWithinThresholds,:);
        ch_detectedPeaksPotential{iChan} = trl_detectedPeaksPotential(tempIndexWithinThresholds,:);
        
        ch_contigSegment{iChan} = trl_contigSegment(tempIndexWithinThresholds);
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
        
        if strcmp(core_cfg.hpfilttype,'but')
            if strcmp(UseFixedFilterOrder_hp,'yes')
                usedFilterOrder_hp = FilterOrder_hp;
            else
                usedFilterOrder_hp = 6;
            end
        end
    end
    
    if ~(strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned'))
        
        usedFilterOrder_bp = NaN;
        bp_hdm.Fs = FrqOfSmpl;
        bp_hdm.Astop1 = NaN;
        bp_hdm.TransitionWidth1 = NaN;
        bp_hdm.F3dB1 = minFreq;
        bp_hdm.F6dB1 = NaN;
        bp_hdm.Fpass1 = NaN;
        bp_hdm.Apass = NaN;
        bp_hdm.Fpass2 = minFreq;
        bp_hdm.F3dB2 = NaN;
        bp_hdm.F6dB2 = NaN;
        bp_hdm.TransitionWidth2 = NaN;
        bp_hdm.Astop2 = NaN;
        if strcmp(core_cfg.bpfilttype,'but')
            if strcmp(UseFixedFilterOrder_bp,'yes')
                usedFilterOrder_bp = FilterOrder_bp;
            else
                usedFilterOrder_bp = 4;
            end
        end
    end
    
    
    bp_f_type_detail = '';
    switch core_cfg.bpfilttype
        case 'but'
            bp_f_type_detail = 'IIR_Butterworth_ml_butter';
        case 'fir'
            bp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
        case 'FIRdesigned'
            bp_f_type_detail = 'FIR_equiripple_signal_toolbox';
        case 'IIRdesigned'
            bp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
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
    
    
    
    fidf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
    %write header
    fprintf(fidf,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'datasetnum','dataset',...
        'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB',...
        'bp_filter','bp_filter_type','bp_dir_and_passing','usedFilterOrder_bp','bp_Fs_Hz','bp_Astop1_dB','bp_Fstop1_Hz','bp_TransitionWidth1_Hz','bp_F3dB1_Hz','bp_F6dB1_Hz','bp_Fpass1_Hz','bp_Apass_dB','bp_Fpass2_Hz','bp_F3dB2_Hz','bp_F6dB2_Hz','bp_TransitionWidth2_Hz','bp_Fstop2_Hz','bp_Astop2_dB');
    %write content
    fprintf(fidf,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e,%f,%f,%f,%f,%f,%e' '\n'],...
        iData,datasetsPath,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass,...
        core_cfg.bpfilttype,bp_f_type_detail,core_cfg.bpfiltdir,usedFilterOrder_bp,bp_hdm.Fs,bp_hdm.Astop1,bp_hdm.Fstop1,bp_hdm.TransitionWidth1,bp_hdm.F3dB1,bp_hdm.F6dB1,bp_hdm.Fpass1,bp_hdm.Apass,bp_hdm.Fpass2,bp_hdm.F3dB2,bp_hdm.F6dB2,bp_hdm.TransitionWidth2,bp_hdm.Fstop2,bp_hdm.Astop2);
    
    
    fclose(fidf);
    
    
    fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
    fide = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_events_' 'datanum_' num2str(iData) '.csv'],'wt');
    fidp = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_' 'datanum_' num2str(iData) '.csv'],'wt');
    fidt = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    fidpf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_freq_' 'datanum_' num2str(iData) '.csv'],'wt');
    fidtf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_freq_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    fidpa = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_filtered_potential_' 'datanum_' num2str(iData) '.csv'],'wt');
    fidta = fopen([pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_filtered_potential_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    %write header of outputfiles
    fprintf(fidc,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','datasetnum','channel','count','density_per_epoch','mean_duration_seconds','mean_amplitude_trough2peak_potential','mean_frequency_by_mean_pk_trgh_cnt_per_dur','epoch_length_seconds','merged_count','lengths_ROI_seconds','used_threshold_basis','used_factor_for_threshold_basis_begin_end','used_factor_for_threshold_basis_criterion','used_min_detection_pass_or_cutoff_freq','used_max_detection_pass_or_cutoff_freq','used_center_freq','mean_SD_of_filtered_signal','mean_troughs_per_event','mean_peaks_per_event','mean_linear_regression_freq_slope','mean_linear_regression_freq_offset','mean_linear_regression_freq_R_squared');
    
    fprintf(fide,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        'datasetnum','channel','duration_seconds','amplitude_peak2trough_max','frequency_by_mean_pk_trgh_cnt_per_dur','duration_samples','sample_begin','sample_end','sample_peak_max','sample_trough_max','envelope_max'...
        ,'dataset','hypnogram','used_stages_for_detection','used_threshold_basis','used_min_pass_or_cutoff_detection_freq','used_max_detection_pass_or_cutoff_freq','used_center_freq','seconds_begin','seconds_end','seconds_peak_max','seconds_trough_max','seconds_envelope_max','id_within_channel',...
        'stage','stage_alt','stage_alt2',...
        'contig_segment',...
        'SD_of_filtered_signal','merged_count',...
        'linear_regression_freq_slope','linear_regression_freq_offset','linear_regression_freq_R_squared',...
        'number_troughs','number_peaks');
    
    
    fprintf(fidp,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    fprintf(fidt,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    for iPT = 1:(maxPeaksOrTroughsPerSpindle-1)
        fprintf(fidp,'%s,',['extr_' num2str(iPT)]);
        fprintf(fidt,'%s,',['extr_' num2str(iPT)]);
    end
    fprintf(fidp,'%s\n',['extr_' num2str(maxPeaksOrTroughsPerSpindle)]);
    fprintf(fidt,'%s\n',['extr_' num2str(maxPeaksOrTroughsPerSpindle)]);
    
    
    fprintf(fidpf,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    fprintf(fidtf,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    for iPT = 1:(maxPeaksOrTroughsPerSpindle-1)
        fprintf(fidpf,'%s,',['extr_freq_' num2str(iPT)]);
        fprintf(fidtf,'%s,',['extr_freq_' num2str(iPT)]);
    end
    fprintf(fidpf,'%s\n',['extr_freq_' num2str(maxPeaksOrTroughsPerSpindle)]);
    fprintf(fidtf,'%s\n',['extr_freq_' num2str(maxPeaksOrTroughsPerSpindle)]);
    
    
    fprintf(fidpa,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    fprintf(fidta,'%s,%s,%s,','datasetnum','channel','id_within_channel');
    for iPT = 1:(maxPeaksOrTroughsPerSpindle-1)
        fprintf(fidpa,'%s,',['extr_pot_' num2str(iPT)]);
        fprintf(fidta,'%s,',['extr_pot_' num2str(iPT)]);
    end
    fprintf(fidpa,'%s\n',['extr_pot_' num2str(maxPeaksOrTroughsPerSpindle)]);
    fprintf(fidta,'%s\n',['extr_pot_' num2str(maxPeaksOrTroughsPerSpindle)]);
    
    
    
    
    
    for iChan = 1:nChannels
        
        ch = data.label{iChan};
        
        epochs = {};
        for iDet = 1:length(ch_detectedTroughsSamples{iChan})
            tempSample = ch_detectedTroughsSamples{iChan}(iDet);
            tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples));
            if ~any(tempInd)
            	tempInd = ((hypnEpochsBeginsSamples <= tempSample+1) & (tempSample <= hypnEpochsEndsSamples));
            end
             if ~any(tempInd)
            	tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample-1 <= hypnEpochsEndsSamples));
            end
            epochs(iDet,:) = [hypnStages(tempInd,1) ...
                              hypnStages(tempInd,2) ...
                              hypnStages(tempInd,3)];
        end;
        
        
        fprintf(fidc,'%i,',iData);
        fprintf(fidc,'%s,',ch);
        fprintf(fidc,'%i,',ch_nDetected{iChan});
        fprintf(fidc,'%f,',ch_densityPerEpoch{iChan});
        tempLengthMeanSeconds = mean(ch_detectedLengthSamples{iChan}'/FrqOfSmpl);
        fprintf(fidc,'%f,',tempLengthMeanSeconds);
        fprintf(fidc,'%e,',mean(ch_detectedPeak2Peaks{iChan}));
        tempNtroughsMean = (length(find(~isnan(ch_detectedSignalTroughsSamples{iChan}))))/size(ch_detectedSignalTroughsSamples{iChan},1);
        tempNpeaksMean = (length(find(~isnan(ch_detectedSignalPeaksSamples{iChan}))))/size(ch_detectedSignalPeaksSamples{iChan},1);
        fprintf(fidc,'%f,', ((tempNtroughsMean + tempNpeaksMean)/2)/tempLengthMeanSeconds);
        fprintf(fidc,'%f,',epochLength);
        fprintf(fidc,'%i,',ch_nMerged{iChan});
        fprintf(fidc,'%f,',lengthsAcrossROIsSeconds);
        fprintf(fidc,'%e,',ch_SDfrqBndPssSignal{iChan});
        fprintf(fidc,'%f,',dataset_factorThresholdBeginEnd);
        fprintf(fidc,'%f,',dataset_factorThresholdCriterion);
        fprintf(fidc,'%f,',minFreq);
        fprintf(fidc,'%f,',maxFreq);
        fprintf(fidc,'%f,',centerFreqFilter);

        fprintf(fidc,'%e,',mean(ch_detectedSDofFilteredSignal{iChan}));
        
        fprintf(fidc,'%f,', tempNtroughsMean);
        fprintf(fidc,'%f,', tempNpeaksMean);
        
        fprintf(fidc,'%f,', mean(ch_detected_linear_regression_freq_slope{iChan}));
        fprintf(fidc,'%f,', mean(ch_detected_linear_regression_freq_offset{iChan}));
        fprintf(fidc,'%f\n', mean(ch_detected_linear_regression_freq_R_squared{iChan}));


        
        
        
        if ch_nDetected{iChan} > 0
            output = cell(ch_nDetected{iChan},11);
            
            output(:,1) = cellstr(repmat(ch, ch_nDetected{iChan}, 1));
            output(:,2) = num2cell(ch_detectedLengthSamples{iChan}');
            output(:,3) = num2cell(ch_detectedBeginSample{iChan}');
            output(:,4) = num2cell(ch_detectedEndSample{iChan}');
            output(:,5) = num2cell(ch_detectedPeaksSamples{iChan}');
            output(:,6) = num2cell(ch_detectedTroughsSamples{iChan}');
            
            output(:,7) = num2cell(ch_detectedPeak2Peaks{iChan});
            output(:,8) = num2cell(ch_detectedEnvelopeMaxs{iChan});
            
            output(:,9) = cellstr(repmat(datasetsPath, ch_nDetected{iChan}, 1));
            output(:,10) = cellstr(repmat(hypnogramPath, ch_nDetected{iChan}, 1));
            output(:,11) = cellstr(repmat(strjoin(sleepStagesOfInterest,' '), ch_nDetected{iChan}, 1));
            output(:,12) = num2cell(repmat(ch_SDfrqBndPssSignal{iChan}, ch_nDetected{iChan}, 1));
            
            output(:,13) = num2cell(repmat(minFreq, ch_nDetected{iChan}, 1));
            output(:,14) = num2cell(repmat(maxFreq, ch_nDetected{iChan}, 1));
            output(:,15) = num2cell(repmat(centerFreqFilter, ch_nDetected{iChan}, 1));
            
            output(:,16) = num2cell(ch_detectedLengthSamples{iChan}'/FrqOfSmpl);
            output(:,17) = num2cell(ch_detectedBeginSample{iChan}'/FrqOfSmpl);
            output(:,18) = num2cell(ch_detectedEndSample{iChan}'/FrqOfSmpl);
            output(:,19) = num2cell(ch_detectedPeaksSamples{iChan}'/FrqOfSmpl);
            output(:,20) = num2cell(ch_detectedTroughsSamples{iChan}'/FrqOfSmpl);
            output(:,21) = num2cell(ch_detectedEnvelopeMaxSamples{iChan}'/FrqOfSmpl);
            
            output(:,22) = num2cell((1:ch_nDetected{iChan})');
            output(:,23) = cellstr(epochs(:,1));
            output(:,24) = cellstr(epochs(:,2));
            output(:,25) = cellstr(epochs(:,3));
            
            output(:,26) = num2cell(ch_detectedSDofFilteredSignal{iChan}');
            output(:,27) = num2cell(ch_detectedMergeCount{iChan}');
            
            output(:,28) = num2cell(ch_detected_linear_regression_freq_slope{iChan});
            output(:,29) = num2cell(ch_detected_linear_regression_freq_offset{iChan});
            output(:,30) = num2cell(ch_detected_linear_regression_freq_R_squared{iChan});
            
            output(:,31) = num2cell(ch_contigSegment{iChan});
            
            %31
            %32
            %33
            
            
            for iLine=1:(size(output,1))
                fprintf(fide,'%i,',iData);
                fprintf(fide,'%s,',output{iLine,1});
                tempLengthSeconds = output{iLine,16};
                fprintf(fide,'%f,',tempLengthSeconds);
                fprintf(fide,'%e,',output{iLine,7});
                tempNtroughs = length(find(~isnan(ch_detectedSignalTroughsSamples{iChan}(iLine,:))));
                tempNpeaks = length(find(~isnan(ch_detectedSignalPeaksSamples{iChan}(iLine,:))));
                fprintf(fide,'%f,',((tempNtroughs + tempNpeaks)/2)/tempLengthSeconds );

                fprintf(fide,'%i,',output{iLine,2});
                fprintf(fide,'%i,',output{iLine,3});
                fprintf(fide,'%i,',output{iLine,4});
                fprintf(fide,'%i,',output{iLine,5});
                fprintf(fide,'%i,',output{iLine,6});
                fprintf(fide,'%e,',output{iLine,8});

                fprintf(fide,'%s,',output{iLine,9});
                fprintf(fide,'%s,',output{iLine,10});
                fprintf(fide,'%s,',output{iLine,11});
                fprintf(fide,'%e,',output{iLine,12});
                fprintf(fide,'%f,',output{iLine,13});
                fprintf(fide,'%f,',output{iLine,14});
                fprintf(fide,'%f,',output{iLine,15});
                

                fprintf(fide,'%f,',output{iLine,17});
                fprintf(fide,'%f,',output{iLine,18});
                fprintf(fide,'%f,',output{iLine,19});
                fprintf(fide,'%f,',output{iLine,20});
                fprintf(fide,'%f,',output{iLine,21});
                fprintf(fide,'%i,',output{iLine,22});
                fprintf(fide,'%s,',output{iLine,23});
                fprintf(fide,'%s,',output{iLine,24});
                fprintf(fide,'%s,',output{iLine,25});
                
                fprintf(fide,'%i,',output{iLine,31});
                
                fprintf(fide,'%e,',output{iLine,26});
                fprintf(fide,'%i,',output{iLine,27});
                fprintf(fide,'%f,',output{iLine,28});
                fprintf(fide,'%f,',output{iLine,29});
                fprintf(fide,'%f,',output{iLine,30});

                fprintf(fide,'%i,',tempNtroughs);
                fprintf(fide,'%i\n',tempNpeaks);
            end
            
            
            for iLine=1:(size(output,1))
                fprintf(fidp,'%i,',iData);
                fprintf(fidp,'%s,',ch);
                fprintf(fidp,'%i,',output{iLine,20});
                fprintf(fidp,'%f,',ch_detectedSignalPeaksSamples{iChan}(iLine,1:end-1)/FrqOfSmpl );
                fprintf(fidp,'%f\n',ch_detectedSignalPeaksSamples{iChan}(iLine,end)/FrqOfSmpl );
            end
            
            for iLine=1:(size(output,1))
                fprintf(fidt,'%i,',iData);
                fprintf(fidt,'%s,',ch);
                fprintf(fidt,'%i,',output{iLine,20});
                fprintf(fidt,'%f,',ch_detectedSignalTroughsSamples{iChan}(iLine,1:end-1)/FrqOfSmpl );
                fprintf(fidt,'%f\n',ch_detectedSignalTroughsSamples{iChan}(iLine,end)/FrqOfSmpl );
            end
            
            for iLine=1:(size(output,1))
                fprintf(fidpf,'%i,',iData);
                fprintf(fidpf,'%s,',ch);
                fprintf(fidpf,'%i,',output{iLine,20});
                fprintf(fidpf,'%f,',ch_detected_inst_freq_peaks{iChan}(iLine,1:end-1));
                fprintf(fidpf,'%f\n',ch_detected_inst_freq_peaks{iChan}(iLine,end));
            end
            
            for iLine=1:(size(output,1))
                fprintf(fidtf,'%i,',iData);
                fprintf(fidtf,'%s,',ch);
                fprintf(fidtf,'%i,',output{iLine,20});
                fprintf(fidtf,'%f,',ch_detected_inst_freq_troughs{iChan}(iLine,1:end-1));
                fprintf(fidtf,'%f\n',ch_detected_inst_freq_troughs{iChan}(iLine,end));
            end
            
            for iLine=1:(size(output,1))
                fprintf(fidpa,'%i,',iData);
                fprintf(fidpa,'%s,',ch);
                fprintf(fidpa,'%i,',output{iLine,20});
                fprintf(fidpa,'%e,',ch_detectedPeaksPotential{iChan}(iLine,1:end-1));
                fprintf(fidpa,'%e\n',ch_detectedPeaksPotential{iChan}(iLine,end));
            end
            
            for iLine=1:(size(output,1))
                fprintf(fidta,'%i,',iData);
                fprintf(fidta,'%s,',ch);
                fprintf(fidta,'%i,',output{iLine,20});
                fprintf(fidta,'%e,',ch_detectedTroughsPotential{iChan}(iLine,1:end-1));
                fprintf(fidta,'%e\n',ch_detectedTroughsPotential{iChan}(iLine,end));
            end
        end
        
        
    end
    fclose(fidc);
    fclose(fide);
    fclose(fidp);
    fclose(fidt);
    fclose(fidpf);
    fclose(fidtf);
    fclose(fidpa);
    fclose(fidta);
    
    data = [];%clear
    chData = [];%clear
    
end

%aggregate all results from datasets
temp_fidf_all = [];
temp_fidc_all = [];
temp_fide_all = [];
temp_fidp_all = [];
temp_fidt_all = [];
temp_fidpf_all = [];
temp_fidtf_all = [];
temp_fidpa_all = [];
temp_fidta_all = [];
delimiter = ',';
if ~strcmp(AggregationOfDatasetOutputsOfDetections,'no')
    
    fprintf('Aggregate results of all datasets\n');
    for iData = iDatas
        
        temp_fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        temp_fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
            temp_fide = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_events_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidp = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidt = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidpf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_freq_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidtf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_freq_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidpa = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_filtered_potential_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
            temp_fidta = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_filtered_potential_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        end
        
        if iData == iDatas(1)
            temp_fidf_all = temp_fidf;
            temp_fidc_all = temp_fidc;
            if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
                temp_fide_all = temp_fide;
                temp_fidp_all = temp_fidp;
                temp_fidt_all = temp_fidt;
                temp_fidpf_all = temp_fidpf;
                temp_fidtf_all = temp_fidtf;
                temp_fidpa_all = temp_fidpa;
                temp_fidta_all = temp_fidta;
            end
        else
            if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
                
                maxpeaks_prev = size(temp_fidp_all,2)-3;
                maxpeaks_next = size(temp_fidp,2)-3;
                maxpeaks = max(maxpeaks_prev,maxpeaks_next);
                nAddEmptyPeaks_prev = maxpeaks - maxpeaks_prev;
                nAddEmptyPeaks_next = maxpeaks - maxpeaks_next;
                
                NaN_prev = zeros(size(temp_fidp_all,1),1);
                NaN_next = zeros(size(temp_fidp,1),1);
                NaN_prev(:,:) = NaN;
                NaN_next(:,:) = NaN;
                
                
                for iAddNaN = 1:nAddEmptyPeaks_prev
                    NaN_prevds = mat2dataset(NaN_prev,'VarNames',{['extr_' num2str(maxpeaks_prev+iAddNaN)]});
                    temp_fidp_all = cat(2,temp_fidp_all,NaN_prevds);
                    temp_fidt_all = cat(2,temp_fidt_all,NaN_prevds);
                    NaN_prevds = mat2dataset(NaN_prev,'VarNames',{['extr_freq_' num2str(maxpeaks_prev+iAddNaN)]});
                    temp_fidpf_all = cat(2,temp_fidpf_all,NaN_prevds);
                    temp_fidtf_all = cat(2,temp_fidtf_all,NaN_prevds);
                    NaN_prevds = mat2dataset(NaN_prev,'VarNames',{['extr_pot_' num2str(maxpeaks_prev+iAddNaN)]});
                    
                    temp_fidpa_all = cat(2,temp_fidpa_all,NaN_prevds);
                    temp_fidta_all = cat(2,temp_fidta_all,NaN_prevds);
                end
                
                for iAddNaN = 1:nAddEmptyPeaks_next
                    NaN_nextds = mat2dataset(NaN_next,'VarNames',{['extr_' num2str(maxpeaks_next+iAddNaN)]});
                    temp_fidp = cat(2,temp_fidp,NaN_nextds);
                    temp_fidt = cat(2,temp_fidt,NaN_nextds);
                    NaN_nextds = mat2dataset(NaN_next,'VarNames',{['extr_freq_' num2str(maxpeaks_next+iAddNaN)]});
                    temp_fidpf = cat(2,temp_fidpf,NaN_nextds);
                    temp_fidtf = cat(2,temp_fidtf,NaN_nextds);
                    NaN_nextds = mat2dataset(NaN_next,'VarNames',{['extr_pot_' num2str(maxpeaks_next+iAddNaN)]});
                    temp_fidpa = cat(2,temp_fidpa,NaN_nextds);
                    temp_fidta = cat(2,temp_fidta,NaN_nextds);
                end
                
            end
            temp_fidf_all = cat(1,temp_fidf_all,temp_fidf);
            temp_fidc_all = cat(1,temp_fidc_all,temp_fidc);
            if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
                
                temp_fide_all = cat(1,temp_fide_all,temp_fide);
                temp_fidp_all = cat(1,temp_fidp_all,temp_fidp);
                temp_fidt_all = cat(1,temp_fidt_all,temp_fidt);
                temp_fidpf_all = cat(1,temp_fidpf_all,temp_fidpf);
                temp_fidtf_all = cat(1,temp_fidtf_all,temp_fidtf);
                temp_fidpa_all = cat(1,temp_fidpa_all,temp_fidpa);
                temp_fidta_all = cat(1,temp_fidta_all,temp_fidta);
            end
        end
        
    end
    export(temp_fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    export(temp_fidc_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
        export(temp_fide_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_events_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidp_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidt_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidpf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_freq_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidtf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_freq_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidpa_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_peaks_filtered_potential_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
        export(temp_fidta_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'spindle_troughs_filtered_potential_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    end
end

res_filters = temp_fidf_all;
res_channels = temp_fidc_all;
res_events = temp_fide_all;
res_peaks = temp_fidp_all;
res_troughs = temp_fidt_all;
res_peaks_freq = temp_fidpf_all;
res_troughs_freq = temp_fidtf_all;
res_peaks_filtered_potential = temp_fidpa_all;
res_troughs_filtered_potential = temp_fidta_all;


fprintf('SPD function finished\n');
toc
memtoc
end