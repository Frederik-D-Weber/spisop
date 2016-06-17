function [res_filters, res_channels, res_events] = spisop_sod_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% Slow oscillation and deltawave detection
% Copyright Frederik D. Weber

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
if strcmp(CenterFrequenciesFileName,'indcenterfreq');
    if exist([pathInputFolder filesep CenterFrequenciesFileName],'file') ~= 2
        error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep CenterFrequenciesFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
end

FrequencyBoundaryMethod = getParam('FrequencyBoundaryMethod',listOfParameters); %generalminmaxfreq or indcenterfreq

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz

PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz
PreDetectionLowPassFilterFreq_FpassRight = str2num(getParam('PreDetectionLowPassFilterFreq_FpassRight',listOfParameters));%in Hz

preCenterFreqFilterTo_FpassLeft = str2num(getParam('preCenterFreqFilterTo_FpassLeft',listOfParameters)); % in Hz
postCenterFreqFilterTo_FpassRight = str2num(getParam('postCenterFreqFilterTo_FpassRight',listOfParameters)); % in Hz

MinDetectionFreq = str2num(getParam('MinDetectionFreq',listOfParameters));%in seconds
MaxDetectionFreq = str2num(getParam('MaxDetectionFreq',listOfParameters));%in seconds


MinAbsoluteDownToUpPeakPotential = str2num(getParam('MinAbsoluteDownToUpPeakPotential',listOfParameters));
MaxAbsoluteDownToUpPeakPotential = str2num(getParam('MaxAbsoluteDownToUpPeakPotential',listOfParameters));
MinAbsoluteUpPeakPotential = str2num(getParam('MinAbsoluteUpPeakPotential',listOfParameters));
MaxAbsoluteDownPeakPotential = str2num(getParam('MaxAbsoluteDownPeakPotential',listOfParameters));

MeanFactorNegativePeak = str2num(getParam('MeanFactorNegativePeak',listOfParameters));
MeanFactorPeak2Peak = str2num(getParam('MeanFactorPeak2Peak',listOfParameters));

epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
%sleepStagesOfInterst = {'S3','S4'};
%sleepStagesOfInterst = {'SWS','S2'};
sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));

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
listOfCenterFrequencies = [];
if strcmp(FrequencyBoundaryMethod,'indcenterfreq');
    listOfCenterFrequencies = dlmread([pathInputFolder filesep CenterFrequenciesFileName],',')';
end

if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
    if strcmp(FrequencyBoundaryMethod,'indcenterfreq') && ~(all(size(listOfDatasetsPaths) == size(listOfCenterFrequencies')))
        error('files or number of Datasetspaths and CenterFrequnecies are invalid or do not aggree')
    end;
    error('files or number of Datasetspaths Hypnogramsfiles ChannelsOfInterest are invalid or do not aggree')
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end

if epochLength < (1/PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff)
    error(['Parameter epochLength ' num2str(epochLength) 's must not be greater in order to support the maximum of PreDetectionHighPassFilter_FpassLeft of ' num2str(PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff) ' Hz!'])
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

core_cfg.lpfilttype    = getParam('ft_cfg_lpfilttype',listOfCoreParameters);
core_cfg.lpfiltdir     = getParam('ft_cfg_lpfiltdir',listOfCoreParameters);
core_cfg.lpinstabilityfix = getParam('ft_cfg_lpinstabilityfix',listOfCoreParameters);

core_cfg.hpfilttype    = getParam('ft_cfg_hpfilttype',listOfCoreParameters);
core_cfg.hpfiltdir     = getParam('ft_cfg_hpfiltdir',listOfCoreParameters);
core_cfg.hpinstabilityfix = getParam('ft_cfg_hpinstabilityfix',listOfCoreParameters);

% if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
%     error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
% end

if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
    error(['low pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

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
AstopRight = str2num(getParam('AstopRight',listOfCoreParameters)); %Attenuation of right stop band (>FstopRight) ripples in db

% Apass_bp = Apass;
% AstopLeft_bp = AstopLeft;
% AstopRight_bp = AstopRight;

Apass_lp = Apass;
AstopRight_lp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;



% StopToPassTransitionWidth_bp = str2num(getParam('StopToPassTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
% PassToStopTransitionWidth_bp = str2num(getParam('PassToStopTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

PassToStopTransitionWidth_lp = str2num(getParam('PassToStopTransitionWidth_lp',listOfCoreParameters)); %for low pass filter frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

StopToPassTransitionWidth_hp = str2num(getParam('StopToPassTransitionWidth_hp',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2

StopToPassTransitionWidth_hp_predownsample = str2num(getParam('StopToPassTransitionWidth_hp_predownsample',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2



% UseFixedFilterOrder_bp = str2num(getParam('UseFixedFilterOrder_bp',listOfCoreParameters));

UseFixedFilterOrder_lp = (getParam('UseFixedFilterOrder_lp',listOfCoreParameters));

UseFixedFilterOrder_hp = (getParam('UseFixedFilterOrder_hp',listOfCoreParameters));

% UseTwoPassAttenuationCorrection_bp = (getParam('UseTwoPassAttenuationCorrection_bp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_lp = (getParam('UseTwoPassAttenuationCorrection_lp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_hp = (getParam('UseTwoPassAttenuationCorrection_hp',listOfCoreParameters));

% FilterOrder_bp = str2num(getParam('FilterOrder_bp',listOfCoreParameters));

FilterOrder_lp = str2num(getParam('FilterOrder_lp',listOfCoreParameters));

FilterOrder_hp = str2num(getParam('FilterOrder_hp',listOfCoreParameters));

MaximizeFilterOrderIfFixedFilterOrderIsUsed = str2num(getParam('MaximizeFilterOrderIfFixedFilterOrderIsUsed',listOfCoreParameters));

% useTwoPassFiltering_bp = 'no';

useTwoPassFiltering_lp = 'no';

useTwoPassFiltering_hp = 'no';

% if ~isempty(strfind(core_cfg.bpfiltdir,'two'))
%     useTwoPassFiltering_bp = 'yes';
% end

if ~isempty(strfind(core_cfg.lpfiltdir,'two'))
    useTwoPassFiltering_lp = 'yes';
end

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

if strcmp(UseFixedFilterOrder_lp,'yes') && strcmp(useTwoPassFiltering_lp,'yes') && (FilterOrder_lp > maxFilterOrder)
    error(['filter order for low pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_lp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_lp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end


% if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
%     error('band pass filter order must be an even number')
% end

if strcmp(UseFixedFilterOrder_lp,'yes') && logical(mod(FilterOrder_lp,2))
    error('low pass order must be an even number')
end

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

if strcmp(useTwoPassFiltering_lp,'yes') && strcmp(UseTwoPassAttenuationCorrection_lp,'yes')
    Apass_lp = Apass_lp/2;
    AstopRight_lp = AstopRight_lp/2;
end

if strcmp(useTwoPassFiltering_hp,'yes') && strcmp(UseTwoPassAttenuationCorrection_hp,'yes')
    Apass_hp = Apass_hp/2;
    AstopLeft_hp = AstopLeft_hp/2;
end


if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(core_cfg.hpfilttype,'FIRdesigned')
    error('UseFixedFilterOrder_hp not allowed for high pass filters of type FIRdesigned')
end

tic
memtic
fprintf('SOD function initialized\n');
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

    minFreq = -1;
    maxFreq = -1;
    
    if strcmp(FrequencyBoundaryMethod,'indcenterfreq');
        centerFreqFilter = listOfCenterFrequencies(iData);
        minFreq = centerFreqFilter - preCenterFreqFilterTo_FpassLeft;
        maxFreq = centerFreqFilter + postCenterFreqFilterTo_FpassRight;
    elseif strcmp(FrequencyBoundaryMethod,'generalminmaxfreq');
        minFreq = MinDetectionFreq;
        maxFreq = MaxDetectionFreq;
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
    cfg.hpfilter = 'yes';
    FpassLeft = PreDetectionHighPassFilter_FpassLeft_or_F3dBcutoff; %left pass frequency in Hz
    FstopLeft = FpassLeft - StopToPassTransitionWidth_hp; %left stop frequency in Hz
    
    usedFilterOrder_hp = NaN;
    hp_hdm = NaN;
    if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned')
        
        hp_d = [];
        hp_hd = [];
        if strcmp(UseFixedFilterOrder_hp,'yes')
            hp_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        else
            hp_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
        end
        fprintf('dataset %i: designing high pass filter for filtering SO band \n',iData);
        if strcmp(core_cfg.hpfilttype,'IIRdesigned')
            hp_hd = design(hp_d,'butter'); %isstable(hp_hd)
        elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
            hp_hd = design(hp_d,'equiripple','MinOrder', 'even');
        else
            error(['highpass filter type of ' core_cfg.hpfilttype ' unknown or not allowed'])
        end
        usedFilterOrder_hp = hp_hd.order;
        cfg.hpfilterdesign = hp_hd;
        hp_hdm = measure(hp_hd);
    end
    if strcmp(UseFixedFilterOrder_hp,'yes')
        cfg.hpfiltord     = FilterOrder_hp;
    end
    cfg.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function
    cfg.feedback = core_cfg.feedback;
    fprintf('dataset %i: reprocess and apply high pass filter for SO band\n',iData);
    data = ft_fw_preprocessing(cfg,data);
    
    
    
    
    
    cfg = [];
    cfg = core_cfg;
    cfg.lpfilter = 'yes';
    FpassRight = PreDetectionLowPassFilterFreq_FpassRight; %right pass frequency in Hz
    FstopRight = FpassRight + PassToStopTransitionWidth_lp; %right stop frequency in Hz
    usedFilterOrder_lp = NaN;
    lp_hdm = NaN;
    if strcmp(core_cfg.lpfilttype,'IIRdesigned') || strcmp(core_cfg.lpfilttype,'FIRdesigned')
        lp_d = [];
        lp_hd = [];
        fprintf('dataset %i: designing low pass filter for filtering SO band \n',iData);
        if strcmp(UseFixedFilterOrder_lp,'yes')
            lp_d = fdesign.lowpass('N,Fp,Fst',FilterOrder_lp,FpassRight,FstopRight,FrqOfSmpl);
            lp_hd = design(lp_d,'equiripple');
        else
            lp_d = fdesign.lowpass('Fp,Fst,Ap,Ast',FpassRight,FstopRight,Apass_lp,AstopRight_lp,FrqOfSmpl);
            lp_hd = design(lp_d,'equiripple','MinOrder', 'even');
        end
        usedFilterOrder_lp = lp_hd.order;
        cfg.lpfilterdesign = lp_hd;
        lp_hdm = measure(lp_hd);
    end
    if strcmp(UseFixedFilterOrder_lp,'yes')
        cfg.lpfiltord     = FilterOrder_lp;
    end
    cfg.lpfreq        = [FpassRight];%dummy values are overwritten by low level function
    cfg.feedback = core_cfg.feedback;
    fprintf('dataset %i: reprocess and apply low pass filter for SO band\n',iData);
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
    
    
    smplsMinDetectionLength  = round(FrqOfSmpl/maxFreq);
    smplsMaxDetectionLength  = round(FrqOfSmpl/minFreq);
    
    minFreqPostFreqBorderBufferLength = round(FrqOfSmpl/minFreq);
    
    
    nChannels = length(data.label);
    
    ch_detectedLengthSamples = [];
    ch_detectedBeginSample = [];
    ch_detectedEndSample = [];
    ch_detectedSignalMin = [];
    ch_detectedSignalMax = [];
    ch_detectedPeak2Peaks = [];
    ch_detectedPeaksSamples  = [];
    ch_detectedTroughsSamples  = [];
    ch_detectedMaxSlopes  = [];
    ch_detectedMaxSlopesSamples  = [];
    ch_detectedMaxDownSlopes  = [];
    ch_detectedMaxDownSlopesSamples  = [];
    ch_detectedTroughToZeroCrossingSlopes = [];
    ch_detectedTroughToZeroCrossingSlopesSamples = [];
    ch_detectedZeroCrossingToTroughSlopes = [];
    %ch_detectedZeroCrossingToTroughSlopesSamples = [];
    ch_detectedSDofFilteredSignal = [];
    ch_nDetected = [];
    ch_densityPerEpoch = [];
    ch_meanNegativePeak =  [];
    ch_meanPeak2Peak = [];
    
    for iChan = 1:nChannels
        %iChan = 1;
        
        fprintf('dataset %i: process channel %s\n',iData,data.label{iChan});
        %iChan = 1;
        cfg = [];
        cfg.channel = ft_channelselection(data.label{iChan}, data.label);
        cfg.feedback = core_cfg.feedback;
        chData = ft_selectdata(cfg,data);
        
        
        trl_detectedLengthSamples = [];
        trl_detectedBeginSample = [];
        trl_detectedEndSample = [];
        trl_detectedSignalMin = [];
        trl_detectedSignalMax = [];
        trl_detectedPeak2Peaks = [];
        trl_detectedPeaksSamples  = [];
        trl_detectedTroughsSamples  = [];
        trl_detectedMaxSlopes  = [];
        trl_detectedMaxSlopesSamples  = [];
        trl_detectedMaxDownSlopes  = [];
        trl_detectedMaxDownSlopesSamples  = [];
        trl_detectedTroughToZeroCrossingSlopes = [];
        trl_detectedTroughToZeroCrossingSlopesSamples = [];
        trl_detectedZeroCrossingToTroughSlopes = [];
        %trl_detectedZeroCrossingToTroughSlopesSamples = [];
        trl_detectedSDofFilteredSignal = [];
        trl_nDetected = 0;
        
        for iTr = 1:size(chData.trial,2)
            %iTr = 4;
            fprintf('dataset %i: channel %s, subpart %i, preselect events of appropriate zero-crossings\n',iData,data.label{iChan},iTr);
            
            rawDataSampleOffset = trlSampleBeginsAndEnds(iTr,1) - 1;
            
            frqBndPssSignal = chData.trial{iTr};
            frqBndPssSignal_hilbert = hilbert(frqBndPssSignal);

            
            %thresholdForDetection = ch_SDfrqBndPssSignal*factorSD;
            
            lengthSignal = length(frqBndPssSignal);
            
            [begins, ends] = getBeginsAndCorrespondingEndsIndicesBelowThreshold(frqBndPssSignal,0);
            
            firstSmplsMinFreqPostFreqBorderBufferLength = minFreqPostFreqBorderBufferLength;
            lastSmplsMinFreqPostFreqBorderBufferLength = lengthSignal - minFreqPostFreqBorderBufferLength;
            indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects filter
            
            begins = begins(indicesValidSamples);
            ends = ends(indicesValidSamples);
            
            nDetected = length(begins)-1;
            if (nDetected > 0)
                ends = begins(2:end)';
                begins = begins(1:end-1)';
                
                tempCandidatesLengths = ends - begins + 1;
                
                indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
                
                nDetected = length(indicesCandiates);
            end
            
            if (nDetected > 0)
                detectedLengthSamples = tempCandidatesLengths(indicesCandiates);
                detectedBeginSample = begins(indicesCandiates);
                detectedEndSample = ends(indicesCandiates);
                
                minPeakDistanceSamples = ceil(((1/(maxFreq)) * FrqOfSmpl)/2); % half of max freq of interest in samples
                
                detectedSignalMin = zeros(1,nDetected);
                detectedSignalMax = zeros(1,nDetected);
                detectedPeak2Peaks = zeros(1,nDetected);
                detectedPeaksSamples  = zeros(1,nDetected);
                detectedTroughsSamples  = zeros(1,nDetected);
                detectedMaxSlopes  = zeros(1,nDetected);
                detectedMaxSlopesSamples  = zeros(1,nDetected);
                detectedMaxDownSlopes  = zeros(1,nDetected);
                detectedMaxDownSlopesSamples  = zeros(1,nDetected);
                detectedTroughToZeroCrossingSlopes = zeros(1,nDetected);
                detectedTroughToZeroCrossingSlopesSamples = zeros(1,nDetected);
                detectedZeroCrossingToTroughSlopes  = zeros(1,nDetected);
                %detectedZeroCrossingToTroughSlopesSamples  = zeros(1,nDetected);
                
                detectedSDofFilteredSignal = zeros(1,nDetected);
                
                fprintf('dataset %i: channel %s, subpart %i, annotate events\n',iData,data.label{iChan},iTr);
                for iIterCand = 1:nDetected
                    %iIterCand = 1
                    currentRawDataSampleOffset = rawDataSampleOffset + detectedBeginSample(iIterCand) - 1;
                    candSignal = frqBndPssSignal(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));
                    candSignal_hilbert = frqBndPssSignal_hilbert(detectedBeginSample(iIterCand):detectedEndSample(iIterCand));

                    
                    tempCandSignalmaxSample = find(max(candSignal) == candSignal);
                    tempCandSignalminSample = find(min(candSignal(1:tempCandSignalmaxSample)) == candSignal(1:tempCandSignalmaxSample));
                    
                    
                    candSignalminSample =  currentRawDataSampleOffset + tempCandSignalminSample;
                    candSignalmaxSample = currentRawDataSampleOffset + tempCandSignalmaxSample;
                    
                    
                    candSignalMaxSlope = max(diff(candSignal(tempCandSignalminSample:tempCandSignalmaxSample)))*FrqOfSmpl;% in potential/s
                    candSignalMaxSlopeSample = find(diff(candSignal)*FrqOfSmpl == candSignalMaxSlope);
                    candSignalMaxSlopeSample = intersect(candSignalMaxSlopeSample,tempCandSignalminSample:tempCandSignalmaxSample);
                    candSignalMaxSlopeSample = candSignalMaxSlopeSample(find(abs(candSignal(candSignalMaxSlopeSample)) == min(abs(0-(candSignal(candSignalMaxSlopeSample)))),1,'first'));
                    
                    candSignalMaxSlopeSample =  currentRawDataSampleOffset + candSignalMaxSlopeSample;
                    
                    
                    
                    candSignalMaxDownSlope = min(diff(candSignal(1:tempCandSignalminSample)))*FrqOfSmpl;% in potential/s
                    if isempty(candSignalMaxDownSlope)
                        candSignalMaxDownSlope = NaN;
                        candSignalMaxDownSlopeSample = NaN;
                    else
                        candSignalMaxDownSlopeSample = find(diff(candSignal)*FrqOfSmpl == candSignalMaxDownSlope);
                        candSignalMaxDownSlopeSample = intersect(candSignalMaxDownSlopeSample,1:tempCandSignalminSample);
                        candSignalMaxDownSlopeSample = candSignalMaxDownSlopeSample(find(abs(candSignal(candSignalMaxDownSlopeSample)) == min(abs(0-(candSignal(candSignalMaxDownSlopeSample)))),1,'first'));
                        candSignalMaxDownSlopeSample =  currentRawDataSampleOffset + candSignalMaxDownSlopeSample;
                    end
                    
                    
                    
                    candSignalmin = candSignal(tempCandSignalminSample);
                    candSignalmax = candSignal(tempCandSignalmaxSample);
                    candPeak2Peak = candSignalmax - candSignalmin;
                    
                    tempCandSignalMinMaxZeroCrossingSample = tempCandSignalminSample + find( abs(0-candSignal(tempCandSignalminSample:tempCandSignalmaxSample)) == min( abs(0-candSignal(tempCandSignalminSample:tempCandSignalmaxSample)) ),1,'first' ) - 1;
                    candSignalMinMaxZeroCrossingSample = currentRawDataSampleOffset + tempCandSignalMinMaxZeroCrossingSample;
                    candSignalTroughToZeroCrossingSlope = abs(candSignalmin)/((tempCandSignalMinMaxZeroCrossingSample - tempCandSignalminSample)/FrqOfSmpl); % in potential/s
                    
                    
                    %tempCandSignalDownMinMaxZeroCrossingSample = 0 + find( abs(0-candSignal(1:tempCandSignalminSample)) == min( abs(0-candSignal(1:tempCandSignalminSample)) ),1,'first' ) - 1;
                    %candSignalDownMinMaxZeroCrossingSample = currentRawDataSampleOffset + tempCandSignalDownMinMaxZeroCrossingSample;
                    candSignalDownZeroCrossingToTroughSlope = -abs(candSignalmin)/((tempCandSignalminSample - 0)/FrqOfSmpl); % in potential/s
                    
                    
%                      %inst_freq = (diff(unwrap(angle(hilbert(candSignal))))/(2*pi))*FrqOfSmpl;
%                     temp_inst_freq = (diff(unwrap(angle(candSignal_hilbert)))/(2*pi))*FrqOfSmpl;
%                     temp_time = (1:(length(candSignal)))./FrqOfSmpl;
%                     temp_time_freq_regression = fitlm(temp_time(2:end),temp_inst_freq);                   
%                     
%                     temp_linear_regression_freq_slope = temp_time_freq_regression.Coefficients.Estimate(2);
%                     temp_linear_regression_freq_offset = temp_time_freq_regression.Coefficients.Estimate(1);
%                     temp_linear_regression_freq_R_squared = temp_time_freq_regression.Rsquared.Ordinary;
%                     
% %                     temp_inst_freq_troughs = temp_inst_freq(tempCandSignalTroughsSamples);
% %                     temp_inst_freq_peaks = temp_inst_freq(tempCandSignalPeaksSamples);
% 
%                     figure
%                     subplot(2,1,1);
%                     plot(temp_time(2:end),temp_inst_freq); hold on;
%                     plot(temp_time(2:end),temp_time_freq_regression.Coefficients.Estimate(2) * temp_time(2:end) + temp_time_freq_regression.Coefficients.Estimate(1))
%                     hold off;
%                     
%                     subplot(2,1,2);
%                     plot(temp_time(2:end),candSignal(2:end))
                    
%                                         figure
%                                         plot((currentRawDataSampleOffset+1):(currentRawDataSampleOffset+(ends(iCand)-begins(iCand))+1),candSignal)
%                                         hold on
%                                         %plot(candSignalPeaksSamples,candSignalPeaks,'g*')
%                                         %plot(begins(iCand):minPeakDistanceSamples:ends(iCand),0,'r*')
%                                         plot(candSignalminSample,candSignalmin,'b*',candSignalmaxSample,candSignalmax,'b*')
%                                         plot(candSignalTroughsSamples,candSignalTroughs,'r*',candSignalPeaksSamples,candSignalPeaks,'g*')
                     
                    detectedSignalMin(iIterCand) = candSignalmin;
                    detectedSignalMax(iIterCand) = candSignalmax;
                    detectedPeak2Peaks(iIterCand) = candPeak2Peak;
                    detectedPeaksSamples(iIterCand) = candSignalmaxSample;
                    detectedTroughsSamples(iIterCand) = candSignalminSample;
                    detectedMaxSlopes(iIterCand) = candSignalMaxSlope;
                    detectedMaxSlopesSamples(iIterCand) = candSignalMaxSlopeSample;
                    detectedMaxDownSlopes(iIterCand) = candSignalMaxDownSlope;
                    detectedMaxDownSlopesSamples(iIterCand) = candSignalMaxDownSlopeSample;
                    detectedTroughToZeroCrossingSlopes(iIterCand) = candSignalTroughToZeroCrossingSlope;
                    detectedTroughToZeroCrossingSlopesSamples(iIterCand) = candSignalMinMaxZeroCrossingSample;
                    detectedZeroCrossingToTroughSlopes(iIterCand) = candSignalDownZeroCrossingToTroughSlope;
                    %detectedZeroCrossingToTroughSlopesSamples(iIterCand) = candSignalDownMinMaxZeroCrossingSample;
                    %detectedSignalTroughsSamples(iIterCand,1:nCandSignalTroughs) = candSignalTroughsSamples;
                    %detectedSignalPeaksSamples(iIterCand,1:nCandSignalPeaks) = candSignalPeaksSamples;
                    detectedSDofFilteredSignal(iIterCand) = std(candSignal);
                    
                end
                
                trl_detectedLengthSamples = cat(2,trl_detectedLengthSamples,detectedLengthSamples);
                trl_detectedBeginSample = cat(2,trl_detectedBeginSample,detectedBeginSample + rawDataSampleOffset - 1);
                trl_detectedEndSample = cat(2,trl_detectedEndSample,detectedEndSample + rawDataSampleOffset - 1);
                
                trl_detectedSignalMin = cat(2,trl_detectedSignalMin,detectedSignalMin);
                trl_detectedSignalMax = cat(2,trl_detectedSignalMax,detectedSignalMax);
                trl_detectedPeak2Peaks = cat(2,trl_detectedPeak2Peaks,detectedPeak2Peaks);
                trl_detectedPeaksSamples = cat(2,trl_detectedPeaksSamples,detectedPeaksSamples);
                trl_detectedTroughsSamples  = cat(2,trl_detectedTroughsSamples,detectedTroughsSamples);
                trl_detectedMaxSlopes = cat(2,trl_detectedMaxSlopes ,detectedMaxSlopes);
                trl_detectedMaxSlopesSamples = cat(2,trl_detectedMaxSlopesSamples,detectedMaxSlopesSamples);
                trl_detectedMaxDownSlopes = cat(2,trl_detectedMaxDownSlopes ,detectedMaxDownSlopes);
                trl_detectedMaxDownSlopesSamples = cat(2,trl_detectedMaxDownSlopesSamples,detectedMaxDownSlopesSamples);
                trl_detectedTroughToZeroCrossingSlopes = cat(2,trl_detectedTroughToZeroCrossingSlopes,detectedTroughToZeroCrossingSlopes);
                trl_detectedTroughToZeroCrossingSlopesSamples = cat(2,trl_detectedTroughToZeroCrossingSlopesSamples,detectedTroughToZeroCrossingSlopesSamples);
                trl_detectedZeroCrossingToTroughSlopes = cat(2,trl_detectedZeroCrossingToTroughSlopes,detectedZeroCrossingToTroughSlopes);
                %trl_detectedZeroCrossingToTroughSlopesSamples = cat(2,trl_detectedZeroCrossingToTroughSlopesSamples,detectedZeroCrossingToTroughSlopesSamples);
                %trl_detectedSignalTroughsSamples = cat(1,trl_detectedSignalTroughsSamples,detectedSignalTroughsSamples);
                %trl_detectedSignalPeaksSamples = cat(1,trl_detectedSignalPeaksSamples,detectedSignalPeaksSamples);
                trl_detectedSDofFilteredSignal = cat(2,trl_detectedSDofFilteredSignal,detectedSDofFilteredSignal);
                
                trl_nDetected = trl_nDetected + nDetected;
            end
        end
        
        fprintf('dataset %i: channel %s, select events\n',iData,data.label{iChan});
        tempIndexAboveMeanThreshold_premean = find((trl_detectedPeak2Peaks >= MinAbsoluteDownToUpPeakPotential) ...
            & (trl_detectedPeak2Peaks <= MaxAbsoluteDownToUpPeakPotential) ...
            & (trl_detectedSignalMax >= MinAbsoluteUpPeakPotential) ...
            & (trl_detectedSignalMin <= MaxAbsoluteDownPeakPotential));
        
        ch_meanNegativePeak{iChan} =  mean(trl_detectedSignalMin(tempIndexAboveMeanThreshold_premean));
        ch_meanPeak2Peak{iChan} =  mean(trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold_premean));
        
        tempIndexAboveMeanThreshold_temp = find((trl_detectedTroughsSamples >= (ch_meanNegativePeak{iChan}*MeanFactorNegativePeak)) ...
            & (trl_detectedPeak2Peaks >= (ch_meanPeak2Peak{iChan}*MeanFactorPeak2Peak)));
        
        tempIndexAboveMeanThreshold = intersect(tempIndexAboveMeanThreshold_temp,tempIndexAboveMeanThreshold_premean);
        
        ch_detectedLengthSamples{iChan} = trl_detectedLengthSamples(tempIndexAboveMeanThreshold);
        ch_detectedBeginSample{iChan} = trl_detectedBeginSample(tempIndexAboveMeanThreshold);
        ch_detectedEndSample{iChan} = trl_detectedEndSample(tempIndexAboveMeanThreshold);
        ch_detectedSignalMin{iChan} = trl_detectedSignalMin(tempIndexAboveMeanThreshold);
        ch_detectedSignalMax{iChan} = trl_detectedSignalMax(tempIndexAboveMeanThreshold);
        ch_detectedPeak2Peaks{iChan} = trl_detectedPeak2Peaks(tempIndexAboveMeanThreshold);
        ch_detectedPeaksSamples{iChan} = trl_detectedPeaksSamples(tempIndexAboveMeanThreshold);
        ch_detectedTroughsSamples{iChan}  = trl_detectedTroughsSamples(tempIndexAboveMeanThreshold);
        ch_detectedMaxSlopes{iChan} = trl_detectedMaxSlopes(tempIndexAboveMeanThreshold);
        ch_detectedMaxSlopesSamples{iChan} = trl_detectedMaxSlopesSamples(tempIndexAboveMeanThreshold);
        ch_detectedMaxDownSlopes{iChan} = trl_detectedMaxDownSlopes(tempIndexAboveMeanThreshold);
        ch_detectedMaxDownSlopesSamples{iChan} = trl_detectedMaxDownSlopesSamples(tempIndexAboveMeanThreshold);
        
        ch_detectedTroughToZeroCrossingSlopes{iChan} = trl_detectedTroughToZeroCrossingSlopes(tempIndexAboveMeanThreshold);
        ch_detectedTroughToZeroCrossingSlopesSamples{iChan} = trl_detectedTroughToZeroCrossingSlopesSamples(tempIndexAboveMeanThreshold);
        ch_detectedZeroCrossingToTroughSlopes{iChan} = trl_detectedZeroCrossingToTroughSlopes(tempIndexAboveMeanThreshold);
        %ch_detectedZeroCrossingToTroughSlopesSamples{iChan} = trl_detectedZeroCrossingToTroughSlopesSamples(tempIndexAboveMeanThreshold);
        
        %ch_detectedSignalTroughsSamples{iChan} = trl_detectedSignalTroughsSamples;
        %ch_detectedSignalPeaksSamples{iChan} = trl_detectedSignalPeaksSamples;
        ch_detectedSDofFilteredSignal{iChan} = trl_detectedSDofFilteredSignal(tempIndexAboveMeanThreshold);
        
        ch_nDetected{iChan} = length(tempIndexAboveMeanThreshold);
        ch_densityPerEpoch{iChan} = length(tempIndexAboveMeanThreshold)/(lengthsAcrossROIsSeconds/epochLength);
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
    
    if ~(strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned'))
        
        usedFilterOrder_lp = NaN;
        lp_hdm.Fs = FrqOfSmpl;
        lp_hdm.Astop = NaN;
        lp_hdm.Fstop = NaN;
        lp_hdm.F6dB = NaN;
        lp_hdm.F3dB = FpassRight;
        lp_hdm.TransitionWidth = NaN;
        lp_hdm.Fpass = NaN;
        lp_hdm.Apass = NaN;
        
        if strcmp(core_cfg.lpfilttype,'but')
            if strcmp(UseFixedFilterOrder_lp,'yes')
                usedFilterOrder_lp = FilterOrder_lp;
            else
                usedFilterOrder_lp = 6;
            end
        end
        
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
        
        usedFilterOrder_hp = NaN;
        hp_hdm.Fs = FrqOfSmpl;
        hp_hdm.Astop = NaN;
        hp_hdm.Fstop = NaN;
        hp_hdm.F6dB = NaN;
        hp_hdm.F3dB = FpassLeft;
        hp_hdm.TransitionWidth = NaN;
        hp_hdm.Fpass = NaN;
        hp_hdm.Apass = NaN;
        
        if strcmp(core_cfg.hpfilttype,'but')
            if strcmp(UseFixedFilterOrder_hp,'yes')
                usedFilterOrder_hp = FilterOrder_hp;
                usedFilterOrder_hp_preDS = FilterOrder_hp;
            else
                usedFilterOrder_hp = 6;
                usedFilterOrder_hp_preDS = 6;
            end
        end
    end
    
    
    
    
    lp_f_type_detail = '';
    switch core_cfg.lpfilttype
        case 'but'
            lp_f_type_detail = 'IIR_Butterworth_ml_butter';
        case 'fir'
            lp_f_type_detail = 'FIR_window_Hamming_ml_fir1';
        case 'FIRdesigned'
            lp_f_type_detail = 'FIR_equiripple_signal_toolbox';
        case 'IIRdesigned'
            lp_f_type_detail = 'IIR_Butterworth_signal_toolbox';
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
    
    
    
    fidf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'so_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
    %write header
    fprintf(fidf,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'datasetnum','dataset',...
        'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB',...
        'hp_filter','hp_filter_type','hp_dir_and_passing','usedFilterOrder_hp','hp_Fs_Hz','hp_Astop_dB','hp_Fstop_Hz','hp_F6dB_Hz','hp_F3dB_Hz','hp_TransitionWidth_Hz','hp_Fpass_Hz','hp_Apass_dB',...
        'lp_filter','lp_filter_type','lp_dir_and_passing','usedFilterOrder_lp','lp_Fs_Hz','lp_Astop_dB','lp_Fstop_Hz','lp_F6dB_Hz','lp_F3dB_Hz','lp_TransitionWidth_Hz','lp_Fpass_Hz','lp_Apass_dB');
    %write content
    fprintf(fidf,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' '\n'],...
        iData,datasetsPath,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp,hp_hdm.Fs,hp_hdm.Astop,hp_hdm.Fstop,hp_hdm.F6dB,hp_hdm.F3dB,hp_hdm.TransitionWidth,hp_hdm.Fpass,hp_hdm.Apass,...
        core_cfg.lpfilttype,lp_f_type_detail,core_cfg.lpfiltdir,usedFilterOrder_lp,lp_hdm.Fs,lp_hdm.Astop,lp_hdm.Fstop,lp_hdm.F6dB,lp_hdm.F3dB,lp_hdm.TransitionWidth,lp_hdm.Fpass,lp_hdm.Apass);
    
    
    fclose(fidf);
    
    
    
    
    %open output files
    fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'so_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
    fide = fopen([pathOutputFolder filesep ouputFilesPrefixString 'so_events_' 'datanum_' num2str(iData) '.csv'],'wt');
    %fidp = fopen([pathOutputFolder filesep 'spindle_peaks_' 'datanum_' num2str(iData) '.csv'],'wt');
    %fidt = fopen([pathOutputFolder filesep 'spindle_troughs_' 'datanum_' num2str(iData) '.csv'],'wt');
    
    %write header of ouptufiles
    fprintf(fidc,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','datasetnum','channel','count','density_per_epoch',...
        'mean_duration_seconds','mean_amplitude_peak2trough_potential',...
        'mean_slope_to_trough_min_potential_per_second','mean_slope_zerocrossing_to_trough_potential_per_second',...
        'mean_slope_trough_to_up_max_potential_per_second','mean_slope_trough_to_zerocrossing_potential_per_second',...
        'mean_frequency_by_duration','mean_frequency_by_trough_to_peak_latency',...
        'epoch_length_seconds','lengths_ROI_seconds',...
        'used_meanNegPeak_potential','used_meanP2T_potential',...
        'mean_SD_of_filtered_signal','mean_negative_peak_potential','mean_positive_peak_potential');
    
    fprintf(fide,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        'datasetnum','channel','duration_seconds','amplitude_peak2trough_max',...
        'slope_to_trough_min_potential_per_second','slope_zerocrossing_to_trough_potential_per_second',...
        'slope_trough_to_up_max_potential_per_second','slope_trough_to_zerocrossing_potential_per_second',...
        'frequency_by_duration','frequency_by_trough_to_peak_latency','duration_samples','sample_begin','sample_end','sample_peak_max','sample_trough_max'...
        ,'dataset','hypnogram','used_stages_for_detection','used_meanNegPeak','used_meanP2T','seconds_begin','seconds_end','seconds_peak_max','seconds_trough_max','id_within_channel',...
        'sample_slope_to_trough_min',...
        'sample_slope_trough_to_up_max','sample_zerocrossig_of_trough_to_zerocrossig_slope_trough_to_up_max',...
        'seconds_slope_to_trough_min',...
        'seconds_slope_trough_to_up_max','seconds_zerocrossig_of_trough_to_zerocrossig_slope_trough_to_up_max'...
        ,'negative_peak_potential','positive_peak_potential','stage','stage_alt','stage_alt2','SD_of_filtered_signal');
    
    
    for iChan = 1:nChannels
        
        ch = data.label{iChan};
        
        epochs = {};
        for iDet = 1:length(ch_detectedTroughsSamples{iChan})
            tempSample = ch_detectedTroughsSamples{iChan}(iDet);
            epochs(iDet,:) = [hypnStages(((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples)),1) ...
                hypnStages(((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples)),2) ...
                hypnStages(((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples)),3)];
        end;
        
        
        fprintf(fidc,'%i,',iData);
        fprintf(fidc,'%s,',ch);
        fprintf(fidc,'%i,',ch_nDetected{iChan});
        fprintf(fidc,'%f,',ch_densityPerEpoch{iChan});
        
        tempLengthMeanSeconds = mean(ch_detectedLengthSamples{iChan}'/FrqOfSmpl);
        fprintf(fidc,'%e,',tempLengthMeanSeconds);
        fprintf(fidc,'%e,',mean(ch_detectedPeak2Peaks{iChan}));
        
        fprintf(fidc,'%e,',mean(ch_detectedMaxDownSlopes{iChan}));
        fprintf(fidc,'%e,',mean(ch_detectedZeroCrossingToTroughSlopes{iChan}));
        
        fprintf(fidc,'%e,',mean(ch_detectedMaxSlopes{iChan}));
        fprintf(fidc,'%e,',mean(ch_detectedTroughToZeroCrossingSlopes{iChan}));
        
        fprintf(fidc,'%f,', (1/tempLengthMeanSeconds));
        tempTroughPeakLengthMean = mean((ch_detectedPeaksSamples{iChan} - ch_detectedTroughsSamples{iChan} )'/FrqOfSmpl);
        fprintf(fidc,'%f,', (1/(tempTroughPeakLengthMean*2)));
        
        fprintf(fidc,'%f,',epochLength);
        fprintf(fidc,'%f,',lengthsAcrossROIsSeconds);
        
        fprintf(fidc,'%e,',ch_meanNegativePeak{iChan});
        fprintf(fidc,'%e,',ch_meanPeak2Peak{iChan});
        
        
        
        fprintf(fidc,'%e,',mean(ch_detectedSDofFilteredSignal{iChan}));
        fprintf(fidc,'%e,',mean(ch_detectedSignalMin{iChan}));
        fprintf(fidc,'%e\n',mean(ch_detectedSignalMax{iChan}));
        

        
        
        if ch_nDetected{iChan} > 0
            output = cell(ch_nDetected{iChan},11);
            
            output(:,1) = cellstr(repmat(ch, ch_nDetected{iChan}, 1));
            output(:,2) = num2cell(ch_detectedLengthSamples{iChan}');
            output(:,3) = num2cell(ch_detectedBeginSample{iChan}');
            output(:,4) = num2cell(ch_detectedEndSample{iChan}');
            output(:,5) = num2cell(ch_detectedPeaksSamples{iChan}');
            output(:,6) = num2cell(ch_detectedTroughsSamples{iChan}');
            
            output(:,7) = num2cell(ch_detectedPeak2Peaks{iChan});
            
            output(:,8) = cellstr(repmat(datasetsPath, ch_nDetected{iChan}, 1));
            output(:,9) = cellstr(repmat(hypnogramPath, ch_nDetected{iChan}, 1));
            output(:,10) = cellstr(repmat(strjoin(sleepStagesOfInterest,' '), ch_nDetected{iChan}, 1));
            output(:,11) = num2cell(repmat(ch_meanNegativePeak{iChan}, ch_nDetected{iChan}, 1));
            output(:,12) = num2cell(repmat(ch_meanPeak2Peak{iChan}, ch_nDetected{iChan}, 1));
            
            output(:,13) = num2cell(ch_detectedLengthSamples{iChan}'/FrqOfSmpl);
            output(:,14) = num2cell(ch_detectedBeginSample{iChan}'/FrqOfSmpl);
            output(:,15) = num2cell(ch_detectedEndSample{iChan}'/FrqOfSmpl);
            output(:,16) = num2cell(ch_detectedPeaksSamples{iChan}'/FrqOfSmpl);
            output(:,17) = num2cell(ch_detectedTroughsSamples{iChan}'/FrqOfSmpl);
            
            output(:,18) = num2cell((1:ch_nDetected{iChan})');
            
            output(:,19) = num2cell(ch_detectedMaxSlopes{iChan}');
            output(:,20) = num2cell(ch_detectedTroughToZeroCrossingSlopes{iChan}');
            output(:,21) = num2cell(ch_detectedMaxSlopesSamples{iChan}');
            output(:,22) = num2cell(ch_detectedTroughToZeroCrossingSlopesSamples{iChan}');
            output(:,23) = num2cell(ch_detectedMaxSlopesSamples{iChan}'/FrqOfSmpl);
            output(:,24) = num2cell(ch_detectedTroughToZeroCrossingSlopesSamples{iChan}'/FrqOfSmpl);
            
            output(:,25) = num2cell(ch_detectedSignalMin{iChan}');
            output(:,26) = num2cell(ch_detectedSignalMax{iChan}');
            output(:,27) = cellstr(epochs(:,1));
            output(:,28) = cellstr(epochs(:,2));
            output(:,29) = cellstr(epochs(:,3));
            output(:,30) = num2cell(ch_detectedSDofFilteredSignal{iChan});
            
            output(:,31) = num2cell(ch_detectedMaxDownSlopes{iChan}');
            output(:,32) = num2cell(ch_detectedZeroCrossingToTroughSlopes{iChan}');
            output(:,33) = num2cell(ch_detectedMaxDownSlopesSamples{iChan}');
            output(:,34) = num2cell(ch_detectedMaxDownSlopesSamples{iChan}'/FrqOfSmpl);
            
            
            
            for iLine=1:(size(output,1))
                fprintf(fide,'%i,',iData);
                fprintf(fide,'%s,',output{iLine,1});
                tempLengthSeconds = output{iLine,13};
                fprintf(fide,'%f,',tempLengthSeconds);
                fprintf(fide,'%e,',output{iLine,7});
                fprintf(fide,'%e,',output{iLine,31});
                fprintf(fide,'%e,',output{iLine,32});
                fprintf(fide,'%e,',output{iLine,19});
                fprintf(fide,'%e,',output{iLine,20});
                tempTroughPeakLength = (output{iLine,16} - output{iLine,17});
                fprintf(fide,'%f,', (1/tempLengthSeconds));
                fprintf(fide,'%f,', (1/(tempTroughPeakLength*2)));
                fprintf(fide,'%i,',output{iLine,2});
                fprintf(fide,'%i,',output{iLine,3});
                fprintf(fide,'%i,',output{iLine,4});
                fprintf(fide,'%i,',output{iLine,5});
                fprintf(fide,'%i,',output{iLine,6});
                fprintf(fide,'%s,',output{iLine,8});
                fprintf(fide,'%s,',output{iLine,9});
                fprintf(fide,'%s,',output{iLine,10});
                fprintf(fide,'%e,',output{iLine,11});
                fprintf(fide,'%e,',output{iLine,12});
                
                
                fprintf(fide,'%f,',output{iLine,14});
                fprintf(fide,'%f,',output{iLine,15});
                
                fprintf(fide,'%f,',output{iLine,16});
                fprintf(fide,'%f,',output{iLine,17});
                fprintf(fide,'%i,',output{iLine,18});

                fprintf(fide,'%i,',output{iLine,33});                
                fprintf(fide,'%i,',output{iLine,21});
                fprintf(fide,'%i,',output{iLine,22});
                
                fprintf(fide,'%f,',output{iLine,34});
                fprintf(fide,'%f,',output{iLine,23});
                fprintf(fide,'%f,',output{iLine,24});
                
                fprintf(fide,'%e,',output{iLine,25});
                fprintf(fide,'%e,',output{iLine,26});
                fprintf(fide,'%s,',output{iLine,27});
                fprintf(fide,'%s,',output{iLine,28});
                fprintf(fide,'%s,',output{iLine,29});
                fprintf(fide,'%e\n',output{iLine,30});
            end
        end
        
        
        
        
    end
    fclose(fidc);
    fclose(fide);
    data = [];%clear
    chData = [];%clear
    
end


%aggregate all results from datasets
temp_fidf_all = [];
temp_fidc_all = [];
temp_fide_all = [];
delimiter = ',';

if ~strcmp(AggregationOfDatasetOutputsOfDetections,'no')
    fprintf('Aggregate results of all datasets\n');
    for iData = iDatas
        
        
        temp_fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'so_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        temp_fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'so_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
            temp_fide = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'so_events_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        end
        if iData == iDatas(1)
            temp_fidf_all = temp_fidf;
            temp_fidc_all = temp_fidc;
            if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
                temp_fide_all = temp_fide;
            end
        else
            temp_fidf_all = cat(1,temp_fidf_all,temp_fidf);
            temp_fidc_all = cat(1,temp_fidc_all,temp_fidc);
            if strcmp(AggregationOfDatasetOutputsOfDetections,'full')
                temp_fide_all = cat(1,temp_fide_all,temp_fide);
            end
        end
        
    end
    export(temp_fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'so_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    export(temp_fidc_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'so_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    export(temp_fide_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'so_events_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    
end

res_filters = temp_fidf_all;
res_channels = temp_fidc_all;
res_events = temp_fide_all;


fprintf('SOD function finished\n');

toc
memtoc
end