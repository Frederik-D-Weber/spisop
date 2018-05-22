function [res_summary, res_channels, res_events] = spisop_remsMaAd_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% REMs and related phenomena detection
% Copyright Frederik D. Weber

functionName = 'remsMaAd';
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


FrqOfSmplWished = 100;

% EOG_low_pass_filter_freq = 5; % Hz
% EOG_high_pass_filter_freq = 0.5; % Hz
% basic_deviation_per_second_threshold = 261; %?V/s
% relaxed_deviation_per_second_threshold = 165; %?V/s

DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);

ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
%AVGoverChannels = getParam('AVGoverChannels',listOfParameters);

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


ExcludeMAepochs = getParam('ExcludeMAepochs',listOfParameters);%
RemoveSignalThatIsNotREM = getParam('RemoveSignalThatIsNotREM',listOfParameters);%

RelaxedThresholdFactorOfBasicThreshold = 0.66;
BasicThresholdFactorOfEOGnoise = 0.7;

try 
    RelaxedThresholdFactorOfBasicThreshold = str2num(getParam('RelaxedThresholdFactorOfBasicThreshold',listOfParameters)); %
catch e
end

try 
    BasicThresholdFactorOfEOGnoise = str2num(getParam('BasicThresholdFactorOfEOGnoise',listOfParameters)); %
catch e
end

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz


epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds

%sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));

%FrqOfSmplWished = str2num(getParam('FrqOfSmplWished',listOfParameters));%samples per second / Hz

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

core_cfg.lpfilttype    = getParam('ft_cfg_lpfilttype',listOfCoreParameters);
core_cfg.lpfiltdir     = getParam('ft_cfg_lpfiltdir',listOfCoreParameters);
core_cfg.lpinstabilityfix = getParam('ft_cfg_lpinstabilityfix',listOfCoreParameters);

core_cfg.hpfilttype    = getParam('ft_cfg_hpfilttype',listOfCoreParameters);
core_cfg.hpfiltdir     = getParam('ft_cfg_hpfiltdir',listOfCoreParameters);
core_cfg.hpinstabilityfix = getParam('ft_cfg_hpinstabilityfix',listOfCoreParameters);


if (strcmp(core_cfg.bpfilttype,'FIRdesigned') || strcmp(core_cfg.bpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.bpinstabilityfix,'no'))
    error(['band pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

if (strcmp(core_cfg.lpfilttype,'FIRdesigned') || strcmp(core_cfg.lpfilttype,'IIRdesigned')) && (~strcmp(core_cfg.lpinstabilityfix,'no'))
    error(['low pass filter instability fix not supported for FIRdesigned or IIRdesigned'])
end

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

Apass_lp = Apass;
AstopRight_lp = AstopRight;

Apass_hp = Apass;
AstopLeft_hp = AstopLeft;


StopToPassTransitionWidth_bp = str2num(getParam('StopToPassTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
PassToStopTransitionWidth_bp = str2num(getParam('PassToStopTransitionWidth_bp',listOfCoreParameters)); %frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

PassToStopTransitionWidth_lp = str2num(getParam('PassToStopTransitionWidth_lp',listOfCoreParameters)); %for low pass filter frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25

StopToPassTransitionWidth_hp = str2num(getParam('StopToPassTransitionWidth_hp',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2


StopToPassTransitionWidth_hp_predownsample = str2num(getParam('StopToPassTransitionWidth_hp_predownsample',listOfCoreParameters)); %for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2



UseFixedFilterOrder_bp = (getParam('UseFixedFilterOrder_bp',listOfCoreParameters));

UseFixedFilterOrder_lp = (getParam('UseFixedFilterOrder_lp',listOfCoreParameters));

UseFixedFilterOrder_hp = (getParam('UseFixedFilterOrder_hp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_bp = (getParam('UseTwoPassAttenuationCorrection_bp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_lp = (getParam('UseTwoPassAttenuationCorrection_lp',listOfCoreParameters));

UseTwoPassAttenuationCorrection_hp = (getParam('UseTwoPassAttenuationCorrection_hp',listOfCoreParameters));

FilterOrder_bp = str2num(getParam('FilterOrder_bp',listOfCoreParameters));

FilterOrder_lp = str2num(getParam('FilterOrder_lp',listOfCoreParameters));

FilterOrder_hp = str2num(getParam('FilterOrder_hp',listOfCoreParameters));

MaximizeFilterOrderIfFixedFilterOrderIsUsed = str2num(getParam('MaximizeFilterOrderIfFixedFilterOrderIsUsed',listOfCoreParameters));

useTwoPassFiltering_bp = 'no';

useTwoPassFiltering_lp = 'no';

useTwoPassFiltering_hp = 'no';

if ~isempty(strfind(core_cfg.bpfiltdir,'two'))
    useTwoPassFiltering_bp = 'yes';
end

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

if strcmp(UseFixedFilterOrder_bp,'yes') && strcmp(useTwoPassFiltering_bp,'yes') && (FilterOrder_bp > maxFilterOrder)
    error(['filter order for band pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_bp < maxFilterOrder) && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_bp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_lp,'yes') && strcmp(useTwoPassFiltering_lp,'yes') && (FilterOrder_lp > maxFilterOrder)
    error(['filter order for low pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_lp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_lp = maxFilterOrder;
end

if strcmp(UseFixedFilterOrder_hp,'yes') && strcmp(useTwoPassFiltering_hp,'yes') && (FilterOrder_hp > maxFilterOrder)
    error(['filter order for high pass not optimal, use a maximal filter order three times of samples in one epoch, i.e. filter order of ' num2str(maxFilterOrder)])
elseif (FilterOrder_hp < maxFilterOrder)  && strcmp(MaximizeFilterOrderIfFixedFilterOrderIsUsed,'yes')
    FilterOrder_hp = maxFilterOrder;
end


if strcmp(UseFixedFilterOrder_bp,'yes') && logical(mod(FilterOrder_bp,2))
    error('band pass filter order must be an even number')
end

if strcmp(UseFixedFilterOrder_lp,'yes') && logical(mod(FilterOrder_lp,2))
    error('low pass order must be an even number')
end

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
fprintf('REMsMaAd function initialized\n');
conseciDatas = 1:length(iDatas);
for conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    
    %iData = 1
    FrqOfSmplWishedPar = FrqOfSmplWished;
    
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    
    channelsOfInterest = listOfChannelsOfInterest(iData,:);
    channelsOfInterest = channelsOfInterest(~(cellfun(@isempty,channelsOfInterest)));
    channel_label_EOG_left = channelsOfInterest{1};
    channel_label_EOG_right = channelsOfInterest{2};
    
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
    
    % [roiBegins, roiEnds] = getROIsByHypnogram(hypnogramPath,epochLengthSamples,sleepStagesOfInterest);
    
    % if length(roiEnds) < 1
    %     error(['no ROI in data left for analysis']);
    % end
    
    if (signalOffsetSamples ~= 0)
        signalOffsetSeconds = signalOffsetSamples/preDownsampleFreq;
        %roiBegins = roiBegins + signalOffsetSamples;
        %roiEnds = roiEnds + signalOffsetSamples;
    end
    
    %     indexLastIncludedROIinData = length(roiBegins);
    %     nSampleLength = -1;
    %     if strcmp(IgnoreDataSetHeader,'no')
    %         nSampleLength = hdr.nSamples*hdr.nTrials + hdr.nSamplesPre;
    %         if (roiEnds(end) > nSampleLength)
    %             nMissingSamples  = roiEnds(end) - nSampleLength;
    %             warning ([ num2str(nMissingSamples) ' scored (hypnogram) samples (' num2str(nMissingSamples/preDownsampleFreq) ' seconds) missing in the dataset ' num2str(iData)]);
    %
    %             indexLastIncludedROIinData = find(roiBegins < nSampleLength,1,'last');
    %             roiBegins = roiBegins(1:indexLastIncludedROIinData);
    %             roiEnds = roiEnds(1:indexLastIncludedROIinData);
    %
    %             if (roiEnds(end) > nSampleLength)
    %                 roiEnds(indexLastIncludedROIinData) = nSampleLength;
    %             end;
    %
    %         end
    %     end
    
    %     if length(roiEnds) < 1
    %         error(['no ROI in data left for analysis']);
    %     end
    
    
    if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
        
        cfg = [];
        %cfg = core_cfg;
        cfg.feedback = core_cfg.feedback;
        cfg.precision = core_cfg.precision;
        %cfg.roiBegins = roiBegins;
        %cfg.roiEnds = roiEnds;
        %cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
        cfg.feedback = core_cfg.feedback;
        %cfg = ft_definetrial(cfg);
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
        %cfg.roiBegins = roiBegins;
        %cfg.roiEnds = roiEnds;
        %cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
        cfg.feedback = core_cfg.feedback;
        %cfg = ft_definetrial(cfg);
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        cfg.dataset = datasetsPath;
        data = ft_fw_preprocessing(cfg);
    end
    
    %     if strcmp(AVGoverChannels,'yes')
    %         for iTr = 1:size(data.trial,2)
    %             data.trial{iTr} = mean(data.trial{iTr},1);
    %
    %         end;
    %         data.label = {'meanOverChannels'};
    %     end
    
    
    %FrqOfSmplWished = 100;
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
    
    
    
    
    %signalOffsetSamples = 100
    %signalOffsetSeconds = 1;
    if (signalOffsetSamples ~= 0) %&& ~strcmp(DoEpochData,'yes')
        signalOffsetSamples_new = round(signalOffsetSeconds*FrqOfSmpl);
        if (signalOffsetSamples_new ~= 0)
            cfg = [];
            if (signalOffsetSamples_new > 0) && (signalOffsetSamples_new+1+epochLengthSamples < size(data.trial{1},2))
                cfg.begsample = signalOffsetSamples_new+1;
                cfg.endsample = numel(data.time{1});
                
                data = ft_redefinetrial(cfg,data);
                cfg = [];
                cfg.offset = signalOffsetSamples_new+1;
                data = ft_redefinetrial(cfg,data);
                data.sampleinfo = [1 numel(data.time{1})];
                warning(['IMPORTANT: for dataset ' num2str(iData) ': JUST CUT ' num2str(signalOffsetSamples_new) ' samples at the beginning corresponding to : ' num2str(signalOffsetSamples_new/FrqOfSmpl) ' seconds because of hypnogram/data offset!']);
            elseif signalOffsetSamples_new < 0
                for iTrTr = 1:numel(data.trial)
                    cfg.padtype = 'zero';
                    data.trial{iTrTr} = ft_preproc_padding(data.trial{iTrTr}, cfg.padtype, -signalOffsetSamples_new, 0);
                end
                data.time{1} = (0:(size(data.trial{1},2)-1))/data.fsample;
                data.sampleinfo = [1 numel(data.time{1})];
                warning(['IMPORTANT: for dataset ' num2str(iData) ': JUST ADDED ' num2str(signalOffsetSamples_new) ' samples at the beginning corresponding to : ' num2str(signalOffsetSamples_new/FrqOfSmpl) ' seconds because of hypnogram/data offset!']);
            end
        end
    end
    
    
    
    
    
    
    %%% Alorythm by Marek_Adamczyk %%%
    
    
    
    usedFixedSamplingRate = 100;
    epochLengthSamples = epochLength*usedFixedSamplingRate;
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    
    hip = hypn(:,1)';
    hipMA = hypn(:,2)';

    if strcmp(ExcludeMAepochs,'yes')
        hip((hipMA > 0) && (hip == 5)) = 51;        
    end
    
    keep_EpochsBeginsSamples = hypnEpochsBeginsSamples((hip == 5));
    keep_EpochsEndsSamples = hypnEpochsEndsSamples((hip == 5));
    
    keep_samples = zeros(1,size(data.trial{1},2));
    for iKeep = 1:numel(keep_EpochsBeginsSamples)
        keep_samples(:,keep_EpochsBeginsSamples(iKeep):keep_EpochsEndsSamples(iKeep)) = 1;
    end
    
     
    if strcmp(RemoveSignalThatIsNotREM,'yes')
        cut_EpochsBeginsSamples = hypnEpochsBeginsSamples(~(hip == 5));
        cut_EpochsEndsSamples = hypnEpochsEndsSamples(~(hip == 5));
        for iCut = 1:numel(cut_EpochsBeginsSamples)
            data.trial{1}(:,cut_EpochsBeginsSamples(iCut):cut_EpochsEndsSamples(iCut)) = 0;
        end
    end
    
    
    if isempty(find(hip == 5, 1))
        
        error(['No usable REM sleep in dataset number ' num2str(iData) ''])
        
        remDens = zeros(1,length(hip));
        meanREMSpeed = zeros(1,length(hip));
        meanREMWidth = zeros(1,length(hip));
        exactlyMarkedREMsVctSampling = usedFixedSamplingRate;
        exactlyMarkedREMs = zeros(length(hip)*epochLength*exactlyMarkedREMsVctSampling, 1);
        
    else
        
        
        sampling = FrqOfSmpl;
        
        % %           [hip, signal, REMsAddress, REMsData] = takeREMDataNotTwins(edfsPath, currFileName, false, false, 'nieee', sigInEDF);
        %
        %            %[hip, signal] = takeREMDataNotTwins(edfsPath, currFileName, false, false, 'nieee', sigInEDF);
        %            signal = takeEOGsig(edfFileAddress, sigInEDF);
        
        cfg = [];
        cfg.feedback = core_cfg.feedback;
        cfg.channel = channel_label_EOG_left;
        data_EOG_left = ft_selectdata(cfg,data);
        cfg.channel = channel_label_EOG_right;
        data_EOG_right = ft_selectdata(cfg,data);
        
        data = [];
        
        signal = [data_EOG_left.trial{1}' data_EOG_right.trial{1}'];
        
        data_EOG_right = [];
        data_EOG_left = [];
        
        
        fprintf('dataset %i: excluding artifacts\n',iData);
        eogClass =    eogClinicClassifyFlexThr(signal(:, 1), signal(:, 2), sampling, hip, false);%looks for artifacts!
        fprintf('dataset %i: detect REMs\n',iData);
        [remDens, exactlyMarkedREMs, meanREMSpeed, meanREMWidth, EOG1, EOG2, EOG1_lp, EOG2_lp, EOG1_lphp, EOG2_lphp, EOG1_lp_dev, EOG2_lp_dev, EOG1_lphp_dev, EOG2_lphp_dev] ...
            = paperSciRemDetectorAdjustableThreshold(signal(:, 1), signal(:, 2), sampling, epochLength, eogClass, false, 100000000000, RelaxedThresholdFactorOfBasicThreshold, BasicThresholdFactorOfEOGnoise);
        
        EOG1_lp_dev = [0; EOG1_lp_dev];
        EOG2_lp_dev = [0; EOG2_lp_dev];
        EOG1_lphp_dev = [0; EOG1_lphp_dev];
        EOG2_lphp_dev = [0; EOG2_lphp_dev];
        eogClass = [];
        signal = [];
        
        [begins, ends] = getBeginsAndCorrespondingEndsIndicesAboveThreshold(exactlyMarkedREMs*2,1);
        
        
         
         begins = begins(1:numel(ends));
         
         if isempty(begins)
             indicesValid = [];
         else
             indicesValid = find(keep_samples(begins) == 1 | (keep_samples(ends) == 1)); %consider border effects filter
         end
         begins = begins(indicesValid);
         ends = ends(indicesValid);
        
         nDetected = length(ends);
        %         indicesValidSamples = find((begins >= firstSmplsMinFreqPostFreqBorderBufferLength) & (ends <= lastSmplsMinFreqPostFreqBorderBufferLength)); %consider border effects filter
        %
        %         begins = begins(indicesValidSamples);
        %         ends = ends(indicesValidSamples);
        %
        
        fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_channels_' 'datanum_' num2str(iData) '.csv'],'wt');
        fide = fopen([pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_events_' 'datanum_' num2str(iData) '.csv'],'wt');
        
        %write header of ouptufiles
        fprintf(fidc,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            'datasetnum','channel',...
            'count','density_per_REM_epoch','density_per_REM_without_MA_epoch',...
            'mean_duration_seconds','mean_speed_uV_per_second',...
            'basic_threshold_factor_of_EOG_noise','relaxed_threshold_factor_of_basic_threshold');
        
        
        fprintf(fide,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
            'datasetnum','channel',...
            'duration_seconds',...
            'speed_uV_per_second_seconds',...
            'dataset','hypnogram','used_stages_for_detection',...
            'seconds_begin','seconds_end','seconds_center',...
            'stage','stage_alt','stage_alt2');
        
        fprintf(fidc,'%i,',iData);
        fprintf(fidc,'%s,',[channel_label_EOG_left '+' channel_label_EOG_right]);
        fprintf(fidc,'%i,',nDetected);
        
        nEpochs_REM = sum(hip == 5)+sum(hip == 51);
        nEpochs_REM_withoutMA = sum(hip == 5)-sum((hipMA == 1) & (hip == 5));
        fprintf(fidc,'%f,',nDetected/nEpochs_REM);
        fprintf(fidc,'%f,',nDetected/nEpochs_REM_withoutMA);
        
       
        if (nDetected > 0)
            
           
            detectedLengthSamples = ends - begins + 1;
            
            %             indicesCandiates = find((tempCandidatesLengths >= smplsMinDetectionLength) & (tempCandidatesLengths <= smplsMaxDetectionLength));
            %
            %             nDetected = length(indicesCandiates);
            
            detected = [];
            detected.id = (1:nDetected)';
            detected.seconds_begin = begins(1:nDetected)/usedFixedSamplingRate;
            detected.seconds_center = (begins(1:nDetected)+floor(detectedLengthSamples/2))/usedFixedSamplingRate;
            detected.seconds_end = ends(1:nDetected)/usedFixedSamplingRate;
            detected.seconds_duration = detectedLengthSamples(1:nDetected)/usedFixedSamplingRate;
            detected.REMSpeed_uV_per_second = zeros(nDetected,1);
            epochs = {};
            for iDetected = 1:nDetected
                detected.REMSpeed_uV_per_second(iDetected) = mean(abs(EOG1_lp_dev(begins:ends))) + mean(abs(EOG2_lp_dev(begins:ends)));
                
                detected.REMSpeed_uV_per_second(iDetected) = detected.REMSpeed_uV_per_second(iDetected)/detected.seconds_duration(iDetected);
                
                tempSample = begins(iDetected);
                tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample <= hypnEpochsEndsSamples));
                if ~any(tempInd)
                    tempInd = ((hypnEpochsBeginsSamples <= tempSample+1) & (tempSample <= hypnEpochsEndsSamples));
                end
                if ~any(tempInd)
                    tempInd = ((hypnEpochsBeginsSamples <= tempSample) & (tempSample-1 <= hypnEpochsEndsSamples));
                end
                epochs(iDetected,:) = [hypnStages(tempInd,1) ...
                    hypnStages(tempInd,2) ...
                    hypnStages(tempInd,3)];
            end
            
            detected.stage = epochs(:,1);
            detected.stage_alt = epochs(:,2);
            detected.stage_alt2 = epochs(:,3);
            epochs = [];
            
            detected.channel_label_EOG_left = cellstr(repmat(channel_label_EOG_left, nDetected, 1));
            detected.channel_label_EOG_right = cellstr(repmat(channel_label_EOG_right, nDetected, 1));
            detected.channel = cellstr(repmat([channel_label_EOG_left '+' channel_label_EOG_right], nDetected, 1));
            
            
            fprintf(fidc,'%f,',mean(detected.seconds_duration));
            fprintf(fidc,'%f,',mean(detected.REMSpeed_uV_per_second));
            fprintf(fidc,'%f,',BasicThresholdFactorOfEOGnoise);
            fprintf(fidc,'%f',RelaxedThresholdFactorOfBasicThreshold);

            %fprintf(fidc,'\n');
            
            
            
            for iDetected = 1:nDetected
                fprintf(fide,'%i,',iData);
                fprintf(fide,'%s,',[channel_label_EOG_left '+' channel_label_EOG_right]);
                
                fprintf(fide,'%f,',detected.seconds_duration(iDetected));
                fprintf(fide,'%f,',detected.REMSpeed_uV_per_second(iDetected));
                
                
                fprintf(fide,'%s,',datasetsPath);
                fprintf(fide,'%s,',hypnogramPath);
                fprintf(fide,'%s,','REM');
                
                fprintf(fide,'%f,',detected.seconds_begin(iDetected));
                fprintf(fide,'%f,',detected.seconds_end(iDetected));
                fprintf(fide,'%f,',detected.seconds_center(iDetected));
                
                fprintf(fide,'%s,',detected.stage{iDetected});
                fprintf(fide,'%s,',detected.stage_alt{iDetected});
                fprintf(fide,'%s',detected.stage_alt2{iDetected});
                if iDetected ~= nDetected
                    fprintf(fide,'\n');
                end
                
            end
        else
            
            fprintf(fidc,'%f,',NaN);
            fprintf(fidc,'%f,',NaN);
            fprintf(fidc,'%f,',BasicThresholdFactorOfEOGnoise);
            fprintf(fidc,'%f',RelaxedThresholdFactorOfBasicThreshold);
            %fprintf(fidc,'\n');
            
        end
        
        fclose(fidc);
        fclose(fide);
        
        
        
    end
    
    % change computed REM density in matlab format into easy readable csv
    % format
    
    
    path_file_humanreadable = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_density_speed_width_per_epoch_' 'datanum_' num2str(iData) '.csv'];
    saveREMsForPeople(path_file_humanreadable, remDens, meanREMSpeed./(meanREMWidth./usedFixedSamplingRate), meanREMWidth./usedFixedSamplingRate, hip, hipMA);
    
    path_file_exactlyMarkedREMs = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_exactlyMarkedREMs_' 'datanum_' num2str(iData) '.mat'];
    path_file_standardREMdens = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_standardREMdens_' 'datanum_' num2str(iData) '.mat'];
    
   % save(path_file_exactlyMarkedREMs,'exactlyMarkedREMs');
   % save(path_file_standardREMdens,'remDens');
    
    MarkedREMsSampling = usedFixedSamplingRate;
    
    [standardREMdensity1stCycle funName REMepochsNb1stCycle] = ...
        firstCycleREM_densFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    REM3Count1stCycle   = firstCycleMiniEp3WithREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1Count1stCycle =   firstCycleMiniEp1WithREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    AllREMsCount1stCycle =  firstCycleAllREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    allREMepochsNb =         allNightREMepochs(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REMlatency =                giveREMlatency(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    standardREMdensityAllNight = allNightREM_densFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REM3Count           = allNightMiniEp3WithREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1Count   =         allNightMiniEp1WithREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    AllREMsCount =                allNightAllREMsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    REMsCountInBurst =      allNightREMsInMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    burstsNumber =     allNightAllREMsNbOfMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REMsInBurstsPrc =    allNightREMsPrcInMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    
    REM1CountInBurst =       allNightREM1InMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1burstsNumber =     allNightREM1NbOfMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    REM1InBurstsPrc =     allNightREM1PrcInMidBurstsFromExactVct(hip, epochLength, exactlyMarkedREMs, MarkedREMsSampling);
    
    
    fid = fopen([pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_' 'datanum_' num2str(iData) '.txt'],'wt');
    
    fprintf(fid, '---------- First cycle data ---------\n\n');
    fprintf(fid, 'REM epochs number: %d\n', REMepochsNb1stCycle);
    fprintf(fid, '\n++++++ REM density +++++++\n\n');
    
    fprintf(fid, '%2.1f sec miniep: %2.2f\n', epochLength/10, standardREMdensity1stCycle);
    fprintf(fid, '1 sec miniep: %2.2f\n', REM1Count1stCycle/REMepochsNb1stCycle);
    fprintf(fid, 'all REMs    : %2.2f\n', AllREMsCount1stCycle/REMepochsNb1stCycle);
    fprintf(fid, '\n++++++ REM activity ++++++\n\n');
    fprintf(fid, '%2.1f sec miniep: %d\n', epochLength/10, REM3Count1stCycle);
    fprintf(fid, '1 sec miniep: %d\n', REM1Count1stCycle);
    fprintf(fid, 'all REMs    : %d\n', AllREMsCount1stCycle);
    
    fprintf(fid, '\n---------- All night data ---------\n\n');
    fprintf(fid, 'REM latency epochs number: %d\n', REMlatency);
    fprintf(fid, 'REM epochs number: %d\n', allREMepochsNb);
    fprintf(fid, '\n++++++ REM density +++++++\n\n');
    fprintf(fid, '%2.1f sec miniep: %2.2f\n', epochLength/10, standardREMdensityAllNight);
    fprintf(fid, '1 sec miniep: %2.2f\n', REM1Count/allREMepochsNb);
    fprintf(fid, 'all REMs    : %2.2f\n', AllREMsCount/allREMepochsNb);
    fprintf(fid, '\n++++++ REM activity ++++++\n\n');
    fprintf(fid, '%2.1f sec miniep: %d\n', epochLength/10, REM3Count);
    fprintf(fid, '1 sec miniep: %d\n', REM1Count);
    fprintf(fid, 'all REMs    : %d\n', AllREMsCount);
    fprintf(fid, '\n++++++  REM bursts  ++++++\n\n');
    fprintf(fid, '>>> all REMs <<<\n');
    fprintf(fid, 'activity: %d\n', REMsCountInBurst);
    fprintf(fid, 'density : %2.2f\n', REMsCountInBurst/allREMepochsNb);
    fprintf(fid, 'REM bursts number: %d\n', burstsNumber);
    fprintf(fid, 'Percentage of REMs in bursts: %2.2f%%\n', REMsInBurstsPrc);
    fprintf(fid, 'Nb of REMs in average burst : %2.2f\n', REMsCountInBurst/burstsNumber);
    
    fprintf(fid, '\n>>> 1 sec miniep <<<\n');
    fprintf(fid, 'activity: %d\n', REM1CountInBurst);
    fprintf(fid, 'density : %2.2f\n', REM1CountInBurst/allREMepochsNb);
    fprintf(fid, 'REM bursts number: %d\n', REM1burstsNumber);
    fprintf(fid, 'Percentage of REMs in bursts: %2.2f%%\n', 100*REM1InBurstsPrc);
    fprintf(fid, 'Nb of REMs in average burst : %2.2f\n', REM1CountInBurst/REM1burstsNumber);
    
    fclose(fid);
    
    path_to_summary_file = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_summarize_' 'datanum_' num2str(iData) '.csv'];
    fid2 = fopen(path_to_summary_file,'wt');
    
    
    fprintf(fid2, strcat('datasetnum,REMepochsNb1stCycle,standardREMdensity1stCycle,1sREMdensity1stCycle,allREMdensity1stCycle,', ...
        'standardREMactivity1stCycle,1sREMactivity1stCycle,allREMactivity1stCycle,', ...
        'REMlatencyEpNb,REMepochsNb,standardREMdensityAllNight,1sREMdensityAllNight,allREMdensityAllNight,', ...
        'standardREMactivityAllNight,1sREMactivityAllNight,allREMactivityAllNight,', ...
        'REMactivityInBurstAllNight,REMdensityInBurstAllNight,burstsNumberAllNight,', ...
        'REMsInBurstsPrc,avgREMsInBurst\n'));
    
    
    fprintf(fid2, '%d,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f\n', ...
        iData, REMepochsNb1stCycle, standardREMdensity1stCycle, REM1Count1stCycle/REMepochsNb1stCycle, ...
        AllREMsCount1stCycle/REMepochsNb1stCycle, REM3Count1stCycle, REM1Count1stCycle, AllREMsCount1stCycle, ...
        REMlatency, allREMepochsNb, standardREMdensityAllNight, REM1Count/allREMepochsNb, AllREMsCount/allREMepochsNb, ...
        REM3Count, REM1Count, AllREMsCount, REMsCountInBurst, REMsCountInBurst/allREMepochsNb, ...
        burstsNumber, REMsInBurstsPrc, REMsCountInBurst/burstsNumber);
    
    fclose(fid2);
    
    
    EOG1 = [];EOG2 = []; EOG1_lp = []; EOG2_lp = []; EOG1_lphp = []; EOG2_lphp = []; EOG1_lp_dev = []; EOG2_lp_dev = []; EOG1_lphp_dev = []; EOG2_lphp_dev = [];
    
    %%% END Alorythm by Marek_Adamczyk %%%
    
end % if regexpi(currFileName, '.edf')
%end





%aggregate all results from datasets
temp_fidc_all = [];
temp_fide_all = [];
temp_fids_all = [];

delimiter = ',';
if ~strcmp(AggregationOfDatasetOutputsOfDetections,'no')
    
    fprintf('Aggregate results of all datasets\n');
    for iData = iDatas
        
        temp_fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_channels_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        temp_fide = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_events_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
        
        path_to_summary_file = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_summarize_' 'datanum_' num2str(iData) '.csv'];
        
        temp_fids = dataset('File',path_to_summary_file,'Delimiter',delimiter);
        
        if iData == iDatas(1)
            temp_fidc_all = temp_fidc;
            temp_fide_all = temp_fide;
            temp_fids_all = temp_fids;
            
        else
            
            temp_fidc_all = cat(1,temp_fidc_all,temp_fidc);
            temp_fide_all = cat(1,temp_fide_all,temp_fide);
            temp_fids_all = cat(1,temp_fids_all,temp_fids);
            
        end
        
    end
    
    export(temp_fidc_all,'File',[pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_channels_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    export(temp_fide_all,'File',[pathOutputFolder filesep ouputFilesPrefixString 'remsMaAd_events_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
    
    path_to_summary_file_all = [pathOutputFolder filesep ouputFilesPrefixString 'REMsMaAd_summarize_' 'datanum_' 'all_recent' '.csv'];
    export(temp_fids_all,'file',path_to_summary_file_all,'Delimiter',delimiter);
    
end

res_summary = temp_fids_all;
res_channels = temp_fidc_all;
res_events = temp_fide_all;

fprintf('REMsMaAd function finished\n');
toc
memtoc
end
