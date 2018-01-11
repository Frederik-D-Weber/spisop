function [res_filters] = spisop_preprocout_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% preprocess files and output them
% Copyright Frederik D. Weber



DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);
ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
AVGoverChannels = getParam('AVGoverChannels',listOfParameters);

CutDataAtEndHypnogram = 'no'; %default
try
CutDataAtEndHypnogram = getParam('CutDataAtEndHypnogram',listOfParameters);% either yes or no default no
catch e
    
end
if ~(strcmp(CutDataAtEndHypnogram,'yes') || strcmp(CutDataAtEndHypnogram,'no'))
    error(['CutDataAtEndHypnogram ' CutDataAtEndHypnogram ' in parameters is unknown, use either yes or no, default no'])
end


if exist([pathInputFolder filesep DataSetPathsFileName],'file') ~= 2
    error(['DataSetPathsFileName file ' [pathInputFolder filesep DataSetPathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
if ~strcmp(IgnoreDataSetHeader,'yes')
    if exist([pathInputFolder filesep DataSetHeaderPathsFileName],'file') ~= 2
        error(['DataSetHeaderPathsFileName file ' [pathInputFolder filesep DataSetHeaderPathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
end

if exist([pathInputFolder filesep ChannelsOfInterestFileName],'file') ~= 2
    error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep ChannelsOfInterestFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end


PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz



DoEpochData = getParam('DoEpochData',listOfParameters);%If the data should be epoched like given in hypnograms and sleep stages of interest
if strcmp(DoEpochData,'yes')
    if exist([pathInputFolder filesep HypnogramsFileName],'file') ~= 2
        error(['HypnogramsFileName file ' [pathInputFolder filesep HypnogramsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
end
DefaultOutputUnit = getParam('DefaultOutputUnit',listOfParameters);
OutputDataformat = getParam('OutputDataformat',listOfParameters);
IncludePostiveMarkerAtBeginning = getParam('IncludePostiveMarkerAtBeginning',listOfParameters);


epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
%sleepStagesOfInterst = {'S3','S4'};
%sleepStagesOfInterst = {'SWS','S2'};
sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));


FrqOfSmplWished = str2num(getParam('FrqOfSmplWished',listOfParameters));%samples per second / Hz
FrqOfSmplWishedPreRedefine = str2num(getParam('FrqOfSmplWishedPreRedefine',listOfParameters));%samples per second / Hz

if FrqOfSmplWishedPreRedefine < FrqOfSmplWished
    error('FrqOfSmplWishedPreRedefine need to be same or bigger than FrqOfSmplWished')
end

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
listOfHypnogramPaths = {};
if strcmp(DoEpochData,'yes') || strcmp(CutDataAtEndHypnogram,'yes')
    listOfHypnogramPaths = read_mixed_csv([pathInputFolder filesep HypnogramsFileName],',');
end
listOfChannelsOfInterest = read_mixed_csv([pathInputFolder filesep ChannelsOfInterestFileName],',');

if strcmp(DoEpochData,'yes') || strcmp(CutDataAtEndHypnogram,'yes')
    if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
        error('files or number of Datasetspaths Hypnogramsfiles ChannelsOfInterest are invalid or do not aggree')
    end
else
    if ~((size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
        error('files or number of Datasetspaths nd ChannelsOfInterest are invalid or do not aggree')
    end
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end



% if epochLength < (1/MinDetectionFrequency_FpassLeft)
%     error(['Parameter epochLength ' num2str(epochLength) 's must not be greater in order to support the maximum of MinDetectionFrequency_FpassLeft of ' num2str(MinDetectionFrequency_FpassLeft) ' Hz!'])
% end

pretestHeaderForPersistentSampleFrequencies(IgnoreDataSetHeader,iDatas,listOfDatasetHeaderPaths,listOfChannelsOfInterest,FrqOfSmplWished);
pretestHeaderForPersistentSampleFrequencies(IgnoreDataSetHeader,iDatas,listOfDatasetHeaderPaths,listOfChannelsOfInterest,FrqOfSmplWishedPreRedefine);



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

ApplyFilterSettings = getParam('ApplyFilterSettings',listOfParameters);%either yes or no
FiltersSettingsDefinitionsFileName = getParam('FiltersSettingsDefinitionsFileName',listOfParameters);
listOfFilterSettingsFiles = {};
if strcmp(ApplyFilterSettings,'yes')
    if exist([pathInputFolder filesep FiltersSettingsDefinitionsFileName],'file') ~= 2
        error(['FiltersSettingsDefinitionsFileName file ' [pathInputFolder filesep FiltersSettingsDefinitionsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
    end
    listOfFilterSettingsFiles = read_mixed_csv([pathInputFolder filesep FiltersSettingsDefinitionsFileName],',');
    if ~(all(size(listOfDatasetsPaths) == size(listOfFilterSettingsFiles)))
        error('files or number of Datasetspaths listOfFilterSettingsFiles are invalid or do not aggree')
    end
    
    for iDefFiles = 1:size(listOfFilterSettingsFiles,1)
        if exist([pathInputFolder filesep listOfFilterSettingsFiles{iDefFiles}],'file') ~= 2
            error(['The fitler settings file listed in FiltersSettingsDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' does not exist'])
        end
        
        try
            readtable([pathInputFolder filesep listOfFilterSettingsFiles{iDefFiles}],'Delimiter',',');
        catch err
            error(['The fitler settings file listed in FiltersSettingsDefinitionsFileName for dataset number ' num2str(iDefFiles)  ' is not readable'])
        end
    end
    
end

DoWriteData = getParam('DoWriteData',listOfParameters);%either yes or no



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


if ~strcmp(core_cfg.bpfilttype,'FIRdesigned')
    error(['filter type for band pass not supported, only FIRdesigned allowed'])
end

if ~strcmp(core_cfg.lpfilttype,'FIRdesigned')
    error(['filter type for low pass not supported, only FIRdesigned allowed'])
end

if ~(strcmp(core_cfg.hpfilttype,'FIRdesigned') || strcmp(core_cfg.hpfilttype,'IIRdesigned') )
    error(['filter type for band pass not supported, only FIRdesigned or IIRdesigned allowed'])
end

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


% % Note that a one- or two-pass filter has consequences for the
% % strength of the filter, i.e. a two-pass filter with the same filter
% % order will attenuate the signal twice as strong.
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
fprintf('PREPROCOUT function initialized\n');
conseciDatas = 1:length(iDatas);
parfor conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 2
    
    FrqOfSmplWishedPar = FrqOfSmplWished;
    FrqOfSmplWishedParPreRedefine = FrqOfSmplWishedPreRedefine;
    
    datasetsPath = listOfDatasetsPaths{iData};
    if strcmp(DoEpochData,'yes') || strcmp(CutDataAtEndHypnogram,'yes')
        hypnogramPath = listOfHypnogramPaths{iData};
    end
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
    
    if strcmp(DoEpochData,'yes')
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
    end
    
    
  
    
    %     for second read in of hypnogramm after downsampling
    %     if (signalOffsetSamples ~= 0)
    %         signalOffsetSamples_downsampled = floor(signalOffsetSeconds*FrqOfSmpl);
    %         roiBegins = roiBegins + signalOffsetSamples_downsampled;
    %         roiEnds = roiEnds + signalOffsetSamples_downsampled;
    %     end
    
    
    if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
        
        cfg = [];
        %cfg = core_cfg;
        cfg.feedback = core_cfg.feedback;
        cfg.precision = core_cfg.precision;
        cfg.feedback = core_cfg.feedback;
        if strcmp(DoEpochData,'yes')
            cfg.roiBegins = roiBegins;
            cfg.roiEnds = roiEnds;
            cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
            cfg = ft_definetrial(cfg);
        end
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        cfg.dataset = datasetsPath;
        cfg.channel = 'all';
        fprintf('dataset %i: preprocess and pre filter data\n',iData);
        data = ft_fw_preprocessing(cfg);
        
        %FrqOfSmplWished = 400;
        
        if (FrqOfSmplWishedParPreRedefine < data.fsample)
            fprintf('dataset %i: resample data from %i to %i Hz\n',iData,data.fsample,FrqOfSmplWishedParPreRedefine);
            data_sel = {};
            temp_all_chans = data.label;
            for iChan = 1:length(temp_all_chans)
                cfg_sel = [];
                cfg_sel.channel = temp_all_chans(iChan);
                data_sel{iChan} = ft_selectdata(cfg_sel, data);
                cfg = [];
                cfg.resamplefs = FrqOfSmplWishedParPreRedefine;%frequency at which the data will be resampled (default = 256 Hz)
                cfg.detrend = 'no';
                cfg.feedback = core_cfg.feedback;
                data_sel{iChan} = ft_resampledata(cfg,data_sel{iChan});
                
                %cfg_sel2 = [];
                %cfg_sel2.channel = data.label(~strcmp(data.label,temp_all_chans(1:iChan)));
                %data = ft_selectdata(cfg_sel2,data);
                if iChan < length(temp_all_chans)
                    for iTr = 1:length(data.trial)
                        %data.trial{iTr}(1,:) = [];
                        data.trial{iTr} = data.trial{iTr}(2:end,:);
                    end
                    data.label(1) = [];
                end
            end
            data = [];
            cfg_append = [];
            data = ft_appenddata(cfg_append, data_sel{:});
        end
        
        
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
            montageTableChanNames = dataset('File',[pathInputFolder filesep linearDeviationMontageFile],'Delimiter',DelimiterLinearDeviationMontage,'ReadVarNames',false,'ReadObsNames',true);
            oldChanNames = {};
            for iCol = 1:length(get(montageTable,'VarNames'))
                oldChanNames{iCol} = montageTableChanNames{1,iCol};
            end
            
            montage = [];
            %montage.labelorg = get(montageTable,'VarNames');
            montage.labelorg = oldChanNames;
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
        cfg.feedback = core_cfg.feedback;
        if strcmp(DoEpochData,'yes')
            cfg.roiBegins = roiBegins;
            cfg.roiEnds = roiEnds;
            cfg.trialfun = 'trialfun_spd_ROIs'; %The cfg.trialfun option is a string containing the name of a function that you wrote yourself and that ft_definetrial will call.
            cfg = ft_definetrial(cfg);
        end
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
    

    if (signalMultiplicator ~= 1)
        data = ft_fw_factorMultiplicationOnSignal(data,'trial',signalMultiplicator);
    end
    
    FrqOfSmpl = data.fsample;%data.hdr.Fs;%samples per second / Hz
    
    %init the filter order variables
    usedFilterOrder_lp = -1;
    usedFilterOrder_hp = -1;
    usedFilterOrder_bp = -1;
    
    used_specified_filtering_once_bp = false;
    used_specified_filtering_once_hp = false;
    used_specified_filtering_once_lp = false;
    
    if strcmp(ApplyFilterSettings,'yes')
        fileFilterSettings = listOfFilterSettingsFiles{iData};
        curr_channel_settings_table = readtable([pathInputFolder filesep fileFilterSettings],'Delimiter',',');
        
        
        %                 [dummy_chan curr_indexOrder_setting] = ismember(data.label,curr_channel_settings_table.channel_label);
        curr_chanOrder = curr_channel_settings_table.channel_order;
        [dummy_chanOrder curr_chanIndexOrder] = sort(curr_chanOrder);
        
        %                 datreoder = data.trial{1};
        %                 datreoder = datreoder(curr_chanIndexOrder,:);
        %                 data.trial{1} = datreoder;
        %                 data.label = data.label(curr_chanIndexOrder);
        
        
        
        fprintf('dataset %i: apply filtering to data\n',iData);
        

        
        data_filt = {};
        iChanCount = 1;
        for iChannelbyOrder = curr_chanIndexOrder'
            %iChannelbyOrder = 2
            
            curr_channel_label = curr_channel_settings_table.channel_label(iChannelbyOrder);
            if ~all(ismember(curr_channel_label,data.label))
                continue
            end
            cfg = [];
            cfg.channel = curr_channel_label;
            data_filt{iChanCount} = ft_selectdata(cfg,data);
            
            curr_filterdefs = curr_channel_settings_table.filter_definitions(iChannelbyOrder);
            curr_filterdefs = strtrim(curr_filterdefs);
            curr_filterdefs = strsplit(char(curr_filterdefs));
            
            
            if (isempty(curr_filterdefs))
                curr_filterdefs = 'no';
            end
            
            curr_filterdefs_filterPos = 1;
            
            curr_bracket_opened = 0;
            %curr_operator_applied = false;
            curr_operator_wait_for_second_term = {};
            %curr_operator = '';
            data_store = {};%stack
            used_specified_filtering_once_bp = false;
            used_specified_filtering_once_hp = false;
            used_specified_filtering_once_lp = false;
            while curr_filterdefs_filterPos <= length(curr_filterdefs)
                curr_filter = char(curr_filterdefs(curr_filterdefs_filterPos));
                FpassLeft = -1;
                FpassRight = -1;
                MultFactor = 1;
                SmoothSeconds = -1;
                conversion_state_success = true;
                useSpecifiedFilterOrder = false;
                overwriteFilterType = false;
                try
                    switch curr_filter
                        case 'normtomean'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'min'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'max'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'ztransform'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'smooth'
                            [SmoothSeconds conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                        case 'envpeaks'
                            [minPeakDistanceSeconds conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;
                        case 'rect'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'env'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case 'mult'
                            [MultFactor conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2;                            
                        case {'bp', 'waveletband'}
                            [FpassLeft conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            [FpassRight conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                            temp_additional_values = 0;
                            if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+3)));
                                if conversion_state_success_potaddval
                                      useSpecifiedFilterOrder = true;
                                      filterOrder = round(pot_additional_value);
                                      temp_additional_values = temp_additional_values + 1; 
                                      
                                      if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+4) && strcmp(curr_filter,'bp')
                                          pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+4));
                                          if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                                overwriteFilterType = true;
                                                used_specified_filtering_once_bp = true;
                                                filterType = pot_filterType;
                                                temp_additional_values = temp_additional_values + 1;
                                          end
                                      end
                                      if strcmp(curr_filter,'waveletband')
                                          used_specified_filtering_once_bp = true;
                                      end
                                end
                            end
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 3 + temp_additional_values;
                        case 'hp'
                            [FpassLeft conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            temp_additional_values = 0;
                            if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                                [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                                if conversion_state_success_potaddval 
                                      useSpecifiedFilterOrder = true;
                                      filterOrder = round(pot_additional_value);
                                      temp_additional_values = temp_additional_values + 1;
                                      if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                          pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                          if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                                overwriteFilterType = true;
                                                used_specified_filtering_once_hp = true;
                                                filterType = pot_filterType;
                                                temp_additional_values = temp_additional_values + 1;
                                          end
                                      end
                                end
                            end
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                        case 'lp'
                            [FpassRight conversion_state_success] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+1)));
                            temp_additional_values = 0;
                            if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+2)
                                [pot_additional_value conversion_state_success_potaddval] = str2num(char(curr_filterdefs(curr_filterdefs_filterPos+2)));
                                if conversion_state_success_potaddval
                                      useSpecifiedFilterOrder = true;
                                      filterOrder = round(pot_additional_value);
                                      temp_additional_values = temp_additional_values + 1;
                                       if (numel(curr_filterdefs) >= curr_filterdefs_filterPos+3)
                                          pot_filterType = char(curr_filterdefs(curr_filterdefs_filterPos+3));
                                          if strcmp(pot_filterType,'iir') || strcmp(pot_filterType,'fir')
                                                overwriteFilterType = true;
                                                used_specified_filtering_once_lp = true;
                                                filterType = pot_filterType;
                                                temp_additional_values = temp_additional_values + 1;
                                          end
                                      end
                                end
                            end
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 2 + temp_additional_values;
                        case 'no'
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                        case {'+' '-'}
                            curr_operator_wait_for_second_term{end+1} = curr_filter;
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                            continue;
                        case '('
                            curr_bracket_opened = curr_bracket_opened + 1;
                            if numel(curr_operator_wait_for_second_term) > 0
                                data_filt{iChanCount} = data_store{end};
                                data_store(end) = [];
                            else
                                if numel(data_store) < 1
                                    data_store{end+1} = data_filt{iChanCount};%stack
                                else
                                    data_store{end+1} = data_store{end};%stack
                                end
                            end
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                            continue;
                        case ')'
                            if numel(curr_operator_wait_for_second_term) > 0
                                switch curr_operator_wait_for_second_term{end}
                                    case '+'
                                        data_filt{iChanCount}.trial{1} = data_term1.trial{1} + data_filt{iChanCount}.trial{1};
                                        data_term1 = [];
                                    case '-'
                                        data_filt{iChanCount}.trial{1} = data_term1.trial{1} - data_filt{iChanCount}.trial{1};
                                        data_term1 = [];
                                    otherwise
                                        error(['filter not well defined for channel ' curr_channel_label])
                                end
                                curr_operator_wait_for_second_term(end) = [];
                            else
                                if (curr_bracket_opened > 0) && ~(numel(curr_operator_wait_for_second_term) > 0)
                                    data_term1 = data_filt{iChanCount};
                                end
                            end
                            curr_bracket_opened = curr_bracket_opened - 1;
                            curr_filterdefs_filterPos = curr_filterdefs_filterPos + 1;
                            continue;
                        otherwise
                            error(['filter not well defined for channel ' curr_channel_label])
                    end
                catch err
                    error(['filter not well defined for channel ' curr_channel_label])
                end
                
                %[numeral conversion_successfull] = str2num(curr_filter);
                %            if conversion_successfull
                %                data_filt{iChanCount}.trial{1} = numeral;
                %            else
                
                
                if ~conversion_state_success
                    error(['converion of filter or multipication parameter failed for filter settings of channel '  curr_channel_label])
                end
                
                cfg = [];
                cfg = core_cfg;
                if curr_filterdefs_filterPos > 1
                    cfg.dftfilter = 'no';
                end
                cfg.channel = curr_channel_label;
                
                switch curr_filter
                    case 'normtomean'
                        data_filt{iChanCount}.trial{1} = (data_filt{iChanCount}.trial{1} - nanmean(data_filt{iChanCount}.trial{1}));
                    case 'min'
                        data_filt{iChanCount}.trial{1} = repmat(min(data_filt{iChanCount}.trial{1}),size(data_filt{iChanCount}.trial{1}));
                    case 'max'
                        data_filt{iChanCount}.trial{1} = repmat(max(data_filt{iChanCount}.trial{1}),size(data_filt{iChanCount}.trial{1}));
                    case 'ztransform'
                        data_filt{iChanCount}.trial{1} = (data_filt{iChanCount}.trial{1} - nanmean(data_filt{iChanCount}.trial{1}))./nanstd(data_filt{iChanCount}.trial{1});
                    case 'smooth'
                        data_filt{iChanCount}.trial{1} = smooth(data_filt{iChanCount}.trial{1},max(1,round(SmoothSeconds*data_filt{iChanCount}.fsample)),'moving')';
                    case 'envpeaks'
                        data_filt{iChanCount}.trial{1} = envpeaks(data_filt{iChanCount}.trial{1},data_filt{iChanCount}.fsample,minPeakDistanceSeconds);
                    case 'rect'
                        data_filt{iChanCount}.trial{1} = abs(data_filt{iChanCount}.trial{1});
                    case 'env'
                        data_filt{iChanCount}.trial{1} = abs(hilbert(data_filt{iChanCount}.trial{1}));
                    case 'mult'
                        data_filt{iChanCount}.trial{1} = data_filt{iChanCount}.trial{1} .* MultFactor;
                    case 'waveletband'
                        usedFilterOrder_bp = NaN;
                        bp_hdm = NaN;
                        
                        cfg.method = 'wavelet';
                        cfg.output = 'pow';
                        FreqSteps = abs(FpassRight-FpassLeft)/10;
                        cfg.foi = [FpassLeft:FreqSteps:FpassRight];
                        if useSpecifiedFilterOrder
                            cfg.width = filterOrder;
                        else
                            cfg.width =  4;
                        end
                        %cfg.pad = 'maxperlen';
                        cfg.toi = data_filt{iChanCount}.time{1};
                        cfg.pad = ( 2*cfg.width*(1/min(cfg.foi)) ) + ( max(cfg.toi)-min(cfg.toi) ) +  10 ;
                        cfg.feedback = core_cfg.feedback;
                        fprintf('dataset %i: reprocess and apply band wavlet filter to %s\n',iData,curr_channel_label{:});
                        data_freq = ft_freqanalysis(cfg,data_filt{iChanCount});
                        temp_freqdat = squeeze(nanmean(data_freq.powspctrm,2))';
                        temp_freqdat(isnan(temp_freqdat)) = 0;
                        data_filt{iChanCount}.trial = {temp_freqdat};
                        temp_freqdat = [];
                        dat_freq = [];
                    case 'bp'
                        cfg.bpfilter = 'yes';
                        FstopLeft = FpassLeft - StopToPassTransitionWidth_bp; %left stop frequency in Hz
                        FstopRight = FpassRight + PassToStopTransitionWidth_bp; %left stop frequency in Hz
                        
                        usedFilterOrder_bp = NaN;
                        bp_hdm = NaN;
                        if strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned') && ~overwriteFilterType
                            bp_d = [];
                            bp_hd = [];
                            fprintf('dataset %i: designing band pass filter\n',iData);
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
                        else
                            cfg.bpinstabilityfix = 'split';
                            switch filterType
                                case 'iir'
                                    cfg.bpfilttype =  'but';   
                                case 'fir'
                                    cfg.bpfilttype =  'fir';
                                otherwise
                                    error(['filter type ' filterType ' unknown'])
                            end
                        end
                        
                        if strcmp(UseFixedFilterOrder_bp,'yes') || useSpecifiedFilterOrder
                            if useSpecifiedFilterOrder 
                              cfg.bpfiltord = filterOrder;  
                            else
                                cfg.bpfiltord = FilterOrder_bp;
                            end
                            usedFilterOrder_bp = cfg.bpfiltord;
                        end
                        cfg.bpfreq        = [FpassLeft FpassRight];%dummy values are overwritten by low level function
                        cfg.feedback = core_cfg.feedback;
                        fprintf('dataset %i: reprocess and apply band filter to %s\n',iData,curr_channel_label{:});
                        data_filt{iChanCount} = ft_fw_preprocessing(cfg,data_filt{iChanCount});
                        
                    case 'hp'
                        cfg.hpfilter = 'yes';
                        FstopLeft = FpassLeft - StopToPassTransitionWidth_hp; %left stop frequency in Hz
                        
                        usedFilterOrder_hp = NaN;
                        hp_hdm = NaN;
                        if strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned') && ~overwriteFilterType
                            hp_d = [];
                            hp_hd = [];
                            if strcmp(UseFixedFilterOrder_hp,'yes')
                                hp_d = fdesign.highpass('N,F3db',FilterOrder_hp,FpassLeft,FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
                            else
                                hp_d = fdesign.highpass('Fst,Fp,Ast,Ap',FstopLeft,FpassLeft,AstopLeft_hp,Apass_hp,FrqOfSmpl);%designmethods(hp_preDS_d); help(hp_preDS_d,'equiripple'); help(hp_preDS_d,'butter')
                            end
                            fprintf('dataset %i: designing high pass filter\n',iData);
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
                        else
                            cfg.hpinstabilityfix = 'split';
                            switch filterType
                                case 'iir'
                                    cfg.hpfilttype =  'but';
                                case 'fir'
                                    cfg.hpfilttype =  'fir';
                                otherwise
                                    error(['filter type ' filterType ' unknown'])
                            end
                        end
                        if strcmp(UseFixedFilterOrder_hp,'yes') || useSpecifiedFilterOrder
                            if useSpecifiedFilterOrder 
                              cfg.hpfiltord = filterOrder;  
                            else
                                cfg.hpfiltord = FilterOrder_hp;
                            end
                            usedFilterOrder_hp = cfg.hpfiltord;
                        end
                        cfg.hpfreq        = [FpassLeft];%dummy values are overwritten by low level function
                        cfg.feedback = core_cfg.feedback;
                        fprintf('dataset %i: reprocess and apply high pass filter to %s\n',iData,curr_channel_label{:});
                        data_filt{iChanCount} = ft_fw_preprocessing(cfg,data_filt{iChanCount});
                        
                    case 'lp'
                        cfg.lpfilter = 'yes';
                        FstopRight = FpassRight + PassToStopTransitionWidth_lp; %right stop frequency in Hz
                        usedFilterOrder_lp = NaN;
                        lp_hdm = NaN;
                        if strcmp(core_cfg.lpfilttype,'IIRdesigned') || strcmp(core_cfg.lpfilttype,'FIRdesigned') && ~overwriteFilterType
                            lp_d = [];
                            lp_hd = [];
                            fprintf('dataset %i: designing low pass filter \n',iData);
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
                        else
                            cfg.lpinstabilityfix = 'split';
                            switch filterType
                                case 'iir'
                                    cfg.lpfilttype =  'but';
                                case 'fir'
                                    cfg.lpfilttype =  'fir';
                                otherwise
                                    error(['filter type ' filterType ' unknown'])
                            end
                        end
                        if strcmp(UseFixedFilterOrder_lp,'yes') || useSpecifiedFilterOrder
                            if useSpecifiedFilterOrder
                                cfg.lpfiltord = filterOrder;
                            else
                                cfg.lpfiltord = FilterOrder_lp;
                            end
                            usedFilterOrder_lp = cfg.lpfiltord;
                        end
                        cfg.lpfreq        = [FpassRight];%dummy values are overwritten by low level function
                        cfg.feedback = core_cfg.feedback;
                        fprintf('dataset %i: reprocess and apply low pass filter to %s\n',iData,curr_channel_label{:});
                        data_filt{iChanCount} = ft_fw_preprocessing(cfg,data_filt{iChanCount});
                    case 'no'
                        cfg.feedback = core_cfg.feedback;
                        fprintf('dataset %i: reprocess without filtering of %s\n',iData,curr_channel_label{:});
                        data_filt{iChanCount} = ft_fw_preprocessing(cfg,data_filt{iChanCount});
                    otherwise
                        error('filter in the filter definitions are not well defined, please check them')
                end
                
            end
            iChanCount = iChanCount + 1;
        end
        
        if length(data_filt) > 1
            data = ft_appenddata([],data_filt{:});
        else
            data = data_filt{:};
        end
        
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
    
    FrqOfSmpl = data.fsample;%data.hdr.Fs;%samples per second / Hz
    
    if strcmp(DoWriteData,'yes')
        
        data.trial = cat(2,data.trial(:));
        data.trial = {cat(2,data.trial{:})};
        data.time = cat(2,data.time(:));
        data.time = {cat(2,data.time{:})};
        
        hdr.nSamples = size(data.trial{1},2);%EEG.Points;
        hdr.nSamplesPre = 0;
        hdr.nTrials = 1;
        % ! potential is assumed to be in microV (µV) not milliVolt(mV)
        
        
        
        hdr.Fs = data.fsample;%EEG.SamplingRate;
        hdr.nChans = length(data.label);%EEG.ChannelNumber;
        hdr.nTrials = 1;
        if strcmp(IgnoreDataSetHeader,'no')
            chantype = {};
            chanunit = {};
            [dummy idat ihdr] = intersect(data.label,hdr.label);
            chantype(ismember(data.label,hdr.label)) = hdr.chantype(ihdr);
            chantype(~ismember(data.label,hdr.label)) = {'eeg'};
            chanunit(ismember(data.label,hdr.label)) = hdr.chanunit(ihdr);
            chanunit(~ismember(data.label,hdr.label)) = {DefaultOutputUnit};%EEG.ChannelUnits;%Nx1 cell-array with the physical units, see FT_CHANUNIT
            hdr.chantype = chantype';
            hdr.chanunit = chanunit';
        else
            hdr.chantype = repmat({'eeg'},hdr.nChans,1);%Nx1 cell-array with the channel type, see FT_CHANTYPE
            hdr.chanunit = repmat({DefaultOutputUnit},hdr.nChans,1);%EEG.ChannelUnits;%Nx1 cell-array with the physical units, see FT_CHANUNIT
        end
        hdr.label = data.label;%EEG.ChannelTitles;
        
        
        if strcmp(IncludePostiveMarkerAtBeginning,'yes')
            sigpositive_data = data.trial{1};
            sigpositive_data(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(sigpositive_data,1),1);
            data.trial{1} = sigpositive_data;
            sigpositive_data = [];
        end
        
        if strcmp(CutDataAtEndHypnogram,'yes') && ~strcmp(DoEpochData,'yes')
            epochLengthSamples = epochLength * data.fsample;
            [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
            if size(data.trial{1},2) > hypnEpochsEndsSamples(end)
                temp_data = data.trial{1};
                temp_data(:,(hypnEpochsEndsSamples(end)+1):end) = [];
                data.trial{1} = temp_data;
                temp_data = [];
                
                temp_time = data.time{1};
                temp_time((hypnEpochsEndsSamples(end)+1):end) = [];
                data.time{1} = temp_time;
                temp_time = [];
                
                %data.hdr = [];
                data.sampleinfo = [1 (hypnEpochsEndsSamples(end))];
            end
           
        end
        
        tempOutputDataformat = OutputDataformat;
        switch OutputDataformat
            case 'brainvision_eeg_int16'
                tempOutputDataformat = 'brainvision_eeg';
                hdr.brainvision_outformat = 'int16';%float32 int16 int32;
            case 'brainvision_eeg_int32'
                tempOutputDataformat = 'brainvision_eeg';
                hdr.brainvision_outformat = 'int32';%float32 int16 int32;
            case 'brainvision_eeg_float32'
                tempOutputDataformat = 'brainvision_eeg';
                hdr.brainvision_outformat = 'float32';%float32 int16 int32;
            case 'edf_autoscale'
                tempOutputDataformat = 'edf';
                hdr.edf_doautoscale = true;
            case 'edf_0.1uV_Ycuttoff'
                tempOutputDataformat = 'edf';
                hdr.edf_doautoscale = false;
                hdr.edf_accuracy = 0.1;
                hdr.edf_docutoff = true;
            case 'edf_0.01uV_Ycuttoff'
                tempOutputDataformat = 'edf';
                hdr.edf_doautoscale = false;
                hdr.edf_accuracy = 0.01;
                hdr.edf_docutoff = true;
            case 'edf_1uV_Ycuttoff'
                tempOutputDataformat = 'edf';
                hdr.edf_doautoscale = false;
                hdr.edf_accuracy = 1;
                hdr.edf_docutoff = true;
        end
        
        data_file_name = [pathOutputFolder filesep ouputFilesPrefixString 'preprocout_' 'datanum_' num2str(iData)];
        ft_write_data([data_file_name], data.trial{:},'dataformat',tempOutputDataformat,'header',hdr);
        
        
    end
    
    
    
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
        hp_preDS_hdm.Fs = preDownsampleFreq
        hp_preDS_hdm.Astop = NaN;
        hp_preDS_hdm.Fstop = NaN;
        hp_preDS_hdm.F6dB = NaN;
        hp_preDS_hdm.F3dB = PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff;
        hp_preDS_hdm.TransitionWidth = NaN;
        hp_preDS_hdm.Fpass = NaN;
        hp_preDS_hdm.Apass = NaN;
        
        usedFilterOrder_hp = NaN;
        hp_hdm.Fs = FrqOfSmpl
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
    
    
    if (usedFilterOrder_hp  == -1)%~exist('usedFilterOrder_hp','var')
        usedFilterOrder_hp = NaN;
        hp_hdm.Fs = FrqOfSmpl;
        hp_hdm.Astop = NaN;
        hp_hdm.Fstop = NaN;
        hp_hdm.F6dB = NaN;
        hp_hdm.F3dB = NaN;
        hp_hdm.TransitionWidth = NaN;
        hp_hdm.Fpass = NaN;
        hp_hdm.Apass = NaN;
    end
    
    if (usedFilterOrder_bp  == -1)%~exist('usedFilterOrder_bp','var')
        usedFilterOrder_bp = NaN;
        bp_hdm.Fs = FrqOfSmpl;
        bp_hdm.Astop1 = NaN;
        bp_hdm.TransitionWidth1 = NaN;
        bp_hdm.F3dB1 = NaN;
        bp_hdm.F6dB1 = NaN;
        bp_hdm.Fpass1 = NaN;
        bp_hdm.Apass = NaN;
        bp_hdm.Fpass2 = NaN;
        bp_hdm.F3dB2 = NaN;
        bp_hdm.F6dB2 = NaN;
        bp_hdm.TransitionWidth2 = NaN;
        bp_hdm.Astop2 = NaN;
        bp_hdm.Fstop1 = NaN;
        bp_hdm.Fstop2 = NaN;
    end
    
    if (usedFilterOrder_lp == -1)%~exist('usedFilterOrder_lp','var')
        usedFilterOrder_lp = NaN;
        lp_hdm.Fs = FrqOfSmpl;
        lp_hdm.Astop = NaN;
        lp_hdm.Fstop = NaN;
        lp_hdm.F6dB = NaN;
        lp_hdm.F3dB = NaN;
        lp_hdm.TransitionWidth = NaN;
        lp_hdm.Fpass = NaN;
        lp_hdm.Apass = NaN;
    end
    
    
    if false
    
    fidf = fopen([pathOutputFolder filesep ouputFilesPrefixString 'preprocout_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
    %write header
    fprintf(fidf,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'...
        ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'...
        ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'...
        ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'datasetnum','dataset',...
        'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB',...
        'hp_filter','hp_filter_type','hp_dir_and_passing','usedFilterOrder_hp','hp_Fs_Hz','hp_Astop_dB','hp_Fstop_Hz','hp_F6dB_Hz','hp_F3dB_Hz','hp_TransitionWidth_Hz','hp_Fpass_Hz','hp_Apass_dB',...
        'bp_filter','bp_filter_type','bp_dir_and_passing','usedFilterOrder_bp','bp_Fs_Hz','bp_Astop1_dB','bp_Fstop1_Hz','bp_TransitionWidth1_Hz','bp_F3dB1_Hz','bp_F6dB1_Hz','bp_Fpass1_Hz','bp_Apass_dB','bp_Fpass2_Hz','bp_F3dB2_Hz','bp_F6dB2_Hz','bp_TransitionWidth2_Hz','bp_Fstop2_Hz','bp_Astop2_dB',...
        'lp_filter','lp_filter_type','lp_dir_and_passing','usedFilterOrder_lp','lp_Fs_Hz','lp_Astop_dB','lp_Fstop_Hz','lp_F6dB_Hz','lp_F3dB_Hz','lp_TransitionWidth_Hz','lp_Fpass_Hz','lp_Apass_dB');
    
    
    
    %write content
    fprintf(fidf,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e'...
        ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e'...
        ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e,%f,%f,%f,%f,%f,%e'...
        ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e'  '\n'],...
        iData,datasetsPath,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp,hp_hdm.Fs,hp_hdm.Astop,hp_hdm.Fstop,hp_hdm.F6dB,hp_hdm.F3dB,hp_hdm.TransitionWidth,hp_hdm.Fpass,hp_hdm.Apass,...
        core_cfg.bpfilttype,bp_f_type_detail,core_cfg.bpfiltdir,usedFilterOrder_bp,bp_hdm.Fs,bp_hdm.Astop1,bp_hdm.Fstop1,bp_hdm.TransitionWidth1,bp_hdm.F3dB1,bp_hdm.F6dB1,bp_hdm.Fpass1,bp_hdm.Apass,bp_hdm.Fpass2,bp_hdm.F3dB2,bp_hdm.F6dB2,bp_hdm.TransitionWidth2,bp_hdm.Fstop2,bp_hdm.Astop2,...
        core_cfg.lpfilttype,lp_f_type_detail,core_cfg.lpfiltdir,usedFilterOrder_lp,lp_hdm.Fs,lp_hdm.Astop,lp_hdm.Fstop,lp_hdm.F6dB,lp_hdm.F3dB,lp_hdm.TransitionWidth,lp_hdm.Fpass,lp_hdm.Apass);
    
    
    fclose(fidf);
    
    end
    
    %     cfg_datbrow = [];
    %     cfg_datbrow.viewmode  = 'vertical';%'butterfly', 'vertical', 'component' for visualizing components e.g. from an ICA (default is 'butterfly')
    %     cfg_datbrow.blocksize = epochLength;
    %     cfg_datbrow.artfctdef.markartifact.artifact = [300 400; 600 900];
    %     cfg_datbrow.selectfeature = {'markartifact'};
    %     cfg_datbrow.selectmode = 'markartifact';
    %     cfg_datbrowres = ft_databrowser(cfg_datbrow, data);
    %     cfg_datbrowres = ft_databrowser(cfg_datbrowres, data);
    
    
    data = [];%clear
    
end

if false
%aggregate all results from datasets
fprintf('Aggregate results of all datasets\n');
fidf_all = [];
delimiter = ',';
for iData = iDatas
    
    fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'preprocout_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    
    
    if iData == iDatas(1)
        fidf_all = fidf;
        
    else
        
        
        fidf_all = cat(1,fidf_all,fidf);
        
    end
    
end
export(fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'preprocout_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);


res_filters = fidf_all;
end; res_filters = [];

fprintf('PREPROCOUT function finished\n');
toc
memtoc
end

function signal = envpeaks(x,fsample,seconds)  
    [pks locs] = findpeaks(x,'MINPEAKDISTANCE',max(1,floor(seconds*fsample)));
    signal = interp1([1 locs numel(x)],[x(1) pks x(end)],1:numel(x),'linear');   
end

