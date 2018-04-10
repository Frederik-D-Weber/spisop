function [res_hypnograms, res_filters] = spisop_confuser_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% generate simulated labes for epoched data.
% Copyright Frederik D. Weber

functionName = 'confuser';
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

DoProcessData = getParam('DoProcessData',listOfParameters);%DoProcessData = 'yes'

HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);

DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
ChannelsOfInterestFileName = getParam('ChannelsOfInterestFileName',listOfParameters);
AVGoverChannels = getParam('AVGoverChannels',listOfParameters);

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

if ~strcmp(DoProcessData,'no')

if exist([pathInputFolder filesep ChannelsOfInterestFileName],'file') ~= 2
    error(['ChannelsOfInterestFileName file ' [pathInputFolder filesep ChannelsOfInterestFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end

end

PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff = str2num(getParam('PreDownSampleHighPassFilter_FpassLeft_or_F3dBcutoff',listOfParameters));%in Hz




epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds
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

if ~strcmp(DoProcessData,'no')
listOfChannelsOfInterest = read_mixed_csv([pathInputFolder filesep ChannelsOfInterestFileName],',');
end

listOfHypnogramPaths = read_mixed_csv([pathInputFolder filesep HypnogramsFileName],',');

if ~strcmp(DoProcessData,'no')

    
if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfChannelsOfInterest,1)))
    error('files or number of Datasetspaths Hypnogramsfiles ChannelsOfInterest are invalid or do not aggree')
end
else
if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)))
    error('files or number of Datasetspaths Hypnogramsfiles are invalid or do not aggree')
end    
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end



IncludeMAepochs = getParam('IncludeMAepochs',listOfParameters); % either yes or no default yes.

includeMAepochsFlag = strcmp(IncludeMAepochs,'yes');

RandomSeed = str2num(getParam('RandomSeed',listOfParameters));
rng(RandomSeed, 'twister');


MaxChangePercentage = str2num(getParam('MaxChangePercentage',listOfParameters));
if (MaxChangePercentage > 100) || (MaxChangePercentage < 0)
    error(['MaxChangePercentage must be between 0 and 100 %.'])
end

MaxSwapPercentage = MaxChangePercentage/2;


if ~strcmp(DoProcessData,'no')
    pretestHeaderForPersistentSampleFrequencies(IgnoreDataSetHeader,iDatas,listOfDatasetHeaderPaths,listOfChannelsOfInterest,FrqOfSmplWished);
end

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





tic
memtic
fprintf('confuser function initialized\n');
conseciDatas = 1:length(iDatas);
for conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    FrqOfSmplWishedPar = FrqOfSmplWished;
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    
if ~strcmp(DoProcessData,'no')
    

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
    %[roiBegins, roiEnds] = getROIsByHypnogramMAoption(hypnogramPath,epochLengthSamples,sleepStagesOfInterest,includeMAepochsFlag);

    [hypn_temp hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);

    roiBegins = hypnEpochsBeginsSamples;
    roiEnds = hypnEpochsEndsSamples;
    
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
    %[roiBegins, roiEnds] = getROIsByHypnogramMAoption(hypnogramPath,epochLengthSamples,sleepStagesOfInterest,includeMAepochsFlag);

    [hypn_temp hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);

    roiBegins = hypnEpochsBeginsSamples;
    roiEnds = hypnEpochsEndsSamples;
    

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
       
%     cfg = [];
%     cfg.length    = NewEpochLength;%single number (in unit of time, typically seconds) of the required snippets
%     cfg.feedback = core_cfg.feedback;
%     data = ft_redefinetrial(cfg,data);
%     NnewSegments = size(data.trial,2);
    
else 
    epochLengthSamples = epochLength;
end
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    tempRepeatNumber = floor(epochLength/epochLength);
%     hypnNew = [];
%     hypnStagesNew = [];
%     for iHyp = 1:size(hypn,1)
%         hypnNew = [hypnNew;repmat(hypn(iHyp,:),tempRepeatNumber,1)];
%         hypnStagesNew = [hypnStagesNew;repmat(hypn(iHyp,:),tempRepeatNumber,1)];
%     end
    
   hypnNew = hypn;
   
if ~strcmp(DoProcessData,'no')
  
    
    hypnEpochsNew = (1:size(hypnNew,1));
    
    newEpochLengthSamples = epochLength * FrqOfSmpl;
    

    hypnEpochsBeginsSamplesNew = (((hypnEpochsNew - 1) * newEpochLengthSamples) + 1)';
    hypnEpochsEndsSamplesNew = (hypnEpochsNew * newEpochLengthSamples)';

end
%     hypnNew
%     hypnEpochsNew
%     hypnStagesNew

    nStages = size(hypnNew,1);
    
    nStages_change = max(0,floor(nStages*MaxSwapPercentage/100));
    
    stage_types = unique(hypnNew(:,1));
    
    switch_pairs_ind = randperm(nStages);
    %switch_pairs_ind_to = switch_pairs_ind(2:end);
    nSwitches = 0;
    for iSt_from = 1:nStages_change
        from_switch_ind = switch_pairs_ind(iSt_from);
        switch_alternatives = switch_pairs_ind;
        switch_alternatives(iSt_from) = [];
        switch_alternatives = switch_alternatives(randperm(numel(switch_alternatives)));%shuffle
        numberSwitchpairsLeft = numel(switch_alternatives);
        for iSt_to = 1:numberSwitchpairsLeft
            to_switch_ind = switch_alternatives(iSt_to);
            if (hypnNew(to_switch_ind,1) ~= hypnNew(from_switch_ind,1)) || (iSt_to == numberSwitchpairsLeft)
                hypnNew(to_switch_ind,:) = hypn(from_switch_ind,:);
                hypnNew(from_switch_ind,:) = hypn(to_switch_ind,:); 
                %switch_pairs_ind_to(iSt_to) = [];
                if (iSt_to ~= numberSwitchpairsLeft)
                    nSwitches = nSwitches + 1;
                end
                break
            end
        end
    end
    fprintf('dataset %i: switched %i and changed %i of %i hypnogram stages (%f percent)\n',iData,nSwitches,nSwitches*2,nStages,100*nSwitches*2/nStages);

    
   [numberEEG numberEOG numberEMG] = getScoringChannelNumbers(data.label);
    
    cfg = [];
    cfg.feedback = core_cfg.feedback;
    cfg.channel = ft_channelselection(data.label{numberEEG}, data.label);
    chData_EEG = ft_selectdata(cfg,data);
    
    cfg = [];
    cfg.feedback = core_cfg.feedback;
    cfg.channel = ft_channelselection(data.label{numberEOG}, data.label);
    chData_EOG = ft_selectdata(cfg,data);
    
    cfg = [];
    cfg.feedback = core_cfg.feedback;
    cfg.channel = ft_channelselection(data.label{numberEMG}, data.label);
    chData_EMG = ft_selectdata(cfg,data);
    
    cfg = [];
    %cfg.feedback = core_cfg.feedback;
    cfg.channel = ft_channelselection(data.label{numberEEG}, data.label);
    cfg.method = 'mtmfft';
    cfg.output = 'pow';
    cfg.keeptrials = 'yes';
    cfg.foilim = [0.5 20];
    cfg.tapsmofrq   = 1;
    chData_EEG = ft_freqanalysis(cfg,data);
   
   
    nChannels = length(data.label);
    %W = (trlSampleLengths./sampleLengthsAcrossROIs); %Nx1
    nextHypColumn = 3+iChan;
    hypnNew(:,nextHypColumn) = 0;
    
    
    fprintf('dataset %i: process channel %s\n',iData,data.label{iChan});
    for iTr = 1:(size(chData_EEG.trial,2))
        
        rawDataSampleOffset = trlSampleBeginsAndEnds(iTr,1) - 1;
        
        tempSignalEEG = chData_EEG.trial{iTr};
        tempSignalEOG = chData_EOG.trial{iTr};
        tempSignalEMG = chData_EMG.trial{iTr};
        
        
        
        
        %              minTempSignal = min(tempSignal);
        %              maxTempSignal = max(tempSignal);
        %
        %              tempSignalMinSample = rawDataSampleOffset + find(minTempSignal == tempSignal);
        %              tempSignalMaxSample = rawDataSampleOffset + find(maxTempSignal == tempSignal);
        %
        %              leftSample = hypnEpochsBeginsSamplesNew(iTr);
        %              rightSample = hypnEpochsEndsSamplesNew(iTr);
        %
        %              if strcmp(MinMaxAllowedPotentialCriterion,'yes')
        %                  tempMinSignalCriterionHits = (minTempSignal < MinAllowedPotential);
        %                  tempMaxSignalCriterionHits = (maxTempSignal > MaxAllowedPotential);
        %                  hypnNew(iTr,nextHypColumn) = hypnNew(iTr,nextHypColumn) + (tempMinSignalCriterionHits || tempMaxSignalCriterionHits);
        %
        %                  leftSample = max(leftSample, min(tempSignalMinSample, tempSignalMaxSample) + hypnEpochsBeginsSamplesNew(iTr) - 1);
        %                  rightSample = min(rightSample, max(tempSignalMinSample, tempSignalMaxSample) + hypnEpochsBeginsSamplesNew(iTr) - 1);
        %
        %              end
        %
        %              leftPaddingReach = leftSample - PreArtifactPaddingSamples;
        %              rightPaddingReach = rightSample + PostArtifactPaddingSamples;
        %              leftPaddingMarkingsIndex = find(hypnNew(hypnEpochsBeginsSamplesNew <= leftPaddingReach),1,'first'):iTr;
        %              rightPaddingMarkingsIndex = iTr:find(hypnNew(hypnEpochsEndsSamplesNew >= rightPaddingReach),1,'first');
        %              hypnNew(leftPaddingMarkingsIndex,nextHypColumn) = hypnNew(leftPaddingMarkingsIndex,nextHypColumn) + 1;
        %              hypnNew(rightPaddingMarkingsIndex,nextHypColumn) = hypnNew(rightPaddingMarkingsIndex,nextHypColumn) + 1;
        %              hypnNew(iTr,nextHypColumn) = hypnNew(iTr,nextHypColumn) - 2; %subtract added too much
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
    
if ~strcmp(DoProcessData,'no')

    fidff = fopen([pathOutputFolder filesep ouputFilesPrefixString 'confuser_filter_' 'datanum_' num2str(iData) '.csv'],'wt');
    %write header
    fprintf(fidff,['%s,%s' ',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'datasetnum','dataset',...
        'hp_preDS_filter','hp_preDS_filter_type','hp_dir_and_passing','usedFilterOrder_hp_preDS','hp_preDS_Fs_Hz','hp_preDS_Astop_dB','hp_preDS_Fstop_Hz','hp_preDS_F6dB_Hz','hp_preDS_F3dB_Hz','hp_preDS_TransitionWidth_Hz','hp_preDS_Fpass_Hz','hp_preDS_Apass_dB');
    %write content
    fprintf(fidff,['%i,%s' ',%s,%s,%s,%i,%i,%e,%f,%f,%f,%f,%f,%e' '\n'],...
        iData,datasetsPath,...
        core_cfg.hpfilttype,hp_f_type_detail,core_cfg.hpfiltdir,usedFilterOrder_hp_preDS,hp_preDS_hdm.Fs,hp_preDS_hdm.Astop,hp_preDS_hdm.Fstop,hp_preDS_hdm.F6dB,hp_preDS_hdm.F3dB,hp_preDS_hdm.TransitionWidth,hp_preDS_hdm.Fpass,hp_preDS_hdm.Apass);
    
    
    fclose(fidff);
end
% 
%     %open output files
%     fidc = fopen([pathOutputFolder filesep ouputFilesPrefixString 'confuser_hypnograminfo_' 'datanum_' num2str(iData) '.csv'],'wt');
%     
%     %write header of ouptufiles
%     fprintf(fidc,'%s,%s,%s,%s,%s\n',...
%         'datasetnum','epoch_length_seconds');
%     
%     fprintf(fidc,'%i,',iData);
%     fprintf(fidc,'%f\n',epochLength);
% 
%     %fprintf(fidc,'%f,',sum(hypnNew(:,3))*epochLength);
%     %fprintf(fidc,'%f\n',lengthsAcrossROIsSeconds);
%     
%     fclose(fidc);

    T = table(hypnNew);
    writetable(T,[pathOutputFolder filesep ouputFilesPrefixString 'confuser_hypn_' 'datanum_' num2str(iData) '.txt'],'Delimiter','\t','WriteVariableNames',false);
   
%     T = table(hypnNew(:,[1 3]));
%     writetable(T,[pathOutputFolder filesep ouputFilesPrefixString 'confounded_hypn_new_' 'datanum_' num2str(iData) '.txt'],'Delimiter','\t','WriteVariableNames',false);

%     temp_ind_art = (hypnNew(:,2) > 0) | (hypnNew(:,3) > 0);
%     hypnNew_consens = hypnNew(:,2);
%     hypnNew_consens(temp_ind_art) = max(hypnNew(temp_ind_art,2:3),[],2);
%     T = table([hypnNew(:,1),hypnNew_consens]);
%     writetable(T,[pathOutputFolder filesep ouputFilesPrefixString 'confounded_hypn_new_conserv_' 'datanum_' num2str(iData) '.txt'],'Delimiter','\t','WriteVariableNames',false);


end

fprintf('Aggregate results of all datasets\n');
if ~strcmp(DoProcessData,'no')

%aggregate all results from datasets
fidff_all = [];
%fidc_all = [];
delimiter = ',';
for iData = iDatas
    
    fidff = dataset('File', [pathOutputFolder filesep ouputFilesPrefixString 'confuser_filter_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    %fidc = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'confounds_hypnograminfo_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    
    if iData == iDatas(1)
        fidff_all = fidff;
        %fidc_all = fidc;
    else
        fidff_all = cat(1,fidff_all,fidff);
       % fidc_all = cat(1,fidc_all,fidc);
    end
    
end
export(fidff_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'confuser_filter_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
%export(fidc_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'confuser_hypnograminfo_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);

res_filters = fidff_all;
%res_confoundsinfo = fidc_all;

end

fprintf('confuser function finished\n');
toc
memtoc
end

