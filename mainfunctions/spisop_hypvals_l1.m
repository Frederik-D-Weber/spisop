function res_hypnvals = spisop_hypvals_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% determine sleep scoring values for sleep table, assumes there is at least
% one valid epoch of S2 or S3 or S4 or REM
% Copyright Frederik D. Weber

DataSetPathsFileName = getParam('DataSetPathsFileName',listOfCoreParameters);
DataSetHeaderPathsFileName = getParam('DataSetHeaderPathsFileName',listOfCoreParameters);
IgnoreDataSetHeader = getParam('IgnoreDataSetHeader',listOfCoreParameters);
HypnogramsFileName = getParam('HypnogramsFileName',listOfCoreParameters);
LightsOffsFileName = getParam('LightsOffsFileName',listOfParameters);

LightsOffsMomentUnit = 'sample';
try
LightsOffsMomentUnit = getParam('LightsOffsMomentUnit',listOfParameters);

if ~(strcmp(LightsOffsMomentUnit,'sample') || strcmp(LightsOffsMomentUnit,'second'))
        error(['LightsOffsMomentUnit parameter must either be sample or second but given was ' LightsOffsMomentUnit ''])
end
catch e
end

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
if exist([pathInputFolder filesep LightsOffsFileName],'file') ~= 2
    error(['LightsOffsFileName file ' [pathInputFolder filesep LightsOffsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end



epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds

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
listOfLightsOffs = load([pathInputFolder filesep LightsOffsFileName]);

SleepOnsetDefinition = getParam('SleepOnsetDefinition',listOfParameters);


ExportHypnogram = getParam('ExportHypnogram',listOfParameters);%either yes or no default no
ExportHypnogramStartOffset = getParam('ExportHypnogramStartOffset',listOfParameters);%either sleeponset or eegonset or lightsoff default sleeponset
ExportHypnogramStartTime = str2num(getParam('ExportHypnogramStartTime',listOfParameters));%time after ExportHypnogramStartOffset to cut in minutes default 0
ExportHypnogramEndTime = str2num(getParam('ExportHypnogramEndTime',listOfParameters));%time after ExportHypnogramStartOffset to cut in minutes default 90


if (strcmp(ExportHypnogram,'yes'))
    if ~(ExportHypnogramEndTime >= ExportHypnogramStartTime)
        error(['Parameter ExportHypnogramEndTime = ' num2str(ExportHypnogramEndTime) ' needs to be at least the value of ExportHypnogramStartTime = ' num2str(ExportHypnogramStartTime)]);
    end
end


if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfLightsOffs,1)))
    error('files or number of Datasetspaths and Hypnogramsfiles and LightsOffs are invalid or do not aggree')
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end


core_cfg = [];
core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);

tempExportPostfix = '_';

if strcmp(ExportHypnogram,'yes')
        if strcmp(ExportHypnogramStartOffset,'sleeponset')
            tempExportPostfix = [tempExportPostfix 'sleeponset'];
        elseif strcmp(ExportHypnogramStartOffset,'eegonset')
            tempExportPostfix = [tempExportPostfix 'eegonset'];
        elseif strcmp(ExportHypnogramStartOffset,'lightsoff')
            tempExportPostfix = [tempExportPostfix 'lightsoff'];            
        else
            error(['wrong parameter for ExportHypnogramStartOffset either sleeponset or eegonset but set is ' ExportHypnogramStartOffset] );
        end
end

tempExportPostfix = [tempExportPostfix '_' num2str(ExportHypnogramStartTime) '_to_' num2str(ExportHypnogramEndTime) '_min'];



sleepOnsetTime = cell(1,numel(iDatas));

S1OnsetTime = cell(1,numel(iDatas));
S2OnsetTime = cell(1,numel(iDatas));

SWSonsetTime = cell(1,numel(iDatas));
S4onsetTime = cell(1,numel(iDatas));
REMonsetTime = cell(1,numel(iDatas));


totalSleepTime = cell(1,numel(iDatas));

S1Time = cell(1,numel(iDatas));
S2Time = cell(1,numel(iDatas));
S3Time = cell(1,numel(iDatas));
S4Time = cell(1,numel(iDatas));
REMtime = cell(1,numel(iDatas));
WakeTime = cell(1,numel(iDatas));
MovementTime = cell(1,numel(iDatas));
SWStime = cell(1,numel(iDatas));
NonREMtime = cell(1,numel(iDatas));


S1TimePreOnset = cell(1,numel(iDatas));
S2TimePreOnset = cell(1,numel(iDatas));
S3TimePreOnset = cell(1,numel(iDatas));
S4TimePreOnset = cell(1,numel(iDatas));
REMtimePreOnset = cell(1,numel(iDatas));
WakeTimePreOnset = cell(1,numel(iDatas));
MovementTimePreOnset = cell(1,numel(iDatas));
SWStimePreOnset = cell(1,numel(iDatas));
NonREMtimePreOnset = cell(1,numel(iDatas));

S1Time_WithoutMA = cell(1,numel(iDatas));
S2Time_WithoutMA = cell(1,numel(iDatas));
S3Time_WithoutMA = cell(1,numel(iDatas));
S4Time_WithoutMA = cell(1,numel(iDatas));
REMtime_WithoutMA = cell(1,numel(iDatas));
WakeTime_WithoutMA = cell(1,numel(iDatas));
MovementTime_WithoutMA = cell(1,numel(iDatas));
SWStime_WithoutMA = cell(1,numel(iDatas));
NonREMtime_WithoutMA = cell(1,numel(iDatas));


totalSleepTime_export = cell(1,numel(iDatas));

S1Time_export = cell(1,numel(iDatas));
S2Time_export = cell(1,numel(iDatas));
S3Time_export = cell(1,numel(iDatas));
S4Time_export = cell(1,numel(iDatas));
REMtime_export = cell(1,numel(iDatas));
WakeTime_export = cell(1,numel(iDatas));
MovementTime_export = cell(1,numel(iDatas));
SWStime_export = cell(1,numel(iDatas));
NonREMtime_export = cell(1,numel(iDatas));


S1Time_WithoutMA_export = cell(1,numel(iDatas));
S2Time_WithoutMA_export = cell(1,numel(iDatas));
S3Time_WithoutMA_export = cell(1,numel(iDatas));
S4Time_WithoutMA_export = cell(1,numel(iDatas));
REMtime_WithoutMA_export = cell(1,numel(iDatas));
WakeTime_WithoutMA_export = cell(1,numel(iDatas));
MovementTime_WithoutMA_export = cell(1,numel(iDatas));
SWStime_WithoutMA_export = cell(1,numel(iDatas));
NonREMtime_WithoutMA_export = cell(1,numel(iDatas));




tic
memtic
fprintf('HypVals function initialized\n');
parfor iData = iDatas
    %iData = 11
    
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    
    lightsOffMoment = listOfLightsOffs(iData);
   
    
    hdr = [];
    preDownsampleFreq = 0;
    if strcmp(IgnoreDataSetHeader,'no')
        headerPath = listOfDatasetHeaderPaths{iData};
        hdr = ft_read_header(headerPath);
        preDownsampleFreq = hdr.Fs;
    elseif strcmp(IgnoreDataSetHeader,'yes')
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
        tempdata = ft_preprocessing(cfg);
        preDownsampleFreq = tempdata.fsample;
        tempdata = [];
    else
        error('wrong parameter for IgnoreDataSetHeader either yes or no');
    end
    
    fprintf('dataset %i: process ROI from hypnogram info\n',iData);
    
    
    if strcmp(LightsOffsMomentUnit,'second')
        lightsOffSample = round(lightsOffMoment*preDownsampleFreq);
    else
        lightsOffSample = lightsOffMoment;
    end
    
    
    sampleFreq = preDownsampleFreq;
    epochLengthSamples = epochLength * preDownsampleFreq;
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    
    
    % Sleep onsset: [S1] S1 S2, but not [S1] S1 X S2, where X is not NonREM
    % S1, starting with first S1 is sleep onset, nothing else
    % S1_NonREM, starting with first S1 followed directly by either S2, S3 or S4,
    %            otherwise with first S2 or S3 or S4
    % S1_XREM, starting with first S1 directly followed by either S2, S3, S4, or REM,
    %            otherwise with first S2 or S3 or S4 or REM
    % NonREM, starting with first one of S2, S3 or S4
    % XREM, starting with first one of S2, S3 or S4 or REM
    
    
    if strcmp(SleepOnsetDefinition,'NonREM')
        onsetCandidate = -1;
        
        for iOnset = 1:(size(hypnStages,1))
            if strcmp(hypnStages(iOnset,3),'NonREM') && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
                onsetCandidate = iOnset;
                break;
            end
        end
        
    elseif strcmp(SleepOnsetDefinition,'XREM')
        onsetCandidate = -1;
        
        for iOnset = 1:(size(hypnStages,1))
            if (strcmp(hypnStages(iOnset,3),'NonREM') ||  strcmp(hypnStages(iOnset,3),'REM')) && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
                onsetCandidate = iOnset;
                break;
            end
        end
        
    elseif strcmp(SleepOnsetDefinition,'S1') || strcmp(SleepOnsetDefinition,'S1_NonREM') || strcmp(SleepOnsetDefinition,'S1_XREM')
        
        onsetCandidate = -1;
        consecS1 = 0;
        hasS1 = logical(0);
        for iOnset = 1:(size(hypnStages,1))
            if strcmp(hypnStages(iOnset,1),'S1') && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
                hasS1 = logical(1);
                consecS1 = consecS1 + 1;
                if ((onsetCandidate+consecS1) ~= iOnset)
                    onsetCandidate = iOnset;
                    consecS1 = 0;
                end
                if strcmp(SleepOnsetDefinition,'S1')
                    break;
                end
            elseif ( strcmp(SleepOnsetDefinition,'S1_XREM') && (strcmp(hypnStages(iOnset,3),'NonREM') || strcmp(hypnStages(iOnset,3),'REM')) ) ...
                    || ( strcmp(SleepOnsetDefinition,'S1_NonREM') && (strcmp(hypnStages(iOnset,3),'NonREM')) ) ...
                    && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
                if ~hasS1
                    onsetCandidate = iOnset;
                end
                break;
            else
                consecS1 = 0;
                hasS1 = logical(0);
            end
        end
        
    end
    
    
    
    sleepOnsetTime{iData} = (hypnEpochsBeginsSamples(onsetCandidate) - lightsOffSample)/sampleFreq;
    
    %preOffsetCandidate = max(find(strcmp(hypnStages(:,1),'S1') | strcmp(hypnStages(:,3),'NonREM')));
    preOffsetCandidate = max(find(strcmp(hypnStages(:,1),'S1') | strcmp(hypnStages(:,3),'NonREM') | strcmp(hypnStages(:,3),'REM') | strcmp(hypnStages(:,3),'MT')));
    
    S1ind = find(strcmp(hypnStages(:,1),'S1'));
    S2ind = find(strcmp(hypnStages(:,1),'S2'));
    SWSind = find(strcmp(hypnStages(:,2),'SWS'));
    S4ind = find(strcmp(hypnStages(:,1),'S4'));
    REMind = find(strcmp(hypnStages(:,1),'REM'));
    S1OnsetTime{iData} = (min(S1ind(S1ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    S2OnsetTime{iData} = (min(S2ind(S2ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    SWSonsetTime{iData} = (min(SWSind(SWSind >= onsetCandidate)) - onsetCandidate)*epochLength;
    S4onsetTime{iData} = (min(S4ind(S4ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    REMonsetTime{iData} = (min(REMind(REMind >= onsetCandidate)) - onsetCandidate)*epochLength;
    
    preOnsetCandidate = onsetCandidate;
    if preOnsetCandidate > 1
        preOnsetCandidate = preOnsetCandidate-1;
    end
    hypnTST = hypn(onsetCandidate:preOffsetCandidate,:);
    hypnStagesTST = hypnStages(onsetCandidate:preOffsetCandidate,:);
    hypnStagesPreSleepOnset = hypnStages(1:preOnsetCandidate,:);
    hypnEpochsTST = hypnEpochs(onsetCandidate:preOffsetCandidate);
    hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(onsetCandidate:preOffsetCandidate,:);
    hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(onsetCandidate:preOffsetCandidate,:);
    
    totalSleepTime{iData} = (length(onsetCandidate:preOffsetCandidate))*epochLength;
    
    S1Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S1')))*epochLength;
    S2Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S2')))*epochLength;
    S3Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S3')))*epochLength;
    S4Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S4')))*epochLength;
    REMtime{iData} = length(find(strcmp(hypnStagesTST(:,1),'REM')))*epochLength;
    WakeTime{iData} = length(find(strcmp(hypnStagesTST(:,1),'Wake')))*epochLength;
    MovementTime{iData} = length(find(strcmp(hypnStagesTST(:,1),'MT')))*epochLength;
    SWStime{iData} = S3Time{iData} + S4Time{iData};
    NonREMtime{iData} = SWStime{iData} + S2Time{iData};
    
    
    S1TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S1')))*epochLength;
    S2TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S2')))*epochLength;
    S3TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S3')))*epochLength;
    S4TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S4')))*epochLength;
    REMtimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'REM')))*epochLength;
    WakeTimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'Wake')))*epochLength;
    MovementTimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'MT')))*epochLength;
    SWStimePreOnset{iData} = S3TimePreOnset{iData} + S4TimePreOnset{iData};
    NonREMtimePreOnset{iData} = SWStimePreOnset{iData} + S2TimePreOnset{iData};
    
    S1Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S1') & (hypnTST(:,2) == 0) ))*epochLength;
    S2Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S2') & (hypnTST(:,2) == 0) ))*epochLength;
    S3Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S3') & (hypnTST(:,2) == 0) ))*epochLength;
    S4Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S4') & (hypnTST(:,2) == 0) ))*epochLength;
    REMtime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'REM') & (hypnTST(:,2) == 0) ))*epochLength;
    WakeTime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'Wake') & (hypnTST(:,2) == 0) ))*epochLength;
    MovementTime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'MT') & (hypnTST(:,2) == 0) ))*epochLength;
    SWStime_WithoutMA{iData} = S3Time_WithoutMA{iData} + S4Time_WithoutMA{iData};
    NonREMtime_WithoutMA{iData} = SWStime_WithoutMA{iData} + S2Time_WithoutMA{iData};
    
    
    if strcmp(ExportHypnogram,'yes')
        tempExportHypnogramStart = NaN;
        tempExpochFactor = 60/epochLength;
        if strcmp(ExportHypnogramStartOffset,'sleeponset')
            tempExportHypnogramStart = onsetCandidate;
        elseif strcmp(ExportHypnogramStartOffset,'eegonset')
            tempExportHypnogramStart = 1;
        elseif strcmp(ExportHypnogramStartOffset,'lightsoff')
            tempExportHypnogramStart = find(hypnEpochsEndsSamples <= lightsOffSample,1,'last')-1;
            
        else
            error(['wrong parameter for ExportHypnogramStartOffset either sleeponset or eegonset but set is ' ExportHypnogramStartOffset] );
        end
        
        
        
        tempExportHypnogramStartTimeEpoch = tempExportHypnogramStart + round(ExportHypnogramStartTime*tempExpochFactor);
        tempExportHypnogramEndTimeEpoch = tempExportHypnogramStart + round(ExportHypnogramEndTime*tempExpochFactor);
        tempExport_hypn = hypn;
        tempExport_hypnStages = hypnStages;
        
        if ((tempExportHypnogramStartTimeEpoch-1) > 0)
            tempExport_hypn(1:(tempExportHypnogramStartTimeEpoch-1),2) = 3;
        end
        
        if (tempExportHypnogramEndTimeEpoch <= size(tempExport_hypn,1))
            tempExport_hypn((tempExportHypnogramEndTimeEpoch):end,2) = 3;
        end
        
        
        [temp_pathstr,temp_name,temp_ext] = fileparts(hypnogramPath);
        fid_export = fopen([ouputFilesPrefixString num2str(iData) '_' temp_name tempExportPostfix '.txt'], 'wt');
        
        for iRow = 1:size(tempExport_hypn,1)
            fprintf(fid_export, '%i\t%i\n', tempExport_hypn(iRow,:));
        end
        fclose(fid_export);
        
        
        tempExport_hypn_index = (tempExport_hypn(:,2) ~= 3);
        tempExport_hypn_calc = tempExport_hypn(tempExport_hypn_index,:);
        tempExport_hypnStages_calc = tempExport_hypnStages(tempExport_hypn_index,:);
        
        totalSleepTime_export{iData} = size(tempExport_hypn_calc,1)*epochLength;
        
        S1Time_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S1')))*epochLength;
        S2Time_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S2')))*epochLength;
        S3Time_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S3')))*epochLength;
        S4Time_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S4')))*epochLength;
        REMtime_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'REM')))*epochLength;
        WakeTime_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'Wake')))*epochLength;
        MovementTime_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'MT')))*epochLength;
        SWStime_export{iData} = S3Time_export{iData} + S4Time_export{iData};
        NonREMtime_export{iData} = SWStime_export{iData} + S2Time_export{iData};
        
        
        S1Time_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S1') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        S2Time_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S2') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        S3Time_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S3') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        S4Time_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'S4') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        REMtime_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'REM') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        WakeTime_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'Wake') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        MovementTime_WithoutMA_export{iData} = length(find(strcmp(tempExport_hypnStages_calc(:,1),'MT') & (tempExport_hypn_calc(:,2) == 0) ))*epochLength;
        SWStime_WithoutMA_export{iData} = S3Time_WithoutMA_export{iData} + S4Time_WithoutMA_export{iData};
        NonREMtime_WithoutMA_export{iData} = SWStime_WithoutMA_export{iData} + S2Time_WithoutMA_export{iData};
        
        
    end
    
    
    
    
    
    
    
end
%open output files

if strcmp(ExportHypnogram,'yes')
    fidh = fopen([pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '_with_export' tempExportPostfix '.csv'],'wt');
else
    fidh = fopen([pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '.csv'],'wt');
end

%write header of ouptufiles


headerOrderSting = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n';

if strcmp(ExportHypnogram,'yes')
    headerOrderSting = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
end

fprintf(fidh,headerOrderSting,'datasetnum','dataset','hypnogram','epoch_length_seconds','Total_sleep_time_min','Sleep_Onset_min','S1_onset_min','S2_onset_min','SWS_onset_min','S4_onset_min','REM_onset_min'...
    ,'S1_min','S2_min','S3_min','S4_min','REM_min','Wake_after_sleep_onset_min','Movement_Time_min','SWS_min','NonREM_without_S1_min'...
    ,'S1_percent','S2_percent','S3_percent','S4_percent','REM_percent','Wake_after_sleep_onset_percent','Movement_Time_percent','SWS_percent','NonREM_without_S1_percent'...
    ,'S1_without_MA_min','S2_without_MA_min','S3_without_MA_min','S4_without_MA_min','REM_without_MA_min','Wake_after_sleep_onset_without_MA_min','Movement_Time_without_MA_min','SWS_without_MA_min','NonREM_without_S1_without_MA_min'...
    ,'S1_without_MA_percent','S2_without_MA_percent','S3_without_MA_percent','S4_without_MA_percent','REM_without_MA_percent','Wake_after_sleep_onset_without_MA_percent','Movement_Time_without_MA_percent','SWS_without_MA_percent','NonREM_without_S1_without_MA_percent'...
    ,'S1_before_sleep_onset_min','S2_before_sleep_onset_min','S3_before_sleep_onset_min','S4_before_sleep_onset_min','REM_before_sleep_onset_min','Wake_before_sleep_onset_min','Movement_before_sleep_onset_Time_min','SWS_before_sleep_onset_min','NonREM_before_sleep_onset_without_S1_min');


if strcmp(ExportHypnogram,'yes')
    fprintf(fidh,[',%s'...
    ',%s,%s,%s,%s,%s,%s,%s,%s,%s'...
    ',%s,%s,%s,%s,%s,%s,%s,%s,%s'...
    ',%s,%s,%s,%s,%s,%s,%s,%s,%s'...
    ',%s,%s,%s,%s,%s,%s,%s,%s,%s' '\n'],...
        'exported_Total_sleep_time_min'...
        ,'exported_S1_min','exported_S2_min','exported_S3_min','exported_S4_min','exported_REM_min','exported_Wake_after_sleep_onset_min','exported_Movement_Time_min','exported_SWS_min','exported_NonREM_without_S1_min'...
        ,'exported_S1_percent','exported_S2_percent','exported_S3_percent','exported_S4_percent','exported_REM_percent','exported_Wake_after_sleep_onset_percent','exported_Movement_Time_percent','exported_SWS_percent','exported_NonREM_without_S1_percent'...
        ,'exported_S1_without_MA_min','exported_S2_without_MA_min','exported_S3_without_MA_min','exported_S4_without_MA_min','exported_REM_without_MA_min','exported_Wake_after_sleep_onset_without_MA_min','exported_Movement_Time_without_MA_min','exported_SWS_without_MA_min','exported_NonREM_without_S1_without_MA_min'...
        ,'exported_S1_without_MA_percent','exported_S2_without_MA_percent','exported_S3_without_MA_percent','exported_S4_without_MA_percent','exported_REM_without_MA_percent','exported_Wake_after_sleep_onset_without_MA_percent','exported_Movement_Time_without_MA_percent','exported_SWS_without_MA_percent','exported_NonREM_without_S1_without_MA_percent');
end








for iData = iDatas
    fprintf(fidh,'%i,',iData);
    fprintf(fidh,'%s,',listOfDatasetsPaths{iData});
    fprintf(fidh,'%s,',listOfHypnogramPaths{iData});
    fprintf(fidh,'%f,',epochLength);
    fprintf(fidh,'%f,',totalSleepTime{iData}/60);
    fprintf(fidh,'%f,',sleepOnsetTime{iData}/60);
    fprintf(fidh,'%f,',S1OnsetTime{iData}/60);
    fprintf(fidh,'%f,',S2OnsetTime{iData}/60);
    fprintf(fidh,'%f,',SWSonsetTime{iData}/60);
    fprintf(fidh,'%f,',S4onsetTime{iData}/60);
    fprintf(fidh,'%f,',REMonsetTime{iData}/60);
    
    fprintf(fidh,'%f,',S1Time{iData}/60);
    fprintf(fidh,'%f,',S2Time{iData}/60);
    fprintf(fidh,'%f,',S3Time{iData}/60);
    fprintf(fidh,'%f,',S4Time{iData}/60);
    fprintf(fidh,'%f,',REMtime{iData}/60);
    fprintf(fidh,'%f,',WakeTime{iData}/60);
    fprintf(fidh,'%f,',MovementTime{iData}/60);
    fprintf(fidh,'%f,',SWStime{iData}/60);
    fprintf(fidh,'%f,',NonREMtime{iData}/60);
    fprintf(fidh,'%f,',100*S1Time{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S2Time{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S3Time{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S4Time{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*REMtime{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*WakeTime{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*MovementTime{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*SWStime{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*NonREMtime{iData}/totalSleepTime{iData});
    
    fprintf(fidh,'%f,',S1Time_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',S2Time_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',S3Time_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',S4Time_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',REMtime_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',WakeTime_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',MovementTime_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',SWStime_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',NonREMtime_WithoutMA{iData}/60);
    fprintf(fidh,'%f,',100*S1Time_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S2Time_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S3Time_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*S4Time_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*REMtime_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*WakeTime_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*MovementTime_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*SWStime_WithoutMA{iData}/totalSleepTime{iData});
    fprintf(fidh,'%f,',100*NonREMtime_WithoutMA{iData}/totalSleepTime{iData});
    
    fprintf(fidh,'%f,',S1TimePreOnset{iData}/60);
    fprintf(fidh,'%f,',S2TimePreOnset{iData}/60);
    fprintf(fidh,'%f,',S3TimePreOnset{iData}/60);
    fprintf(fidh,'%f,',S4TimePreOnset{iData}/60);
    fprintf(fidh,'%f,',REMtimePreOnset{iData}/60);
    fprintf(fidh,'%f,',WakeTimePreOnset{iData}/60);
    fprintf(fidh,'%f,',MovementTimePreOnset{iData}/60);
    fprintf(fidh,'%f,',SWStimePreOnset{iData}/60);
    if ~strcmp(ExportHypnogram,'yes')
        fprintf(fidh,'%f\n',NonREMtimePreOnset{iData}/60);
        
    else
        fprintf(fidh,'%f,',NonREMtimePreOnset{iData}/60);
        
        
        fprintf(fidh,'%f,',totalSleepTime_export{iData}/60);
        
        fprintf(fidh,'%f,',S1Time_export{iData}/60);
        fprintf(fidh,'%f,',S2Time_export{iData}/60);
        fprintf(fidh,'%f,',S3Time_export{iData}/60);
        fprintf(fidh,'%f,',S4Time_export{iData}/60);
        fprintf(fidh,'%f,',REMtime_export{iData}/60);
        fprintf(fidh,'%f,',WakeTime_export{iData}/60);
        fprintf(fidh,'%f,',MovementTime_export{iData}/60);
        fprintf(fidh,'%f,',SWStime_export{iData}/60);
        fprintf(fidh,'%f,',NonREMtime_export{iData}/60);
        fprintf(fidh,'%f,',100*S1Time_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S2Time_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S3Time_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S4Time_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*REMtime_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*WakeTime_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*MovementTime_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*SWStime_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*NonREMtime_export{iData}/totalSleepTime_export{iData});
        
        fprintf(fidh,'%f,',S1Time_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',S2Time_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',S3Time_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',S4Time_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',REMtime_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',WakeTime_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',MovementTime_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',SWStime_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',NonREMtime_WithoutMA_export{iData}/60);
        fprintf(fidh,'%f,',100*S1Time_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S2Time_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S3Time_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*S4Time_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*REMtime_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*WakeTime_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*MovementTime_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f,',100*SWStime_WithoutMA_export{iData}/totalSleepTime_export{iData});
        fprintf(fidh,'%f\n',100*NonREMtime_WithoutMA_export{iData}/totalSleepTime_export{iData});
        
    end
    
end
fclose(fidh);
if strcmp(ExportHypnogram,'yes')
    res_hypnvals = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '_with_export' tempExportPostfix '.csv'],'Delimiter',',');
else
    res_hypnvals = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '.csv'],'Delimiter',',');
end


fprintf('HypVals function finished\n');

toc
memtoc
end
