function res_nonEvents = spisop_nonEvents_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% discover correspondent non-events according to new events and hypnogram
% Copyright Frederik D. Weber

functionName = 'nonEvents';
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
EventsFilePathsFileName = getParam('EventsFilePathsFileName',listOfParameters);

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
if exist([pathInputFolder filesep EventsFilePathsFileName],'file') ~= 2
    error(['EventsFilePathsFileName file ' [pathInputFolder filesep EventsFilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end

EventsColumn = getParam('EventsColumn',listOfParameters);
ChannelColumn = getParam('ChannelColumn',listOfParameters);

MaxTrialsBeforeIncreaseSDsearchWindowArroundDetectedEvent = str2num(getParam('MaxTrialsBeforeIncreaseSDsearchWindowArroundDetectedEvent',listOfParameters));
MaxTrialsBeforeStopSearch = str2num(getParam('MaxTrialsBeforeStopSearch',listOfParameters));
SDsearchWindowArroundDetectedEventStepsize = str2num(getParam('SDsearchWindowArroundDetectedEventStepsize',listOfParameters)); % in seconds

PreEventBufferTime = str2num(getParam('PreEventBufferTime',listOfParameters)); % in seconds
PostEventBufferTime = str2num(getParam('PostEventBufferTime',listOfParameters)); % in seconds

PreBoundaryNonEventBufferTimeToEvents = str2num(getParam('PreBoundaryNonEventBufferTimeToEvents',listOfParameters)); % in seconds
PostBoundaryNonEventBufferTimeToEvents = str2num(getParam('PostBoundaryNonEventBufferTimeToEvents',listOfParameters)); % in seconds

NonOverlapConsideringAllChannels = getParam('NonOverlapConsideringAllChannels',listOfParameters);
SDsearchWindowArroundDetectedEvent = str2num(getParam('SDsearchWindowArroundDetectedEvent',listOfParameters)); % in seconds

RandomSeed = str2num(getParam('RandomSeed',listOfParameters));

%sleepStagesOfInterst = {'S3','S4'};
%sleepStagesOfInterst = {'SWS','S2'};
sleepStagesOfInterest = strsplit(getParam('sleepStagesOfInterest',listOfParameters));

rng(RandomSeed, 'twister');

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
listOfEventsFilePaths = read_mixed_csv([pathInputFolder filesep EventsFilePathsFileName],',');


if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramPaths)) && (size(listOfDatasetsPaths,1) == size(listOfEventsFilePaths,1)))
    error('files or number of Datasetspaths and Hypnogramsfiles are invalid or do not aggree')
end

iDatas = 1:(length(listOfDatasetsPaths));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end

AggregationOfDatasetOutputsOfDetections = getParam('AggregationOfDatasetOutputsOfDetections',listOfCoreParameters);%If the aggregation of datasetOutputfiles should be skipped either full or fast or no default full


core_cfg = [];
core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);

tic
memtic
fprintf('NonEvents function initialized\n');
conseciDatas = 1:length(iDatas);
parfor conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    
    datasetsPath = listOfDatasetsPaths{iData};
    hypnogramPath = listOfHypnogramPaths{iData};
    eventsPath = listOfEventsFilePaths{iData};
    
    dsEvents = dataset('File',eventsPath,'Delimiter',',');
    
    if any(isnumeric(dsEvents.(ChannelColumn)))
        dsEvents.(ChannelColumn) = cellstr(num2str(dsEvents.(ChannelColumn)));
    end
    
    events = dsEvents.(EventsColumn);
    channelsIndicator = dsEvents.(ChannelColumn);
    nEvents = length(events);
    
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
        tempdata = [];%clear
    else
        error('wrong parameter for IgnoreDataSetHeader either yes or no');
    end
    
    fprintf('dataset %i: process ROI from hypnogram info\n',iData);

    
    %sampleFreq = preDownsampleFreq;
    %epochLengthSamples = epochLength * preDownsampleFreq;
    [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLength);
    columnOfInterestHypnStages = 1;
    if ismember('SWS',sleepStagesOfInterest)
        columnOfInterestHypnStages = 2;
    elseif ismember('NonREM',sleepStagesOfInterest)
        columnOfInterestHypnStages = 3;
    end
    
    epochsOfInterst = hypnEpochs(ismember(hypnStages(:,columnOfInterestHypnStages),sleepStagesOfInterest) & (hypn(:,2) == 0));
    
    
    nEpochsOfInterest = length(epochsOfInterst);
    
    [consecBegins, consecEnds] = consecutiveBeginsAndEnds(epochsOfInterst,1);
    
    roiBegins = hypnEpochsBeginsSamples(consecBegins);
    roiEnds = hypnEpochsEndsSamples(consecEnds);
    hypnStagesOfInterst = hypnStages(consecBegins,:);
    
    if length(roiEnds) < 1
        error(['no ROI in data left for analysis']);
    end
    
    
    hypnEpochsBeginsSeconds = hypnEpochsBeginsSamples - 1;
    hypnEpochsEndsSeconds = hypnEpochsEndsSamples;
    
    roiBeginsSeconds = roiBegins - 1;
    roiEndsSeconds = roiEnds;
    
    ROILengthsSeconds = roiEndsSeconds - roiBeginsSeconds;
    
    indexROIsufficientLength = find(ROILengthsSeconds > (PreEventBufferTime + PostEventBufferTime));
    
    roiBeginsSeconds = roiBeginsSeconds(indexROIsufficientLength);
    roiEndsSeconds = roiEndsSeconds(indexROIsufficientLength);
    
    ROILengthsSeconds_sufficient = roiEndsSeconds - roiBeginsSeconds;
    
    hypnStagesOfInterst = hypnStagesOfInterst(indexROIsufficientLength);
    POI_length_cumsum = cumsum(ROILengthsSeconds_sufficient);
    
    nonEvent_maximum_nonInterfering_withinChannel_timepoints = sum(ROILengthsSeconds_sufficient)/((PreEventBufferTime + PostEventBufferTime)*2) - 2;
    if (nEvents > nonEvent_maximum_nonInterfering_withinChannel_timepoints)
        warning(['Space might not be sufficient to find enough non-overlapping non-events within even one channel in dataset ' num2str(iData)]);
    end
    
    channels = unique(channelsIndicator,'stable');
    ch_all_nonEvents = {};
    ori_event_row_index_indicator = [];
    iEvch = 0;
    
    ft_progress('init', 'text',    ['Dataset ' num2str(iData) ': Please wait...']);

    progress_count = 0;
    for iChan = 1:length(channels)
        %iChan = 1;
        ch = channels(iChan);
        %disp(iChan);
        ch_all_nonEvents{iChan} = [];
        candidate_compare_list = events(strcmp(channelsIndicator,ch));
        if strcmp(NonOverlapConsideringAllChannels,'yes')
            candidate_compare_list = events;
        end
        ev_ind = find(strcmp(channelsIndicator,ch));
        ev_ind_ind = randperm(length(ev_ind));
        ev_ind = ev_ind(ev_ind_ind);
        for iNE = 1:(length(ev_ind))
            %iNE = 2;
            iEv = ev_ind(iNE);
            iEvch = iEvch + 1;
            ori_event_row_index_indicator(iEvch) = iEv;
            %disp(iEv);
            progress_count = progress_count + 1;
            ft_progress(progress_count/nEvents, ['Dataset ' num2str(iData) ': Processing event %d of %d'], progress_count, nEvents);  % show string, x=i/N

            ev = events(iEv);
            evd = dsEvents(iEv,:);
            ev_epoch = hypnStages(((hypnEpochsBeginsSeconds <= ev) & (ev < hypnEpochsEndsSeconds)),columnOfInterestHypnStages);

            nonEventNonOverpapEvents_is_match = 0;
            nonEventNonOverpapNonEvents_is_match = 0;
            nonEvent_candidate_timepoint = -1;
            
               
            tempTryCounter = 1; 
            tempSD = SDsearchWindowArroundDetectedEvent;
            
            while (~nonEventNonOverpapEvents_is_match || ~nonEventNonOverpapNonEvents_is_match) && ~((tempTryCounter > MaxTrialsBeforeStopSearch) && nonEventNonOverpapEvents_is_match)
                detected_spindle_peak_pseudo = ev;
                
                if (detected_spindle_peak_pseudo > max(roiEndsSeconds))
                    timewindow_candidate_pre_index = length(roiEndsSeconds);
                    detected_spindle_peak_pseudo = max(roiEndsSeconds);
                elseif (detected_spindle_peak_pseudo < min(roiBeginsSeconds))
                    timewindow_candidate_pre_index = 1;
                    detected_spindle_peak_pseudo = min(roiBeginsSeconds);
                else
                    timewindow_candidate_pre_index = min(find(detected_spindle_peak_pseudo <= roiEndsSeconds));
                    
                    %timewindow_candidate_pre = hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index,]
                    
                    if(~((detected_spindle_peak_pseudo >= roiBeginsSeconds(timewindow_candidate_pre_index)) && (detected_spindle_peak_pseudo <= roiEndsSeconds(timewindow_candidate_pre_index))))
                        if ((detected_spindle_peak_pseudo - roiEndsSeconds(timewindow_candidate_pre_index-1) ) < ...
                                (roiBeginsSeconds(timewindow_candidate_pre_index) - detected_spindle_peak_pseudo))
                            timewindow_candidate_pre_index = timewindow_candidate_pre_index - 1;
                            detected_spindle_peak_pseudo = roiEndsSeconds(timewindow_candidate_pre_index);
                        else
                            %timewindow_candidate_pre_index = timewindow_candidate_pre_index
                            detected_spindle_peak_pseudo = roiBeginsSeconds(timewindow_candidate_pre_index);
                        end
                        
                    end
                    
                    if(timewindow_candidate_pre_index < 1)
                        timewindow_candidate_pre_index = 1;
                    end
                    %stop('error: timepoint min max = ',detected_spindle_peak_pseudo, ' ', min(roiBeginsSeconds), ' ', max(roiEndsSeconds) );
                end
                
                
                %timewindow_candidate_pre = hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index,]
                
                timewindow_gaps_cumulative_length = roiEndsSeconds(timewindow_candidate_pre_index) - POI_length_cumsum(timewindow_candidate_pre_index);
                
                corrected_center = detected_spindle_peak_pseudo - timewindow_gaps_cumulative_length;
                
                uncorrected_nonEvent_candidate_timepoint = -1;

                while((uncorrected_nonEvent_candidate_timepoint < 0) || (uncorrected_nonEvent_candidate_timepoint > max(POI_length_cumsum))  ) 
                    uncorrected_nonEvent_candidate_timepoint = corrected_center + tempSD*randn(1);
                end
                
                if mod(tempTryCounter,MaxTrialsBeforeIncreaseSDsearchWindowArroundDetectedEvent) == 0
                        tempSD = tempSD + SDsearchWindowArroundDetectedEventStepsize;
                end
                tempTryCounter = tempTryCounter + 1;
                
                
                timewindow_candidate_index = min(find(POI_length_cumsum >= uncorrected_nonEvent_candidate_timepoint));
                timewindow_gaps_cumulative_length_before_candidate = roiEndsSeconds(timewindow_candidate_index) - POI_length_cumsum(timewindow_candidate_index);
                
                
                nonEvent_candidate_timepoint = timewindow_gaps_cumulative_length_before_candidate + uncorrected_nonEvent_candidate_timepoint;
                
                overlapSpindels_index = find(...
                    ( (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (candidate_compare_list - PreBoundaryNonEventBufferTimeToEvents) & ...
                    (nonEvent_candidate_timepoint + PostEventBufferTime) <=  (candidate_compare_list + PostBoundaryNonEventBufferTimeToEvents) ) | ...
                    ( (nonEvent_candidate_timepoint - PreEventBufferTime) >=  (candidate_compare_list - PreBoundaryNonEventBufferTimeToEvents) &  ...
                    (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (candidate_compare_list + PostBoundaryNonEventBufferTimeToEvents) ) | ...
                    ( (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (candidate_compare_list - PreBoundaryNonEventBufferTimeToEvents) & ...
                    (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (candidate_compare_list + PostBoundaryNonEventBufferTimeToEvents) )...
                    );
                nonEventNonOverpapEvents_is_match = (length(overlapSpindels_index) == 0);
                
                overlapNonSpindels_index = find(...
                    ( (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (ch_all_nonEvents{iChan} - PreEventBufferTime) & ...
                    (nonEvent_candidate_timepoint + PostEventBufferTime) <=  (ch_all_nonEvents{iChan} + PostEventBufferTime) ) |...
                    ( (nonEvent_candidate_timepoint - PreEventBufferTime) >=  (ch_all_nonEvents{iChan} - PreEventBufferTime) & ...
                    (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (ch_all_nonEvents{iChan} + PostEventBufferTime) ) |...
                    ( (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (ch_all_nonEvents{iChan} - PreEventBufferTime) & ...
                    (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (ch_all_nonEvents{iChan} + PostEventBufferTime) )...
                    );
                nonEventNonOverpapNonEvents_is_match = (length(overlapNonSpindels_index) == 0);
            end
            ch_all_nonEvents{iChan} = [ch_all_nonEvents{iChan} nonEvent_candidate_timepoint];            
        end
    end
    ft_progress('close');
    
    
%     fide = fopen([pathOutputFolder filesep ouputFilesPrefixString 'non_events_' 'datanum_' num2str(iData) '.csv'],'wt');
%     
%     %write header of ouptufiles
%     fprintf(fide,'%s,%s,%s\n',...
%         'datasetnum', ChannelColumn, EventsColumn);
%     
%     
%     for iChan = 1:length(channels)
%         ch = channels{iChan};
%         %ch = 'C3';
%         
%         tempNewEvents = ch_all_nonEvents{iChan};
%         
%         for iLine=1:(length(tempNewEvents))
%             fprintf(fide,'%i,',iData);
%             fprintf(fide,'%s,',ch);
%             fprintf(fide,'%f\n',tempNewEvents(iLine));
%         end
%         
%     end
%     fclose(fide);
    
    
    varNames = {'datasetnum', ChannelColumn, EventsColumn};
    output = dataset([],[],[],'VarNames',varNames);
    for iChan = 1:length(channels)
        ch = channels{iChan};
        %ch = 'C3';      
        tempLengthEvents = length(ch_all_nonEvents{iChan});
        dummytempLengthEvents = tempLengthEvents;
        if tempLengthEvents == 1
            dummytempLengthEvents = 2;
            output = cat(1,output,dataset(repmat(iData,dummytempLengthEvents,1), repmat(ch,dummytempLengthEvents,1), repmat(ch_all_nonEvents{iChan},dummytempLengthEvents,1),'VarNames',varNames));
            output(size(output,1),:) = [];
        else
           output = cat(1,output,dataset(repmat(iData,dummytempLengthEvents,1), repmat(ch,dummytempLengthEvents,1), ch_all_nonEvents{iChan}','VarNames',varNames));
        end
        
    end

    dsEvents = set(dsEvents,'VarNames',strcat('ori_',get(dsEvents,'VarNames')));
    output = cat(2,output,dsEvents(ori_event_row_index_indicator,:));
    output(ori_event_row_index_indicator,:) = output(:,:);
    fprintf('Dataset %i: write results\n',iData);
    export(output,'file',[pathOutputFolder filesep ouputFilesPrefixString 'non_events_' 'datanum_' num2str(iData) '.csv'],'Delimiter',',');
end



%aggregate all results from datasets
fidf_all = [];
delimiter = ',';

if ~strcmp(AggregationOfDatasetOutputsOfDetections,'no')
fprintf('Aggregate results of all datasets\n');
for iData = iDatas

    fidf = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'non_events_' 'datanum_' num2str(iData) '.csv'],'Delimiter',delimiter);
    
    if iData == iDatas(1)
        fidf_all = fidf;
    else
        fidf_all = cat(1,fidf_all,fidf);
    end

end
export(fidf_all,'file',[pathOutputFolder filesep ouputFilesPrefixString 'non_events_' 'datanum_' 'all_recent' '.csv'],'Delimiter',delimiter);
end
res_nonEvents= fidf_all;




fprintf('NonEvents function finished\n');
toc
memtoc
end


% R template:
% 
% all_eeg_nonEvents_annotation_df = [];
% 
% for (p in all_eeg_spindles_annotation_df_split_by_part_list){
% 
%   
%   p = p[which(as_character(p$is_usable_for_analysis) == 'yes'),]
%   
%   p_split_by_channel = split(p, p$eeg_channel)
%   
%   all_eeg_nonEvents_annotation_df_pc = data_frame();
%   for (pc in p_split_by_channel){
%     %pc = p_split_by_channel[[1]]
%     candidate_compare_list = pc
%     if (NonOverlapConsideringAllChannels){
%       candidate_compare_list = p
%     }
%     
%     pc_1 = pc[1,]
%     p_corresp_epoch_index = which(as_character(all_eeg_hypn_df$vpn) == as_character(pc_1$vpn) & 
%                                     as_character(all_eeg_hypn_df$nightnumber) == as_character(pc_1$nightnumber) & 
%                                     as_character(all_eeg_hypn_df$partnumber) == as_character(pc_1$partnumber))
%     
%     hypn_candidate_epochs = all_eeg_hypn_df[p_corresp_epoch_index,]
%     for (i in 1:nrow(pc)){
%       %i = 1
%       cat(' ');cat(i);
%       se = pc[i,]
%       
%       se_corresp_epoch_index = which(as_character(hypn_candidate_epochs$stage_alt) == as_character(se$stage_alt) & 
%                                        as_character(hypn_candidate_epochs$is_usable_for_analysis) == as_character(se$is_usable_for_analysis))
%       
%       hypn_candidate_epochs_select = hypn_candidate_epochs[se_corresp_epoch_index,]
%       
%       consec_epoch = 1
%       consec_epoch_idicator = c();
%       for (e in 1:nrow(hypn_candidate_epochs_select)){
%         if ((e > 1) && (hypn_candidate_epochs_select[e-1,]$eeg_epoch + 1 != hypn_candidate_epochs_select[e,]$eeg_epoch)){
%             consec_epoch = consec_epoch + 1
%         }
%         consec_epoch_idicator = c(consec_epoch_idicator,consec_epoch)  
%       }
%       hypn_candidate_epochs_select_ext = cbind(hypn_candidate_epochs_select, consec_epoch_idicator)
%     
%       hypn_candidate_epochs_select_timewindows = cbind(aggregate(eeg_epoch_length_seconds ~ consec_epoch_idicator, data = hypn_candidate_epochs_select_ext, sum),
%             epoch_start_seconds = aggregate(epoch_start_seconds ~ consec_epoch_idicator, data = hypn_candidate_epochs_select_ext,min)$epoch_start_seconds,
%             epoch_end_seconds = aggregate(epoch_end_seconds ~ consec_epoch_idicator, data = hypn_candidate_epochs_select_ext,max)$epoch_end_seconds
%             )
%       hypn_candidate_epochs_select_timewindows = hypn_candidate_epochs_select_timewindows[which(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds > (PreEventBufferTime + PostEventBufferTime)),]
%       hypn_candidate_epochs_select_timewindows = within(hypn_candidate_epochs_select_timewindows, {eeg_epoch_length_seconds = eeg_epoch_length_seconds - PreEventBufferTime - PostEventBufferTime
%                                                          epoch_start_seconds = epoch_start_seconds + PreEventBufferTime
%                                                          epoch_end_seconds = epoch_end_seconds - PostEventBufferTime})
%       hypn_candidate_epochs_select_timewindows = cbind(hypn_candidate_epochs_select_timewindows,eeg_epoch_length_seconds_cumsum = cumsum(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds))
%       
%       nonEvent_maximum_nonInterfering_withinChannel_timepoints = sum(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds)/((PreEventBufferTime + PostEventBufferTime)*2) - 2
%       
%       nonEventNonOverpapEvents_is_match = F;
%       nonEventNonOverpapNonEvents_is_match = F;
%       nonEvent_candidate_timepoint = -1;
%       while(!nonEventNonOverpapEvents_is_match || !nonEventNonOverpapNonEvents_is_match){
%         detected_spindle_peak_pseudo = unfactor(se$Timepoints_of_Event_Peaks,as_numeric);
%         
%           if (detected_spindle_peak_pseudo > max(roiEndsSeconds)){
%             timewindow_candidate_pre_index = nrow(hypn_candidate_epochs_select_timewindows);
%             detected_spindle_peak_pseudo = max(roiEndsSeconds)
%           } else if (detected_spindle_peak_pseudo < min(roiBeginsSeconds)) {
%             timewindow_candidate_pre_index = 1;
%             detected_spindle_peak_pseudo = min(roiBeginsSeconds)
%             
%           } else {
%             timewindow_candidate_pre_index = min(which(detected_spindle_peak_pseudo <= roiEndsSeconds));
%             
%             timewindow_candidate_pre = hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index,]
%             
%             if(!((detected_spindle_peak_pseudo >= timewindow_candidate_pre$epoch_start_seconds) && (detected_spindle_peak_pseudo <= timewindow_candidate_pre$epoch_end_seconds))){
%               if ((detected_spindle_peak_pseudo - hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index-1,]$epoch_end_seconds ) <
%                 (timewindow_candidate_pre$epoch_start_seconds - detected_spindle_peak_pseudo)){
%                 timewindow_candidate_pre_index = timewindow_candidate_pre_index - 1
%                 detected_spindle_peak_pseudo = hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index,]$epoch_end_seconds
%               } else {
%                 %timewindow_candidate_pre_index = timewindow_candidate_pre_index
%                 detected_spindle_peak_pseudo = timewindow_candidate_pre$epoch_start_seconds
%               }
%               
%             }
%             
%             if(timewindow_candidate_pre_index < 1){
%               timewindow_candidate_pre_index = 1
%             }
%             %stop('error: timepoint min max = ',detected_spindle_peak_pseudo, ' ', min(roiBeginsSeconds), ' ', max(roiEndsSeconds) );
%           }
% 
%         
%         timewindow_candidate_pre = hypn_candidate_epochs_select_timewindows[timewindow_candidate_pre_index,]
% 
%         timewindow_gaps_cumulative_length = timewindow_candidate_pre$epoch_end_seconds - timewindow_candidate_pre$eeg_epoch_length_seconds_cumsum 
%         %timewindow_gaps_cumulative_length = timewindow_candidate_pre$epoch_start_seconds - hypn_candidate_epochs_select_timewindows[max(which(roiEndsSeconds <= timewindow_candidate$epoch_start_seconds)),]$eeg_epoch_length_seconds_cumsum         
%        
%         
%         corrected_center = detected_spindle_peak_pseudo - timewindow_gaps_cumulative_length
%         
%         uncorrected_nonEvent_candidate_timepoint = -1
%         while((uncorrected_nonEvent_candidate_timepoint < 0) || (uncorrected_nonEvent_candidate_timepoint > max(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds_cumsum)) ){
%           uncorrected_nonEvent_candidate_timepoint = rnorm(1, mean = corrected_center, sd = SDsearchWindowArroundDetectedEvent)
%         }
%         %timewindow_candidate_index = runif(n = 1, min = 0, max = sum(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds))
%         %uncorrected_nonEvent_candidate_timepoint
%         
%         timewindow_candidate = hypn_candidate_epochs_select_timewindows[min(which(hypn_candidate_epochs_select_timewindows$eeg_epoch_length_seconds_cumsum >= uncorrected_nonEvent_candidate_timepoint)),]
%         timewindow_gaps_cumulative_length_before_candidate = timewindow_candidate$epoch_end_seconds - timewindow_candidate$eeg_epoch_length_seconds_cumsum 
% 
%         %timewindow_gaps_cumulative_length_before_candidate = timewindow_candidate$epoch_start_seconds - sum(hypn_candidate_epochs_select_timewindows[which(roiEndsSeconds <= timewindow_candidate$epoch_start_seconds),]$eeg_epoch_length_seconds_cumsum) 
%         
%         nonEvent_candidate_timepoint = timewindow_gaps_cumulative_length_before_candidate + uncorrected_nonEvent_candidate_timepoint
%         %nonEvent_candidate_timepoint = runif(n = 1, min = timewindow_candidate$epoch_start_seconds, max = timewindow_candidate$epoch_end_seconds)
%         
%         overlapSpindels_index = which(
%           ( (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) - PreBoundaryNonEventBufferTimeToEvents) & 
%             (nonEvent_candidate_timepoint + PostEventBufferTime) <=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) + PostBoundaryNonEventBufferTimeToEvents) ) |
%             ( (nonEvent_candidate_timepoint - PreEventBufferTime) >=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) - PreBoundaryNonEventBufferTimeToEvents) & 
%                 (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) + PostBoundaryNonEventBufferTimeToEvents) ) |
%             ( (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) - PreBoundaryNonEventBufferTimeToEvents) & 
%                 (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (unfactor(candidate_compare_list$Timepoints_of_Event_Peaks,as_numeric) + PostBoundaryNonEventBufferTimeToEvents) )
%           ) 
%         nonEventNonOverpapEvents_is_match = length(overlapSpindels_index) == 0
%         %candidate_compare_list[overlapSpindels_index,]
%         
%         overlapNonSpindels_index = which(
%           ( (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) - PreEventBufferTime) & 
%               (nonEvent_candidate_timepoint + PostEventBufferTime) <=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) + PostEventBufferTime) ) |
%             ( (nonEvent_candidate_timepoint - PreEventBufferTime) >=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) - PreEventBufferTime) & 
%                 (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) + PostEventBufferTime) ) |
%             ( (nonEvent_candidate_timepoint - PreEventBufferTime) <=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) - PreEventBufferTime) & 
%                 (nonEvent_candidate_timepoint + PostEventBufferTime) >=  (unfactor(all_eeg_nonEvents_annotation_df_pc$Timepoints_of_Event_Peaks,as_numeric) + PostEventBufferTime) )
%             ) 
%         nonEventNonOverpapNonEvents_is_match = length(overlapSpindels_index) == 0
%       }
%       
%       se_nonEventSegments = se
%       se_nonEventSegments$Timepoints_of_Event_Peaks = nonEvent_candidate_timepoint
%       
%       all_eeg_nonEvents_annotation_df_pc = rbind(all_eeg_nonEvents_annotation_df_pc,se_nonEventSegments)
%     }
%   }
%   all_eeg_nonEvents_annotation_df = rbind(all_eeg_nonEvents_annotation_df,all_eeg_nonEvents_annotation_df_pc)
% }
% 
% 
% if (doWriteFiles){
%   write_csv(all_eeg_nonEvents_annotation_df,file = paste0('/data3_ext/HamburgMEG/','all_eeg_non_spindles_annotation_csv'),row_names = F)
% }
% 
% all_eeg_nonEvents_annotation_df_split_by_part_list = split(all_eeg_nonEvents_annotation_df, all_eeg_nonEvents_annotation_df$filepath)
% if (doWriteFiles){
%   for (p in all_eeg_nonEvents_annotation_df_split_by_part_list){
%     write_csv(p,file = paste0('/data3_ext/HamburgMEG/',p$filepath[1],'/',p$subject[1],'_',p$vpn[1],'_',p$nightnumber[1],'_',p$partnumber[1],'_',p$frequencyInHz[1],'_','all_nonEvents_in_eeg_channels_and_annotation_csv'), row_names = F)
%   }
% }
% 
% 
% all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM = all_eeg_nonEvents_annotation_df[which(all_eeg_nonEvents_annotation_df$is_usable_for_analysis == 'yes' & all_eeg_nonEvents_annotation_df$stage_alt2 == 'NonREM' & all_eeg_nonEvents_annotation_df$id %in% all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]
% 
% all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM = all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]
% 
% 
% all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM_split_by_part_list = split(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM, all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$filepath)
% if (doWriteFiles){
%   for (p in all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM_split_by_part_list){
%     p_split_by_channel = split(p, p$eeg_channel)
%     write_csv(p,file = paste0('/data3_ext/HamburgMEG/',p$filepath[1],'/',p$subject[1],'_',p$vpn[1],'_',p$nightnumber[1],'_',p$partnumber[1],'_',p$frequencyInHz[1],'_','all_usable_nonEvents_in_eeg_channels_csv'), row_names = F)
%     for (pc in p_split_by_channel){
%       write_csv(pc,file = paste0('/data3_ext/HamburgMEG/',pc$filepath[1],'/',pc$subject[1],'_',pc$vpn[1],'_',pc$nightnumber[1],'_',pc$partnumber[1],'_',pc$frequencyInHz[1],'_','all_usable_nonEvents_in_eeg_channel_',pc$eeg_channel[1],'_csv'), row_names = F)
%     }
%   }
% }
% %test for matching with spindles
% all(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$id == 
%   all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$id)
% 
% 
% hist(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks -
%       unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric),
%      breaks=18000,xlim=c(-100,100))
% 
% qqnorm(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks -
%        unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric))
% sd(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks -
%          unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric))
% 
% length(which(abs(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks -
%      unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric)) > 60))
% 
% shapiro_test(
%   sample((all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks -
%          unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric)),5000))
% 
% hist(unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric)
%      ,breaks=100%,xlim=c(-500,500)
%      )
% hist(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks
%      ,breaks=100%,xlim=c(-500,500)
% )
% 
% qqplot(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,
%        unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric))
% 
% %length(unique(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$id))/
% %length(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$id)  
% 
% boxplot(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_nonEvents_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,
%         unfactor(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM[order(all_eeg_spindles_annotation_df_is_usable_for_analysis_and_NonREM$id),]$Timepoints_of_Event_Peaks,as_numeric))
% 
