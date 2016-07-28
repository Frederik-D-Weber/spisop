function [res_match, res_mismatch, res_summary] = spisop_eventCooccurrence_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfParameters)
% discover co-occurrence of test and target events, i.e. if test events
% fall within a defined timewindow arround target events.
% Copyright Frederik D. Weber

EventsTestFilePathsFileName = getParam('EventsTestFilePathsFileName',listOfParameters);
EventsTargetFilePathsFileName = getParam('EventsTargetFilePathsFileName',listOfParameters);

if exist([pathInputFolder filesep EventsTestFilePathsFileName],'file') ~= 2
    error(['EventsTestFilePathsFileName file ' [pathInputFolder filesep EventsTestFilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
if exist([pathInputFolder filesep EventsTargetFilePathsFileName],'file') ~= 2
    error(['EventsTargetFilePathsFileName file ' [pathInputFolder filesep EventsTargetFilePathsFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end

EventsTestTimePointColumn = getParam('EventsTestTimePointColumn',listOfParameters);
EventsTargetTimePointColumn = getParam('EventsTargetTimePointColumn',listOfParameters);

EventsTestCompareColumns = strsplit(getParam('EventsTestCompareColumns',listOfParameters),' ');
EventsTargetCompareColumns = strsplit(getParam('EventsTargetCompareColumns',listOfParameters),' ');

EventsTestGroupSummaryByColumns = strsplit(getParam('EventsTestGroupSummaryByColumns',listOfParameters),' ');

EventTargetTimeWindowOffsetTime = str2num(getParam('EventTargetTimeWindowOffsetTime',listOfParameters)); % in units of EventsTargetTimePointColumn
EventTargetTimeWindowPreOffsetTime = str2num(getParam('EventTargetTimeWindowPreOffsetTime',listOfParameters)); % in units of EventsTargetTimePointColumn
EventTargetTimeWindowPostOffsetTime = str2num(getParam('EventTargetTimeWindowPostOffsetTime',listOfParameters)); % in units of EventsTargetTimePointColumn
EventTestTimePointOffsetTime = str2num(getParam('EventTestTimePointOffsetTime',listOfParameters)); % in units of EventsTestTimePointColumn

UseSecondColumnAndOnlyOffsetsForTimeWindow = getParam('UseSecondColumnAndOnlyOffsetsForTimeWindow',listOfParameters);
EventsTargetTimePointColumn2 = getParam('EventsTargetTimePointColumn2',listOfParameters);
EventTargetTimeWindowOffsetTime2 = str2num(getParam('EventTargetTimeWindowOffsetTime2',listOfParameters)); % in units of EventsTargetTimePointColumn


EventsFilesWhich = getParam('EventsFilesWhich',listOfParameters);%Event files to be processed either all or subset if subset then DataSetsNumbers is used for selection default all
EventsFilesNumbers = str2num(getParam('EventsFilesNumbers',listOfParameters));%The line numbers of the events file to be processed if EventsFilesWhich parameter is set to subset


listOfEventsTestPaths = read_mixed_csv([pathInputFolder filesep EventsTestFilePathsFileName],',');
listOfEventsTargetPaths = read_mixed_csv([pathInputFolder filesep EventsTargetFilePathsFileName],',');

if ~(all(size(EventsTestCompareColumns) == size(EventsTargetCompareColumns)))
    error('number of test events and target events columns do not aggree')
end

if ~(all(size(listOfEventsTestPaths) == size(listOfEventsTargetPaths)))
    error('files or number of test events and target events files paths are invalid or do not aggree')
end

iDatas = 1:(length(listOfEventsTestPaths));

if strcmp(EventsFilesWhich,'subset')
    if ~(ismember(min(EventsFilesNumbers),iDatas) && ismember(max(EventsFilesNumbers),iDatas))
        error('Parameter EventsFilesNumbers contains numbers not matching to any line number, e.g. too less EventsTestFilePaths in EventsTestFilePathsFileName!')
    end
    iDatas = EventsFilesNumbers;
end


FilterValuesSplitString = ' ';
EventsTestFilterForColumn = getParam('EventsTestFilterForColumn',listOfParameters);%variable name for event value in test events files to apply text filter to if nothing is entered it is not filtered e.g. channel default is no value entered 
EventsTargetFilterForColumn = getParam('EventsTargetFilterForColumn',listOfParameters);%variable name for event value in test events files to apply text filter to if nothing is entered it is not filtered e.g. channel default is no value entered
EventsTestFilterValues = strsplit(getParam('EventsTestFilterValues',listOfParameters),FilterValuesSplitString);%variable values for EventsTestFilterForColum in test events files to apply text filter e.g. Cz
EventsTargetFilterValues = strsplit(getParam('EventsTargetFilterValues',listOfParameters),FilterValuesSplitString);%variable values for EventsTargetFilterForColum in test events files to apply text filter e.g. Cz


IDColumnsSplitString = ' ';

MismatchIdenticalEvents = getParam('MismatchIdenticalEvents',listOfParameters);%either yes no default no
MismatchDuplicateTestTargetOrTargetTestMatchingEvents = getParam('MismatchDuplicateTestTargetOrTargetTestMatchingEvents',listOfParameters);%either yes no default no
EventsTestIDColumns =  strsplit(getParam('EventsTestIDColumns',listOfParameters),IDColumnsSplitString);
EventsTargetIDColumns = strsplit(getParam('EventsTargetIDColumns',listOfParameters),IDColumnsSplitString);

ChunkBufferSize = getParam('ChunkBufferSize',listOfParameters);%size for buffer to not let datasets get to big for concatenation. this give the maximal number of lines a dataset can have for concatenation and when low can increase performance for >10000 events. default 100


GroupByConcatString = '#%#';

tic
memtic

fprintf('EventCooccurrence function initialized\n');

overlaps = [];
nonoverlaps = [];
summarys = [];
for iData = iDatas
    overlaps{iData} = [];
    nonoverlaps{iData} = [];
    summarys{iData} = [];
end

conseciDatas = 1:length(iDatas);
parfor conseciData = conseciDatas
    iData = iDatas(conseciData);
    %iData = 1
    
    
    eventsTestPath = listOfEventsTestPaths{iData};
    eventsTargetPath = listOfEventsTargetPaths{iData};
    
    
    dsEventsTest = dataset('File',eventsTestPath,'Delimiter',',');
    dsEventsTarget = dataset('File',eventsTargetPath,'Delimiter',',');
    
    if ~isempty(EventsTestFilterForColumn)
        matchIndicator = zeros(size(dsEventsTest,1),1);
         for iComb = 1:length(EventsTestFilterValues)
            %iComp = 1
            tempCompTest = EventsTestFilterValues{iComb};
            %tempCompTarget = EventsTargetFilterValues{iComb};
            
            if iscell(dsEventsTest.(EventsTestFilterForColumn))
                matchIndicator = matchIndicator | ( strcmp(dsEventsTest.(EventsTestFilterForColumn), tempCompTest) );
            else
                matchIndicator = matchIndicator | ( dsEventsTest.(EventsTestFilterForColumn) ==  tempCompTest);
            end
        end
        
        dsEventsTest = dsEventsTest(matchIndicator,:);
    end
    
    if ~isempty(EventsTargetFilterForColumn)
        matchIndicator = zeros(size(dsEventsTarget,1),1);
         for iComb = 1:length(EventsTargetFilterValues)
            %iComp = 1
            %tempCompTest = EventsTestFilterValues{iComb};
            tempCompTarget = EventsTargetFilterValues{iComb};
            
            if iscell(dsEventsTarget.(EventsTargetFilterForColumn))
                matchIndicator = matchIndicator | ( strcmp(dsEventsTarget.(EventsTargetFilterForColumn), tempCompTarget) );
            else
                matchIndicator = matchIndicator | ( dsEventsTarget.(EventsTargetFilterForColumn) ==  tempCompTarget);
            end
        end
        
        dsEventsTarget = dsEventsTarget(matchIndicator,:);
    end
    

    
    nEventsTest = size(dsEventsTest,1);
    nEventsTarget = size(dsEventsTarget,1);
    
    groupByMapOverlap = containers.Map();
    groupByMapNonOverlap = containers.Map();
    
    groupByMapAllTest = containers.Map();
    groupByMapAllTarget = containers.Map();
    
    if strcmp(MismatchIdenticalEvents,'yes') || strcmp(MismatchDuplicateTestTargetOrTargetTestMatchingEvents,'yes')
        
        IDmergeString = '#';
        temp_test_id = {''};
        temp_target_id = {''};
         for iComb = 1:numel(EventsTestIDColumns)
            tempCompTest = EventsTestIDColumns{iComb};
            tempCompTarget = EventsTargetIDColumns{iComb};
            
            %dsEventsTest.(tempCompTest)
            %dsEventsTarget.(tempCompTarget)

            if iscell(dsEventsTest.(tempCompTest))
                temp_test_id = strcat(temp_test_id,{IDmergeString},dsEventsTest.(tempCompTest));
            else
                temp_test_id = strcat(temp_test_id,{IDmergeString},num2str(dsEventsTest.(tempCompTest),'%-g'));
            end
            
            if iscell(dsEventsTarget.(tempCompTarget))
                temp_target_id = strcat(temp_target_id,{IDmergeString},dsEventsTarget.(tempCompTarget));
            else
                temp_target_id = strcat(temp_target_id,{IDmergeString},num2str(dsEventsTarget.(tempCompTarget),'%-g'));
            end
         end
         
         dsEventsTest.duplication_id = temp_test_id;
         dsEventsTarget.duplication_id = temp_target_id;
         
         temp_test_id = [];
         temp_target_id = [];
         
         duplicateIDMapTest = containers.Map();
         duplicateIDMapTarget = containers.Map();
        
    end
    
    overlap = [];
    nonoverlap =[];
    
    
    
    temp_overlap_collector_iterator = 1;
    temp_overlap_collector = {};
    temp_overlap_collector{temp_overlap_collector_iterator} = [];
    
    temp_nonoverlap_collector_iterator = 1;
    temp_nonoverlap_collector = {};
    temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = [];
    
    ft_progress('init', 'text',    ['EventTestFile ' num2str(iData) ': Please wait...']);
    
    columnNamesTestNew = strcat('test_',get(dsEventsTest,'VarNames'));
    columnNamesTargetNew = strcat('target_',get(dsEventsTarget,'VarNames'));
    
    progress_count = 0;
    for iEvTest = 1:nEventsTest
        %iEvTest = 1
        
        progress_count = progress_count + 1;
        ft_progress(progress_count/nEventsTest, ['EventTestFile ' num2str(iData) ': Processing test event %d (matchchunk %d, mismatchchunk %d) of %d against %d targets'], progress_count,temp_overlap_collector_iterator, temp_nonoverlap_collector_iterator, nEventsTest, nEventsTarget);  % show string, x=i/N
        
        
        eventTest = dsEventsTest(iEvTest,:);
        matchIndicator = ones(nEventsTarget,1);
        
        
        if strcmp(MismatchIdenticalEvents,'yes')
        	matchIndicator = matchIndicator & ~(strcmp(eventTest.duplication_id,dsEventsTarget.duplication_id));
        end
        
        
        for iComb = 1:length(EventsTestCompareColumns)
            %iComp = 1
            tempCompTest = EventsTestCompareColumns{iComb};
            tempCompTarget = EventsTargetCompareColumns{iComb};
            
            if iscell(eventTest.(tempCompTest))
                matchIndicator = matchIndicator & ( strcmp(eventTest.(tempCompTest),dsEventsTarget.(tempCompTarget)) );
            else
                matchIndicator = matchIndicator & ( eventTest.(tempCompTest) == dsEventsTarget.(tempCompTarget) );
            end
        end
        
        
        groupBy = 'group';
        for iGroup = 1:length(EventsTestGroupSummaryByColumns)
            %iComp = 1
            tempCompTest = EventsTestGroupSummaryByColumns{iGroup};
     
            if iscell(eventTest.(tempCompTest))
                groupBy = [groupBy GroupByConcatString eventTest.(tempCompTest){:}];
            else
                groupBy = [groupBy GroupByConcatString num2str(eventTest.(tempCompTest))];
            end
        end
        
        if ~isKey(groupByMapOverlap,{groupBy})
            groupByMapOverlap(groupBy) = 0;
        end
        
        if ~isKey(groupByMapNonOverlap,{groupBy})
            groupByMapNonOverlap(groupBy) = 0;
        end
        
        if ~isKey(groupByMapAllTest,{groupBy})
            groupByMapAllTest(groupBy) = 0;
        end
        
        
        groupByMapAllTarget(groupBy) = sum(matchIndicator);
        
        groupByMapAllTest(groupBy) = groupByMapAllTest(groupBy) + 1;
        
        if strcmp(UseSecondColumnAndOnlyOffsetsForTimeWindow,'yes')
            matchIndicator = matchIndicator & ( (eventTest.(EventsTestTimePointColumn) + EventTestTimePointOffsetTime ) >= (dsEventsTarget.(EventsTargetTimePointColumn) + EventTargetTimeWindowOffsetTime) );
            matchIndicator = matchIndicator & ( (eventTest.(EventsTestTimePointColumn) + EventTestTimePointOffsetTime ) <= (dsEventsTarget.(EventsTargetTimePointColumn2) + EventTargetTimeWindowOffsetTime2) );
        else
            matchIndicator = matchIndicator & ( (eventTest.(EventsTestTimePointColumn) + EventTestTimePointOffsetTime ) <= ((dsEventsTarget.(EventsTargetTimePointColumn) + EventTargetTimeWindowOffsetTime) + EventTargetTimeWindowPostOffsetTime) );
            matchIndicator = matchIndicator & ( (eventTest.(EventsTestTimePointColumn) + EventTestTimePointOffsetTime ) >= ((dsEventsTarget.(EventsTargetTimePointColumn) + EventTargetTimeWindowOffsetTime) - EventTargetTimeWindowPreOffsetTime) );
        end
        %matchIndicator(2:3) = 1;
        

        eventTest = set(eventTest,'VarNames',columnNamesTestNew);

        if any(matchIndicator)
            tempDsEventsTarget = dsEventsTarget(matchIndicator,:);
            if strcmp(MismatchDuplicateTestTargetOrTargetTestMatchingEvents,'yes')
                
                temp_test_target_id_strings = strcat(eventTest.(['test_' 'duplication_id']),dsEventsTarget.duplication_id(matchIndicator));
                temp_target_test_id_strings = strcat(dsEventsTarget.duplication_id(matchIndicator),eventTest.(['test_' 'duplication_id']));
                
                temp_already_contained_test_first = isKey(duplicateIDMapTest,temp_test_target_id_strings);
                temp_already_contained_target_first = isKey(duplicateIDMapTarget,temp_target_test_id_strings);
                
                temp_add_to_overlapp_notcontainted_before_index = ~( temp_already_contained_test_first | temp_already_contained_target_first );

                if any(~temp_already_contained_test_first)
                    temp_add_to_map = temp_test_target_id_strings(~temp_already_contained_test_first);
                    duplicateIDMapTest = [duplicateIDMapTest;containers.Map(temp_add_to_map,ones(numel(temp_add_to_map),1))];
                end
                if any(~temp_already_contained_target_first)
                    temp_add_to_map = temp_target_test_id_strings(~temp_already_contained_target_first);
                    duplicateIDMapTarget = [duplicateIDMapTarget;containers.Map(temp_add_to_map,ones(numel(temp_add_to_map),1))];
                end
                tempDsEventsTarget = tempDsEventsTarget(temp_add_to_overlapp_notcontainted_before_index,:);
            end
            
            tempDsEventsTarget = set(tempDsEventsTarget,'VarNames',columnNamesTargetNew);
            nOverlaps = size(tempDsEventsTarget,1);
            for iTarRow = 1:nOverlaps
                
                if size(temp_overlap_collector{temp_overlap_collector_iterator},1) > 0
                %if size(overlap,1) > 0
                    temp_overlap_collector{temp_overlap_collector_iterator}  = cat(1,temp_overlap_collector{temp_overlap_collector_iterator} ,cat(2,eventTest,tempDsEventsTarget(iTarRow,:)));
                    %overlap = cat(1,overlap,cat(2,eventTest,tempDsEventsTarget(iTarRow,:)));
                else
                    temp_overlap_collector{temp_overlap_collector_iterator} = cat(2,eventTest,tempDsEventsTarget(iTarRow,:));
                    %overlap = cat(2,eventTest,tempDsEventsTarget(iTarRow,:));
                end
                if size(temp_overlap_collector{temp_overlap_collector_iterator},1) > ChunkBufferSize
                    temp_overlap_collector_iterator = temp_overlap_collector_iterator + 1;
                    temp_overlap_collector{temp_overlap_collector_iterator} = [];
                end
                
            end
            groupByMapOverlap(groupBy) = groupByMapOverlap(groupBy) + nOverlaps;
        else
            
            if size(temp_nonoverlap_collector{temp_nonoverlap_collector_iterator},1) > 0
            %if size(nonoverlap,1) > 0
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator}  = cat(1,temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} ,eventTest);
                %nonoverlap = cat(1,nonoverlap,eventTest);
            else
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = eventTest;
                %nonoverlap = eventTest;
            end
            if size(temp_nonoverlap_collector{temp_nonoverlap_collector_iterator},1) > ChunkBufferSize
                temp_nonoverlap_collector_iterator = temp_nonoverlap_collector_iterator + 1;
                temp_nonoverlap_collector{temp_nonoverlap_collector_iterator} = [];
            end
            
            groupByMapNonOverlap(groupBy) = groupByMapNonOverlap(groupBy) + 1;
        end
    end
    
    for iTemp_overlap_collector_iterator = 1:numel(temp_overlap_collector)
        if iTemp_overlap_collector_iterator == 1
            overlap = temp_overlap_collector{iTemp_overlap_collector_iterator};
        else
            overlap = cat(1,overlap,temp_overlap_collector{iTemp_overlap_collector_iterator});
        end
    end
    temp_overlap_collector = [];
    
    
      for iTemp_nonoverlap_collector_iterator = 1:numel(temp_nonoverlap_collector)
        if iTemp_nonoverlap_collector_iterator == 1
            nonoverlap = temp_nonoverlap_collector{iTemp_nonoverlap_collector_iterator};
        else
            nonoverlap = cat(1,nonoverlap,temp_nonoverlap_collector{iTemp_nonoverlap_collector_iterator});
        end
    end
    temp_nonoverlap_collector = [];
    
    
    
%     if isempty(nonoverlap)
%         for iDScol = 1:length(columnNamesTestNew)
%             nonoverlap = cat(2,nonoverlap,dataset([],'VarNames',{num2str(iDScol)}));
%         end
%         nonoverlap = set(nonoverlap,'VarNames',columnNamesTestNew);
%     end
%     
%     if isempty(overlap)
%         for iDScol = 1:length([columnNamesTestNew,columnNamesTargetNew])
%             overlap = cat(2,overlap,dataset([],'VarNames',{num2str(iDScol)}));
%         end
%         overlap = set(overlap,'VarNames',[columnNamesTestNew,columnNamesTargetNew]);
%     end
    
    ft_progress('close');  
    varNames = {'used_2nd_column_and_offset_not_pre_and_post_and_offsets','test_compare_columns', 'target_compare_columns','test_timepoint_column', 'target_timepoint_column', 'target_offset', 'target_pre_offset', 'target_post_offset','test_offset'};
    
    tempEventsTestCompareColumns = EventsTestCompareColumns;
    tempEventsTargetCompareColumns = EventsTargetCompareColumns;
    
    if size(tempEventsTestCompareColumns,2) > 1
        tempEventsTestCompareColumns = {strjoin(EventsTestCompareColumns,' ')};
        tempEventsTargetCompareColumns = {strjoin(EventsTargetCompareColumns,' ')};
    end
    
 
    %fprintf(['error at dataset' num2str(iData) '\n']);
%     if true
        addRightOverlap = dataset(repmat({UseSecondColumnAndOnlyOffsetsForTimeWindow},size(overlap,1),1),repmat(tempEventsTestCompareColumns,size(overlap,1),1),...
            repmat(tempEventsTargetCompareColumns,size(overlap,1),1),...
            repmat({EventsTestTimePointColumn},size(overlap,1),1),repmat({EventsTargetTimePointColumn},size(overlap,1),1),...
            repmat({EventTargetTimeWindowOffsetTime},size(overlap,1),1),repmat({EventTargetTimeWindowPreOffsetTime},size(overlap,1),1),repmat({EventTargetTimeWindowPostOffsetTime},size(overlap,1),1),repmat({EventTestTimePointOffsetTime},size(overlap,1),1),'VarNames',varNames);
%     else
%         addRightOverlap = dataset(repmat({UseSecondColumnAndOnlyOffsetsForTimeWindow},size(overlap,1),1),repmat(tempEventsTestCompareColumns,size(overlap,1),1),...
%             repmat(tempEventsTargetCompareColumns,size(overlap,1),1),...
%             repmat({EventsTestTimePointColumn},size(overlap,1),1),repmat({EventsTargetTimePointColumn},size(overlap,1),1),...
%             {repmat(EventTargetTimeWindowOffsetTime,size(overlap,1),1)},{repmat(EventTargetTimeWindowPreOffsetTime,size(overlap,1),1)},{repmat(EventTargetTimeWindowPostOffsetTime,size(overlap,1),1)},{repmat(EventTestTimePointOffsetTime,size(overlap,1),1)},'VarNames',varNames);
%     end
    
%     if size(nonoverlap,1) > 1
        addRightNonOverlap = dataset(repmat({UseSecondColumnAndOnlyOffsetsForTimeWindow},size(nonoverlap,1),1),repmat(tempEventsTestCompareColumns,size(nonoverlap,1),1),...
            repmat(tempEventsTargetCompareColumns,size(nonoverlap,1),1),...
            repmat({EventsTestTimePointColumn},size(nonoverlap,1),1),repmat({EventsTargetTimePointColumn},size(nonoverlap,1),1),...
            repmat({EventTargetTimeWindowOffsetTime},size(nonoverlap,1),1),repmat({EventTargetTimeWindowPreOffsetTime},size(nonoverlap,1),1),repmat({EventTargetTimeWindowPostOffsetTime},size(nonoverlap,1),1),repmat({EventTestTimePointOffsetTime},size(nonoverlap,1),1),'VarNames',varNames);
%     else
%         addRightNonOverlap = dataset(repmat({UseSecondColumnAndOnlyOffsetsForTimeWindow},size(nonoverlap,1),1),repmat(tempEventsTestCompareColumns,size(nonoverlap,1),1),...
%             repmat(tempEventsTargetCompareColumns,size(nonoverlap,1),1),...
%             repmat(EventsTestTimePointColumn,size(nonoverlap,1),1),repmat(EventsTargetTimePointColumn,size(nonoverlap,1),1),...
%             {repmat(EventTargetTimeWindowOffsetTime,size(nonoverlap,1),1)},{repmat(EventTargetTimeWindowPreOffsetTime,size(nonoverlap,1),1)},{repmat(EventTargetTimeWindowPostOffsetTime,size(nonoverlap,1),1)},{repmat(EventTestTimePointOffsetTime,size(nonoverlap,1),1)},'VarNames',varNames);
%     end
    
    overlap = cat(2,overlap,addRightOverlap);
    nonoverlap = cat(2,nonoverlap,addRightNonOverlap);
    
    overlap = cat(2,dataset(repmat(iData,size(overlap,1),1),repmat({eventsTestPath},size(overlap,1),1),repmat({eventsTargetPath},size(overlap,1),1),'VarNames',{'event_files_number','test_event_file','target_event_file'}),overlap);
    nonoverlap = cat(2,dataset(repmat(iData,size(nonoverlap,1),1),repmat({eventsTestPath},size(nonoverlap,1),1),repmat({eventsTargetPath},size(nonoverlap,1),1),'VarNames',{'event_files_number','test_event_file','target_event_file'}),nonoverlap);
    
   
    
    tempEventsTestGroupSummaryByColumns = EventsTestGroupSummaryByColumns;
    
    if size(tempEventsTestGroupSummaryByColumns,2) > 1
        tempEventsTestGroupSummaryByColumns = {strjoin(EventsTestGroupSummaryByColumns,' ')};
    end
    
    groups = keys(groupByMapOverlap);
    if size(groups,2) < 2
            summary = dataset(repmat('group',2,1),'VarNames',{'group'});
            summary(2,:) = [];

    else
            summary = dataset(repmat('group',size(groups,2),1),'VarNames',{'group'});
    end
    for iGroup = 1:length(EventsTestGroupSummaryByColumns)
        tempGroupTest = EventsTestGroupSummaryByColumns{iGroup};
        summary = cat(2,summary,dataset(repmat(cellstr('group'),size(groups,2),1),'VarNames',{tempGroupTest}));
    end
    
    summary = cat(2,summary,dataset(values(groupByMapOverlap)','VarNames',{'number_match'}));
    summary = cat(2,summary,dataset(values(groupByMapNonOverlap)','VarNames',{'number_mismatch'}));
    summary = cat(2,summary,dataset(values(groupByMapAllTest)','VarNames',{'number_test_events'}));
    summary = cat(2,summary,dataset(values(groupByMapAllTarget)','VarNames',{'number_matching_target_events'}));
    summary = cat(2,summary,dataset(((cell2mat(summary.number_match) + cell2mat(summary.number_mismatch)) - cell2mat(summary.number_test_events)),'VarNames',{'number_test_matches_target_more_than_once'}));
    summary = cat(2,summary,dataset(repmat({tempEventsTestCompareColumns},size(summary,1),1),'VarNames',{'test_compare_columns'}));
    summary = cat(2,summary,dataset(repmat({tempEventsTargetCompareColumns},size(summary,1),1),'VarNames',{'target_compare_columns'}));
    summary = cat(2,summary,dataset(repmat({tempEventsTestGroupSummaryByColumns},size(summary,1),1),'VarNames',{'test_group_by_columns'}));
    summary = cat(2,summary,dataset(repmat({EventsTestFilterForColumn},size(summary,1),1),'VarNames',{'test_filter_column'}));
    summary = cat(2,summary,dataset(repmat({EventsTargetFilterForColumn},size(summary,1),1),'VarNames',{'target_filter_column'}));
    summary = cat(2,summary,dataset(repmat({EventsTestFilterValues},size(summary,1),1),'VarNames',{'test_filter_value'}));
    summary = cat(2,summary,dataset(repmat({EventsTargetFilterValues},size(summary,1),1),'VarNames',{'target_filter_value'}));

    for iGroupKey = 1:size(groups,2)
        groupsSplits = strsplit(groups{iGroupKey},GroupByConcatString);
        for iGroup = 1:length(EventsTestGroupSummaryByColumns)
            tempGroupTest = EventsTestGroupSummaryByColumns{iGroup};
            summary.(tempGroupTest)(iGroupKey) = groupsSplits(iGroup+1);
        end
    end
    
    if length(EventsTestGroupSummaryByColumns) > 0
        summary.group = [];
    end
    
    fprintf(['EventTestFile ' num2str(iData) ': write summary results\n']);
    export(summary,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_summary_' 'events_file_num_' num2str(iData) '.csv'],'Delimiter',',');
    
    fprintf(['EventTestFile ' num2str(iData) ': write match and non-match results\n']);
    export(overlap,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_matching_' 'events_file_num_' num2str(iData) '.csv'],'Delimiter',',');
    export(nonoverlap,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_mismatching_' 'events_file_num_' num2str(iData) '.csv'],'Delimiter',',');
    
    
    overlaps{conseciData} = overlap;
    nonoverlaps{conseciData} = nonoverlap;
    summarys{conseciData} = summary;

end


fprintf('Aggregate results of all datasets\n');
full_overlap = [];
full_nonoverlap = [];
full_summary = [];
for iData = conseciDatas
    %iData = 1
    if isempty(full_overlap)
        full_overlap = overlaps{iData};
    else
        if ~isempty(overlaps{iData})
            full_overlap = cat(1,full_overlap,overlaps{iData});
        end
    end
    
    if isempty(full_nonoverlap)
        full_nonoverlap = nonoverlaps{iData};
    else
        if ~isempty(nonoverlaps{iData})
            full_nonoverlap = cat(1,full_nonoverlap,nonoverlaps{iData});
        end

    end
    
    if isempty(full_summary)
        full_summary = summarys{iData};
    else
        full_summary = cat(1,full_summary,summarys{iData});
    end
end

export(full_overlap,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_matching_' 'all_recent' '.csv'],'Delimiter',',');
export(full_nonoverlap,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_mismatching_' 'all_recent' '.csv'],'Delimiter',',');
export(full_summary,'file',[pathOutputFolder filesep ouputFilesPrefixString 'events_cooccurrance_summary_' 'all_recent' '.csv'],'Delimiter',',');

res_match = full_overlap;
res_mismatch = full_nonoverlap;
res_summary = full_summary;


fprintf('EventCooccurrence function finished\n');
toc
memtoc
end