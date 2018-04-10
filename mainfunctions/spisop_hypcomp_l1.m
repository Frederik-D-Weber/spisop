function [res_stats_table_all, res_stats_table_single_category_all, res_stats_table_single_category_MA_all, res_cc_all_all_norm, res_cc_all_all_MA] = spisop_hypcomp_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters)
% compare sleep scorings against each other
% Copyright Frederik D. Weber

functionName = 'hypcomp';
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


HypnogramLinksFileName = getParam('HypnogramLinksFileName',listOfParameters);
listOfHypnogramLinksFileName = {};
listOfListsHypnogramLinksFileName = {};
if exist([pathInputFolder filesep HypnogramLinksFileName],'file') ~= 2
    error(['HypnogramLinksFileName file ' [pathInputFolder filesep HypnogramLinksFileName] ' does not exist. Check if this file is a correct parameter, if so then check for correct path and if file exists in it.'])
end
listOfHypnogramLinksFileName = read_mixed_csv([pathInputFolder filesep HypnogramLinksFileName],',');
%     if ~(all(size(listOfDatasetsPaths) == size(listOfHypnogramLinksFileName)))
%         error('files or number of Datasetspaths and LinearDeviationMontagePaths are invalid or do not aggree')
%     end

for iHypFiles = 1:size(listOfHypnogramLinksFileName,1)
    if exist([listOfHypnogramLinksFileName{iHypFiles}],'file') ~= 2
        error(['The hypnogram links file given in HypnogramLinksFileName for dataset number ' num2str(iHypFiles)  ' does not exist'])
    end
    
    try
        listOfListsHypnogramLinksFileName{iHypFiles} = read_mixed_csv([listOfHypnogramLinksFileName{iHypFiles}],',');
        for iHypFilesInside = 1:size(listOfListsHypnogramLinksFileName{iHypFiles},1)
            tempHyps = listOfListsHypnogramLinksFileName{iHypFiles};
            if exist([tempHyps{iHypFilesInside}],'file') ~= 2
                error(['The hypnogram file referenced by listed in HypnogramLinksFileName for hypnogram link number ' num2str(iHypFiles)  ' and within listed hypnogram on position ' num2str(iHypFilesInside) ' is not readable'])
            end
        end
    catch err
        error(['The hypnogram link file listed in HypnogramLinksFileName for hypnogram link number ' num2str(iHypFiles)  ' is not readable'])
    end
end



epochLength = str2num(getParam('epochLength',listOfCoreParameters)); % in seconds

DataSetsWhich = getParam('DataSetsWhich',listOfParameters);%Datasets to be processed either all or subset if subset then DataSetsNumbers is used for selection default all
DataSetsNumbers = str2num(getParam('DataSetsNumbers',listOfParameters));%The line numbers of the Datasets to be processed if DataSetsWich parameter is set to subset


SleepOnsetDefinition = getParam('SleepOnsetDefinition',listOfParameters);

iDatas = 1:(length(listOfHypnogramLinksFileName));

if strcmp(DataSetsWhich,'subset')
    if ~(ismember(min(DataSetsNumbers),iDatas) && ismember(max(DataSetsNumbers),iDatas))
        error('Parameter DataSetsNumbers contains numbers not matching to any line number, e.g. too less DataSetPaths in DataSetPathsFile!')
    end
    iDatas = DataSetsNumbers;
end

HypnogramTimeTicks = str2num(getParam('HypnogramTimeTicks',listOfParameters));


SkipStatistics = getParam('SkipStatistics',listOfParameters);%,no,choose if statistics output shall be skipped e.g. for figure creation only either yes or no default no


GenerateHypnogramFigures = getParam('GenerateHypnogramFigures',listOfParameters);%choose if for all hypnograms a Figure shall be created. either yes or no default no
GenerateHypnogramFiguresConsensus = getParam('GenerateHypnogramFiguresConsensus',listOfParameters);%choose if for all hypnograms a Figure shall be created. either yes or no default no
GenerateHypnogramFiguresFormat = getParam('GenerateHypnogramFiguresFormat',listOfParameters);%choose format dimensions (in inches 1in=2.54 cm) of hypnograms. either png or epsc or svg or tiff or pdf or bmp or fig default png
GenerateHypnogramFiguresUnit = getParam('GenerateHypnogramFiguresUnit',listOfParameters);% choose dimension unit (in inches 1in=2.54 cm) of hypnograms. either points or normalized or inches or centimeters or pixels default inches
GenerateHypnogramFiguresFormatWidth = str2num(getParam('GenerateHypnogramFiguresFormatWidth',listOfParameters));% choose format dimensions in inches, (1 in=2.54 cm) of hypnograms. default 20
GenerateHypnogramFiguresFormatHeight = str2num(getParam('GenerateHypnogramFiguresFormatHeight',listOfParameters));%format dimensions in inches (1 in=2.54 cm) of hypnograms. default 8
GenerateHypnogramFiguresFormatResolution = getParam('GenerateHypnogramFiguresFormatResolution',listOfParameters);% choose resolution in pixesl per inches (1 in=2.54 cm) of hypnograms. default 300
GenerateHypnogramFiguresFormatFontSize = str2num(getParam('GenerateHypnogramFiguresFormatFontSize',listOfParameters));% Font size in units stated in Parameter GenerateHypnogramFiguresUnit. default 0.1

referenceOption = getParam('referenceOption',listOfParameters);%which hypnogram should be chosen as a reference for comparison. either firstInList or consensushypnogram default firstInList
statisticsAlphaLevel = str2num(getParam('statisticsAlphaLevel',listOfParameters));%the alpha level for the statistics default 0.05

AlignToSleepOnset = 'no';
try
    AlignToSleepOnset = getParam('AlignToSleepOnset',listOfParameters);%choose if the hypnograms should be aligned (cut) to the sleep onset either yes or no default no
catch e
    
end

core_cfg = [];
core_cfg.feedback = getParam('ft_cfg_feedback',listOfCoreParameters);

dummySampleRate = 100;
lightsOffSample = 0;
epochLengthSamples = epochLength * dummySampleRate;

res_cc_all_all_norm = [];
res_cc_all_all_MA = [];
res_stats_table_all = [];
res_stats_table_single_category_all = [];
res_stats_table_single_category_MA_all = [];

tic
memtic
fprintf('HypComp function initialized\n');
for iData = iDatas
    %iData = 11
    
    datasetsPath = listOfHypnogramLinksFileName{iData};
    hypnogramListPaths = listOfListsHypnogramLinksFileName{iData};
    
    %lightsOffSample = listOfLightsOffs(iData);
    
    fprintf('dataset %i: process all hypnograms info\n',iData);
    
    hypnList = {};
    hypnVecNorm = [];
    hypnVecMA = [];
    
    for iHyp = 1:numel(hypnogramListPaths)
        
        hypnogramPath = hypnogramListPaths{iHyp};
        
        [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
        
        
        if strcmp(AlignToSleepOnset,'yes')
            
            
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
            
            
            
            hypn = hypn(onsetCandidate:end,:);
            hypnStages = hypnStages(onsetCandidate:end,:);
            hypnEpochsBeginsSamples = hypnEpochsBeginsSamples(onsetCandidate:end,:);
            
            
        end
        
        hypnList{iHyp} = hypn;
        
        
        if strcmp(GenerateHypnogramFigures,'yes')
            fprintf('dataset %i: generate hypnogram figure of hypfile %d\n',iData,iHyp);
            
            titleName = sprintf('Hypnogram_datasetnum_%d_hypfile_%d',iData,iHyp);
            
            generateTheHypnogramFigure(hypn,hypnStages,hypnEpochsBeginsSamples,lightsOffSample,SleepOnsetDefinition,epochLengthSamples,epochLength,...
                GenerateHypnogramFiguresFormat,GenerateHypnogramFiguresUnit,GenerateHypnogramFiguresFormatWidth,GenerateHypnogramFiguresFormatHeight,GenerateHypnogramFiguresFormatResolution,GenerateHypnogramFiguresFormatFontSize,...
                dummySampleRate,pathOutputFolder,ouputFilesPrefixString,titleName,HypnogramTimeTicks...
                )
            
        end
        
        if iHyp > 1
            nEpochs = size(hypnVecNorm,2);
            if size(hypn,1) > nEpochs
                missingEpochs =  size(hypn,1) - nEpochs;
                hypnVecNorm = [hypnVecNorm NaN(size(hypnVecNorm,1),missingEpochs)];
                hypnVecMA = [hypnVecMA NaN(size(hypnVecMA,1),missingEpochs)];
            end
            
            if size(hypn,1) < nEpochs
                missingEpochs =  nEpochs - size(hypn,1);
                hypn = [hypn;NaN(missingEpochs,2)];
            end
        end
        hypnVecNorm(iHyp,:) = hypn(:,1)';
        hypnVecMA(iHyp,:) = hypn(:,2)';
        
        
        
    end
    
    
    if ~strcmp(SkipStatistics,'yes')
        fprintf('dataset %i: calculate consensus and statistics\n',iData);
        
        
        hypnVecNorm(isnan(hypnVecNorm)) = -1;
        
        [k_norm,sek_norm,p_norm,z_norm,ci_norm,kj_norm,sekj_norm,zkj_norm,pkj_norm,consensusModusVoteMatrix_norm,consensusModusCountMatrix_norm,cross_comparisons_norm,chi2_cross_norm,p_cross_norm] = muliple_kappa(hypnVecNorm,statisticsAlphaLevel,referenceOption);
        
        hypnVecMA(isnan(hypnVecMA)) = -1;
        
        [k_MA,sek_MA,p_MA,z_MA,ci_MA,kj_MA,sekj_MA,zkj_MA,pkj_MA,consensusModusVoteMatrix_MA,consensusModusCountMatrix_MA,cross_comparisons_MA,chi2_cross_MA,p_cross_MA] = muliple_kappa(hypnVecMA,statisticsAlphaLevel,referenceOption);
        
        temp_nRows = 1;%numel(kj_norm);
        temp_stats = repmat([k_norm,k_MA,sek_norm,sek_MA,p_norm,p_MA,z_norm,z_MA,ci_norm(1),ci_norm(2),ci_MA(1),ci_MA(2)],temp_nRows,1);
        temp_stats_table = array2table(temp_stats,'VariableNames',{'kappa','kappa_MA','kappa_SEM','kappa_SEM_MA','kappa_pvalue','kappa_MA_pvalue','kappa_zvalue','kappa_MA_zvalue','kappa_CI_lower','kappa_CI_higher','kappa_MA_CI_lower','kappa_MA_CI_higher'});
        
        
        
        
        temp_datasetnumber_table = array2table(repmat(iData,temp_nRows,1));
        temp_hypfilenumber_table = cell2table(repmat(cellstr(datasetsPath),temp_nRows,1));
        temp_datasetnumber_table.Properties.VariableNames = {'datasetnum'};
        temp_hypfilenumber_table.Properties.VariableNames = {'hyplinkfilepath'};
        
        temp_stats_table_all = [temp_datasetnumber_table temp_hypfilenumber_table temp_stats_table];
        
        writetable(temp_stats_table_all,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_stats_datasetnum_' num2str(iData) '.csv'],'FileType','text','Delimiter',',');
        
        
        if isempty(res_stats_table_all)
            res_stats_table_all = temp_stats_table_all;
        else
            res_stats_table_all = [res_stats_table_all ; temp_stats_table_all];
        end
        
        
        temp_nRows = numel(kj_norm);
        temp_stats2_norm = [(1:temp_nRows)',kj_norm',repmat(sekj_norm,temp_nRows,1),pkj_norm',zkj_norm'];
        temp_stats_table2_norm = array2table(temp_stats2_norm,'VariableNames',{'single_category','kappa_single','kappa_single_SEM','kappa_single_pvalue','kappa_single_zvalue'});
        
        
        temp_datasetnumber_table = array2table(repmat(iData,temp_nRows,1));
        temp_hypfilenumber_table = cell2table(repmat(cellstr(datasetsPath),temp_nRows,1));
        temp_datasetnumber_table.Properties.VariableNames = {'datasetnum'};
        temp_hypfilenumber_table.Properties.VariableNames = {'hyplinkfilepath'};
        temp_stats_table2_norm = [temp_datasetnumber_table temp_hypfilenumber_table temp_stats_table2_norm];
        writetable(temp_stats_table2_norm,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_stats_category_datasetnum_' num2str(iData) '.csv'],'FileType','text','Delimiter',',');
        
        
        if isempty(res_stats_table_single_category_all)
            res_stats_table_single_category_all = temp_stats_table2_norm;
        else
            res_stats_table_single_category_all = [res_stats_table_single_category_all ; temp_stats_table2_norm];
        end
        
        temp_nRows = numel(kj_MA);
        temp_stats2_MA = [(1:temp_nRows)',kj_MA',repmat(sekj_MA,temp_nRows,1),pkj_MA',zkj_MA'];
        temp_stats_table2_MA = array2table(temp_stats2_MA,'VariableNames',{'single_category','kappa_single_MA','kappa_single_SEM_MA','kappa_single_MA_pvalue','kappa_single_MA_zvalue'});
        
        
        temp_datasetnumber_table = array2table(repmat(iData,temp_nRows,1));
        temp_hypfilenumber_table = cell2table(repmat(cellstr(datasetsPath),temp_nRows,1));
        temp_datasetnumber_table.Properties.VariableNames = {'datasetnum'};
        temp_hypfilenumber_table.Properties.VariableNames = {'hyplinkfilepath'};
        temp_stats_table2_MA = [temp_datasetnumber_table temp_hypfilenumber_table temp_stats_table2_MA];
        writetable(temp_stats_table2_MA,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_stats_category_MA_datasetnum_' num2str(iData) '.csv'],'FileType','text','Delimiter',',');
        
        
        
        if isempty(res_stats_table_single_category_MA_all)
            res_stats_table_single_category_MA_all = temp_stats_table2_MA;
        else
            res_stats_table_single_category_MA_all = [res_stats_table_single_category_MA_all ; temp_stats_table2_MA];
        end
        
        cc_all_norm = concatCrossComp(iData,cross_comparisons_norm,hypnogramListPaths,chi2_cross_norm,p_cross_norm);
        cc_all_MA = concatCrossComp(iData,cross_comparisons_MA,hypnogramListPaths,chi2_cross_MA,p_cross_MA);
        
        if isempty(res_cc_all_all_norm)
            res_cc_all_all_norm = cc_all_norm;
        else
            res_cc_all_all_norm = [res_cc_all_all_norm ; cc_all_norm];
        end
        
        if isempty(res_cc_all_all_MA)
            res_cc_all_all_MA = cc_all_MA;
        else
            res_cc_all_all_MA = [res_cc_all_all_MA ; cc_all_MA];
        end
        
        if ~isempty(cc_all_norm)
            writetable(cc_all_norm,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_cross_comp_datasetnum_' num2str(iData) '.csv'],'FileType','text','Delimiter',',');
        end
        if ~isempty(cc_all_MA)
            writetable(cc_all_MA,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_cross_comp_MA_datasetnum_' num2str(iData) '.csv'],'FileType','text','Delimiter',',');
        end
        
        hypn = [consensusModusVoteMatrix_norm' consensusModusVoteMatrix_MA'];
        [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = getHypnogramParameters(hypn,epochLengthSamples);
        
        titleName = sprintf('hypcomp_consensus_hypnogram_datasetnum_%d',iData);
        fid_export = fopen([ouputFilesPrefixString titleName '.txt'], 'wt');
        
        for iRow = 1:size(hypn,1)
            fprintf(fid_export, '%i\t%i\n', hypn(iRow,:));
        end
        fclose(fid_export);
        
        %if strcmp(GenerateHypnogramFigures,'yes') && strcmp(GenerateHypnogramFiguresConsensus,'yes')
        if strcmp(GenerateHypnogramFiguresConsensus,'yes')
   
            fprintf('dataset %i: generate consensus hypnogram figure\n',iData);
            
            titleName = sprintf('Hypnogram_consensus_datasetnum_%d',iData);
            
            generateTheHypnogramFigure(hypn,hypnStages,hypnEpochsBeginsSamples,lightsOffSample,SleepOnsetDefinition,epochLengthSamples,epochLength,...
                GenerateHypnogramFiguresFormat,GenerateHypnogramFiguresUnit,GenerateHypnogramFiguresFormatWidth,GenerateHypnogramFiguresFormatHeight,GenerateHypnogramFiguresFormatResolution,GenerateHypnogramFiguresFormatFontSize,...
                dummySampleRate,pathOutputFolder,ouputFilesPrefixString,titleName,HypnogramTimeTicks...
                )
            
        end
        
    end
    
    
    
    
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %
    %             if size(hypn,1) < nEpochs
    %                 missingEpochs = nEpochs - size(hypn,1);
    %                 hypn(end+1:end+missingEpochs,:) = [ones(1,missingEpochs,1)*-1 zeros(0,missingEpochs,1)];
    %             end
    %
    %         hypn = [ones(nEpochs,1)*-1 zeros(nEpochs,1)];
    %
    %
    %
    %     sampleFreq = preDownsampleFreq;
    %     epochLengthSamples = epochLength * preDownsampleFreq;
    %     [hypn hypnStages hypnEpochs hypnEpochsBeginsSamples hypnEpochsEndsSamples] = readInSleepHypnogram(hypnogramPath,epochLengthSamples);
    %
    %
    
    %
    %
    %
    %     sleepOnsetTime{iData} = (hypnEpochsBeginsSamples(onsetCandidate) - lightsOffSample)/sampleFreq;
    %
    %     %preOffsetCandidate = max(find(strcmp(hypnStages(:,1),'S1') | strcmp(hypnStages(:,3),'NonREM')));
    %     preOffsetCandidate = max(find(strcmp(hypnStages(:,1),'S1') | strcmp(hypnStages(:,3),'NonREM') | strcmp(hypnStages(:,3),'REM') | strcmp(hypnStages(:,3),'MT')));
    %
    %     S1ind = find(strcmp(hypnStages(:,1),'S1'));
    %     S2ind = find(strcmp(hypnStages(:,1),'S2'));
    %     SWSind = find(strcmp(hypnStages(:,2),'SWS'));
    %     S4ind = find(strcmp(hypnStages(:,1),'S4'));
    %     REMind = find(strcmp(hypnStages(:,1),'REM'));
    %     S1OnsetTime{iData} = (min(S1ind(S1ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    %     S2OnsetTime{iData} = (min(S2ind(S2ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    %     SWSonsetTime{iData} = (min(SWSind(SWSind >= onsetCandidate)) - onsetCandidate)*epochLength;
    %     S4onsetTime{iData} = (min(S4ind(S4ind >= onsetCandidate)) - onsetCandidate)*epochLength;
    %     REMonsetTime{iData} = (min(REMind(REMind >= onsetCandidate)) - onsetCandidate)*epochLength;
    %
    %     preOnsetCandidate = onsetCandidate;
    %     if preOnsetCandidate > 1
    %         preOnsetCandidate = preOnsetCandidate-1;
    %     end
    %     hypnTST = hypn(onsetCandidate:preOffsetCandidate,:);
    %     hypnStagesTST = hypnStages(onsetCandidate:preOffsetCandidate,:);
    %     hypnStagesPreSleepOnset = hypnStages(1:preOnsetCandidate,:);
    %     hypnEpochsTST = hypnEpochs(onsetCandidate:preOffsetCandidate);
    %     hypnEpochsBeginsSamplesTST = hypnEpochsBeginsSamples(onsetCandidate:preOffsetCandidate,:);
    %     hypnEpochsEndsSamplesTST = hypnEpochsEndsSamples(onsetCandidate:preOffsetCandidate,:);
    %
    %     totalSleepTime{iData} = (length(onsetCandidate:preOffsetCandidate))*epochLength;
    %
    %     S1Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S1')))*epochLength;
    %     S2Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S2')))*epochLength;
    %     S3Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S3')))*epochLength;
    %     S4Time{iData} = length(find(strcmp(hypnStagesTST(:,1),'S4')))*epochLength;
    %     REMtime{iData} = length(find(strcmp(hypnStagesTST(:,1),'REM')))*epochLength;
    %     WakeTime{iData} = length(find(strcmp(hypnStagesTST(:,1),'Wake')))*epochLength;
    %     MovementTime{iData} = length(find(strcmp(hypnStagesTST(:,1),'MT')))*epochLength;
    %     SWStime{iData} = S3Time{iData} + S4Time{iData};
    %     NonREMtime{iData} = SWStime{iData} + S2Time{iData};
    %
    %
    %     S1TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S1')))*epochLength;
    %     S2TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S2')))*epochLength;
    %     S3TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S3')))*epochLength;
    %     S4TimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'S4')))*epochLength;
    %     REMtimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'REM')))*epochLength;
    %     WakeTimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'Wake')))*epochLength;
    %     MovementTimePreOnset{iData} = length(find(strcmp(hypnStagesPreSleepOnset(:,1),'MT')))*epochLength;
    %     SWStimePreOnset{iData} = S3TimePreOnset{iData} + S4TimePreOnset{iData};
    %     NonREMtimePreOnset{iData} = SWStimePreOnset{iData} + S2TimePreOnset{iData};
    %
    %     S1Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S1') & (hypnTST(:,2) == 0) ))*epochLength;
    %     S2Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S2') & (hypnTST(:,2) == 0) ))*epochLength;
    %     S3Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S3') & (hypnTST(:,2) == 0) ))*epochLength;
    %     S4Time_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'S4') & (hypnTST(:,2) == 0) ))*epochLength;
    %     REMtime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'REM') & (hypnTST(:,2) == 0) ))*epochLength;
    %     WakeTime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'Wake') & (hypnTST(:,2) == 0) ))*epochLength;
    %     MovementTime_WithoutMA{iData} = length(find(strcmp(hypnStagesTST(:,1),'MT') & (hypnTST(:,2) == 0) ))*epochLength;
    %     SWStime_WithoutMA{iData} = S3Time_WithoutMA{iData} + S4Time_WithoutMA{iData};
    %     NonREMtime_WithoutMA{iData} = SWStime_WithoutMA{iData} + S2Time_WithoutMA{iData};
    %
    %
    %     if strcmp(ExportHypnogram,'yes')
    %         tempExpochFactor = 60/epochLength;
    %         tempExportPostfix = '_';
    %         if strcmp(ExportHypnogramStartOffset,'sleeponset')
    %             tempExportPostfix = [tempExportPostfix 'sleeponset'];
    %             tempExportHypnogramStart = onsetCandidate;
    %         elseif strcmp(ExportHypnogramStartOffset,'eegonset')
    %              tempExportPostfix = [tempExportPostfix 'eegonset'];
    %              tempExportHypnogramStart = 1;
    %         elseif strcmp(ExportHypnogramStartOffset,'lightsoff')
    %              tempExportPostfix = [tempExportPostfix 'lightsoff'];
    %              tempExportHypnogramStart = find(hypnEpochsEndsSamples <= lightsOffSample,1,'last')-1;
    %
    %         else
    %             error(['wrong parameter for ExportHypnogramStartOffset either sleeponset or eegonset but set is ' ExportHypnogramStartOffset] );
    %         end
    %        tempExportPostfix = [tempExportPostfix '_' num2str(ExportHypnogramStartTime) '_to_' num2str(ExportHypnogramEndTime) '_min'];
    %
    %
    %
    %        tempExportHypnogramStartTimeEpoch = tempExportHypnogramStart + round(ExportHypnogramStartTime*tempExpochFactor);
    %        tempExportHypnogramEndTimeEpoch = tempExportHypnogramStart + round(ExportHypnogramEndTime*tempExpochFactor);
    %        tempExport_hypn = hypn;
    %
    %        if ((tempExportHypnogramStartTimeEpoch-1) > 0)
    %            tempExport_hypn(1:(tempExportHypnogramStartTimeEpoch-1),2) = 3;
    %        end
    %
    %        if (tempExportHypnogramEndTimeEpoch <= size(tempExport_hypn,1))
    %            tempExport_hypn((tempExportHypnogramEndTimeEpoch):end,2) = 3;
    %        end
    %
    %
    %        [temp_pathstr,temp_name,temp_ext] = fileparts(hypnogramPath);
    %        fid_export = fopen([ouputFilesPrefixString num2str(iData) '_' temp_name tempExportPostfix '.txt'], 'wt');
    %
    %        for iRow = 1:size(tempExport_hypn,1)
    %            fprintf(fid_export, '%i %i\n', tempExport_hypn(iRow,:));
    %        end
    %        fclose(fid_export);
    %     end
    %
    %
    %
    %
    %
    %
    
end
% %open output files
%     fidh = fopen([pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '.csv'],'wt');
%
%     %write header of ouptufiles
%     fprintf(fidh,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n','datasetnum','dataset','hypnogram','epoch_length_seconds','Total_sleep_time_min','Sleep_Onset_min','S1_onset_min','S2_onset_min','SWS_onset_min','S4_onset_min','REM_onset_min'...
%         ,'S1_min','S2_min','S3_min','S4_min','REM_min','Wake_after_sleep_onset_min','Movement_Time_min','SWS_min','NonREM_without_S1_min'...
%         ,'S1_percent','S2_percent','S3_percent','S4_percent','REM_percent','Wake_after_sleep_onset_percent','Movement_Time_percent','SWS_percent','NonREM_without_S1_percent'...
%         ,'S1_without_MA_min','S2_without_MA_min','S3_without_MA_min','S4_without_MA_min','REM_without_MA_min','Wake_after_sleep_onset_without_MA_min','Movement_Time_without_MA_min','SWS_without_MA_min','NonREM_without_S1_without_MA_min'...
%         ,'S1_without_MA_percent','S2_without_MA_percent','S3_without_MA_percent','S4_without_MA_percent','REM_without_MA_percent','Wake_after_sleep_onset_without_MA_percent','Movement_Time_without_MA_percent','SWS_without_MA_percent','NonREM_without_S1_without_MA_percent'...
%         ,'S1_before_sleep_onset_min','S2_before_sleep_onset_min','S3_before_sleep_onset_min','S4_before_sleep_onset_min','REM_before_sleep_onset_min','Wake_before_sleep_onset_min','Movement_before_sleep_onset_Time_min','SWS_before_sleep_onset_min','NonREM_before_sleep_onset_without_S1_min');
%
%     for iData = iDatas
%         fprintf(fidh,'%i,',iData);
%         fprintf(fidh,'%s,',listOfDatasetsPaths{iData});
%         fprintf(fidh,'%s,',listOfHypnogramPaths{iData});
%         fprintf(fidh,'%f,',epochLength);
%         fprintf(fidh,'%f,',totalSleepTime{iData}/60);
%         fprintf(fidh,'%f,',sleepOnsetTime{iData}/60);
%         fprintf(fidh,'%f,',S1OnsetTime{iData}/60);
%         fprintf(fidh,'%f,',S2OnsetTime{iData}/60);
%         fprintf(fidh,'%f,',SWSonsetTime{iData}/60);
%         fprintf(fidh,'%f,',S4onsetTime{iData}/60);
%         fprintf(fidh,'%f,',REMonsetTime{iData}/60);
%
%         fprintf(fidh,'%f,',S1Time{iData}/60);
%         fprintf(fidh,'%f,',S2Time{iData}/60);
%         fprintf(fidh,'%f,',S3Time{iData}/60);
%         fprintf(fidh,'%f,',S4Time{iData}/60);
%         fprintf(fidh,'%f,',REMtime{iData}/60);
%         fprintf(fidh,'%f,',WakeTime{iData}/60);
%         fprintf(fidh,'%f,',MovementTime{iData}/60);
%         fprintf(fidh,'%f,',SWStime{iData}/60);
%         fprintf(fidh,'%f,',NonREMtime{iData}/60);
%         fprintf(fidh,'%f,',100*S1Time{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S2Time{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S3Time{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S4Time{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*REMtime{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*WakeTime{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*MovementTime{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*SWStime{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*NonREMtime{iData}/totalSleepTime{iData});
%
%         fprintf(fidh,'%f,',S1Time_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',S2Time_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',S3Time_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',S4Time_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',REMtime_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',WakeTime_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',MovementTime_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',SWStime_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',NonREMtime_WithoutMA{iData}/60);
%         fprintf(fidh,'%f,',100*S1Time_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S2Time_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S3Time_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*S4Time_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*REMtime_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*WakeTime_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*MovementTime_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*SWStime_WithoutMA{iData}/totalSleepTime{iData});
%         fprintf(fidh,'%f,',100*NonREMtime_WithoutMA{iData}/totalSleepTime{iData});
%
%         fprintf(fidh,'%f,',S1TimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',S2TimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',S3TimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',S4TimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',REMtimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',WakeTimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',MovementTimePreOnset{iData}/60);
%         fprintf(fidh,'%f,',SWStimePreOnset{iData}/60);
%         fprintf(fidh,'%f\n',NonREMtimePreOnset{iData}/60);
%
%     end
%     fclose(fidh);
%     res_hypnvals = dataset('File',[pathOutputFolder filesep ouputFilesPrefixString 'hypvals_full_' 'datanum_all_selected' '.csv'],'Delimiter',',');

if ~isempty(res_cc_all_all_norm)
    writetable(res_cc_all_all_norm,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_cross_comp_datasetnum_' 'all_recent' '.csv'],'FileType','text','Delimiter',',');
end
if ~isempty(res_cc_all_all_MA)
    writetable(res_cc_all_all_MA,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_cross_comp_MA_datasetnum_' 'all_recent' '.csv'],'FileType','text','Delimiter',',');
end
if ~isempty(res_stats_table_all)
    writetable(res_stats_table_all,[pathOutputFolder filesep ouputFilesPrefixString 'hypcomp_stats_datasetnum_' 'all_recent' '.csv'],'FileType','text','Delimiter',',');
end

fprintf('HypComp function finished\n');

toc
memtoc
end

function onsetCandidateIndex = getSleepOnsetEpoch(hypnStages,hypnEpochsBeginsSamples,lightsOffSample,SleepOnsetDefinition)
% Sleep onsset: [S1] S1 S2, but not [S1] S1 X S2, where X is not NonREM
% S1, starting with first S1 is sleep onset, nothing else
% S1_NonREM, starting with first S1 followed directly by either S2, S3 or S4,
%            otherwise with first S2 or S3 or S4
% S1_XREM, starting with first S1 directly followed by either S2, S3, S4, or REM,
%            otherwise with first S2 or S3 or S4 or REM
% NonREM, starting with first one of S2, S3 or S4
% XREM, starting with first one of S2, S3 or S4 or REM

onsetCandidateIndex = -1;

if strcmp(SleepOnsetDefinition,'NonREM')
    for iOnset = 1:(size(hypnStages,1))
        if strcmp(hypnStages(iOnset,3),'NonREM') && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
            onsetCandidateIndex = iOnset;
            break;
        end
    end
    
elseif strcmp(SleepOnsetDefinition,'XREM')
    for iOnset = 1:(size(hypnStages,1))
        if (strcmp(hypnStages(iOnset,3),'NonREM') ||  strcmp(hypnStages(iOnset,3),'REM')) && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
            onsetCandidateIndex = iOnset;
            break;
        end
    end
    
elseif strcmp(SleepOnsetDefinition,'S1') || strcmp(SleepOnsetDefinition,'S1_NonREM') || strcmp(SleepOnsetDefinition,'S1_XREM')
    
    consecS1 = 0;
    hasS1 = logical(0);
    for iOnset = 1:(size(hypnStages,1))
        if strcmp(hypnStages(iOnset,1),'S1') && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
            hasS1 = logical(1);
            consecS1 = consecS1 + 1;
            if ((onsetCandidateIndex+consecS1) ~= iOnset)
                onsetCandidateIndex = iOnset;
                consecS1 = 0;
            end
            if strcmp(SleepOnsetDefinition,'S1')
                break;
            end
        elseif ( strcmp(SleepOnsetDefinition,'S1_XREM') && (strcmp(hypnStages(iOnset,3),'NonREM') || strcmp(hypnStages(iOnset,3),'REM')) ) ...
                || ( strcmp(SleepOnsetDefinition,'S1_NonREM') && (strcmp(hypnStages(iOnset,3),'NonREM')) ) ...
                && (hypnEpochsBeginsSamples(iOnset) >= lightsOffSample)
            if ~hasS1
                onsetCandidateIndex = iOnset;
            end
            break;
        else
            consecS1 = 0;
            hasS1 = logical(0);
        end
    end
    
end
end




function generateTheHypnogramFigure(hypn,hypnStages,hypnEpochsBeginsSamples,lightsOffSample,SleepOnsetDefinition,epochLengthSamples,epochLength,...
    GenerateHypnogramFiguresFormat,GenerateHypnogramFiguresUnit,GenerateHypnogramFiguresFormatWidth,GenerateHypnogramFiguresFormatHeight,GenerateHypnogramFiguresFormatResolution,GenerateHypnogramFiguresFormatFontSize,...
    dummySampleRate,pathOutputFolder,ouputFilesPrefixString,titleName,HypnogramTimeTicks...
    )


nEpochs = size(hypn,1);
onsetCandidateIndex = getSleepOnsetEpoch(hypnStages,hypnEpochsBeginsSamples,lightsOffSample,SleepOnsetDefinition);
preOffsetCandidate = max(find(strcmp(hypnStages(:,1),'S1') | strcmp(hypnStages(:,3),'NonREM') | strcmp(hypnStages(:,3),'REM') | strcmp(hypnStages(:,3),'MT')));
if isempty(preOffsetCandidate)
    preOffsetCandidate = nEpochs;
end



%%% plot hypnogram figure
plot_MA_offset = -5.5;
[hypn_plot_interpol hypn_plot_interpol_MA] = interpolate_hypn_for_plot(hypn,epochLengthSamples,plot_MA_offset);

hhyp = figure;
hhypfigax = gca;
set(hhyp,'color',[1 1 1]);
set(hhypfigax,'FontUnits',GenerateHypnogramFiguresUnit)
set(hhypfigax,'Fontsize',GenerateHypnogramFiguresFormatFontSize);
temp_max_y = 1.0;

x_time = (1:length(hypn_plot_interpol))/(dummySampleRate);
x_time = x_time/60; % minutes

x_time_hyp = x_time(1:length(hypn_plot_interpol));

axh = hhypfigax;
plot(axh,x_time_hyp,hypn_plot_interpol,'Color',[0 0 0])
hold(axh,'on');
if onsetCandidateIndex ~= -1
    onset_time = (onsetCandidateIndex-0.5)*(epochLength/60);%in minutes
    onset_y_coord_offset = 0.2;
    onset_y_coord = hypn_plot_interpol(find(x_time >=onset_time,1,'first'))+onset_y_coord_offset;
    hold(axh,'on');
    scatter(axh,onset_time,onset_y_coord,'filled','v','MarkerFaceColor',[0 1 0])
end


offset_time = (preOffsetCandidate+0.5)*(epochLength/60);%in minutes
offset_y_coord_offset = 0.2;
offset_y_coord = hypn_plot_interpol(find(x_time <=offset_time,1,'last'))+offset_y_coord_offset;
hold(axh,'on');
scatter(axh,offset_time,offset_y_coord,'filled','^','MarkerFaceColor',[0 0 1])

plot(axh,x_time_hyp,hypn_plot_interpol_MA,'Color',[1 0 0])
xlim(axh,[0 max(x_time)]);
ylabel(axh,'Stages');
ylim(axh,[plot_MA_offset temp_max_y])
yTick = [1 0 -0.5 -1 -2 -3 -4 plot_MA_offset+1 plot_MA_offset+0.5];
yTickLabel = {'?' 'W' 'REM' 'S1' 'S2' 'S3' 'S4' 'MT' 'MA'};
set(axh, 'yTick', flip(yTick));
set(axh, 'yTickLabel', flip(yTickLabel));
set(axh,'TickDir','out');
xTick = [0:HypnogramTimeTicks:max(x_time)];
set(axh, 'xTick', xTick);
set(axh, 'box', 'off')

%     begsample = 0;
%     endsample = 0;
%     x_pos_begin = x_time(begsample);
%     x_pos_end = x_time(endsample);
%     x_pos = [x_pos_begin x_pos_end x_pos_end x_pos_begin];
%     y_pos = [plot_MA_offset plot_MA_offset 1 1];
%     pos_now = patch(x_pos,y_pos,[0.5 0.25 1],'parent',axh);
%     set(pos_now,'FaceAlpha',0.4);
%     set(pos_now,'EdgeColor','none');

%     line([x_pos_begin x_pos_begin],[plot_MA_offset temp_max_y],'color',[0.25 0.125 1],'parent',axh);

%titleName = sprintf('Hypnogram_datasetnum_%d_file_%d',iData,iHyp);
set(hhyp, 'Name', titleName);

hold(axh,'off')

title(titleName,'Interpreter','none');
xlabel('Time [min]');
ylabel('Sleep stage');


figure_width = GenerateHypnogramFiguresFormatWidth;     % Width in inches
figure_height = GenerateHypnogramFiguresFormatHeight;    % Height in inches
pos = get(hhyp, 'Position');

%set(hhyp, 'Position', [pos(1) pos(2) figure_width*str2num(GenerateHypnogramFiguresFormatResolution), figure_height*str2num(GenerateHypnogramFiguresFormatResolution)]); %<- Set size
set(hhyp, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
% Here we preserve the size of the image when we save it.
set(hhyp,'InvertHardcopy','on');
set(hhyp,'PaperUnits', GenerateHypnogramFiguresUnit);

%set(hhyp,'PaperPositionMode','Auto')
set(hhyp,'PaperSize',[figure_width, figure_height])

papersize = get(hhyp, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(hhyp,'PaperPosition', myfiguresize);
set(hhyp,'PaperOrientation', 'portrait');


switch GenerateHypnogramFiguresFormat
    case 'fig'
        saveas(hhyp, [pathOutputFolder filesep ouputFilesPrefixString titleName '.fig']);
    case 'eps'
        print(hhyp,['-d' 'epsc'],['-r' GenerateHypnogramFiguresFormatResolution],[pathOutputFolder filesep ouputFilesPrefixString titleName]);
    otherwise
        print(hhyp,['-d' GenerateHypnogramFiguresFormat],['-r' GenerateHypnogramFiguresFormatResolution],[pathOutputFolder filesep ouputFilesPrefixString titleName]);
end
%saveas(hhyp, [titleName '.png']);
%saveas(hhyp, [titleName '.eps']);

close(hhyp);


%%% plot hypnogram figure end

end

function cc_all = concatCrossComp(iData,cross_comparisons,hypnogramListPaths,chi2_cross,p_cross)
cc_all = [];
for iCrossComp = 1:numel(cross_comparisons)
    %iCrossComp = 2
    if isempty(cross_comparisons{iCrossComp})
        continue;
    end
    temp_cc_nrows = size(cross_comparisons{iCrossComp},1);
    
    temp_datasetnumber_table = array2table(repmat(iData,temp_cc_nrows,1));
    temp_hypfilenumber_table = array2table(repmat(iCrossComp,temp_cc_nrows,1));
    temp_hyps_table = cell2table(repmat(cellstr(hypnogramListPaths{iCrossComp}),temp_cc_nrows,1));
    temp_stage_code_table = cell2table(cross_comparisons{iCrossComp}.Properties.RowNames);
    temp_chi2_table = array2table(repmat(chi2_cross(iCrossComp),temp_cc_nrows,1));
    temp_chi2_pvalue_table = array2table(repmat(p_cross(iCrossComp),temp_cc_nrows,1));
    
    
    temp_datasetnumber_table.Properties.VariableNames = {'datasetnum'};
    temp_hypfilenumber_table.Properties.VariableNames = {'hypfilenum'};
    temp_hyps_table.Properties.VariableNames = {'hypfilepath'};
    temp_stage_code_table.Properties.VariableNames = {'stage_code_test_hyps'};
    temp_chi2_table.Properties.VariableNames = {'test_chi2'};
    temp_chi2_pvalue_table.Properties.VariableNames = {'test_chi2_pvalue'};
    
    
    temp_cc = [temp_datasetnumber_table temp_hypfilenumber_table temp_hyps_table temp_chi2_table temp_chi2_pvalue_table temp_stage_code_table cross_comparisons{iCrossComp}];
    
    temp_cc.Properties.RowNames =  cellstr(arrayfun(@(x) [num2str(iData) '_' num2str(iCrossComp) '_' x], num2str((1:size(temp_cc,1))'),'UniformOutput',false));
    
    if isempty(cc_all)
        cc_all = temp_cc;
    else
        cc_all = [cc_all ; temp_cc];
    end
end
end
