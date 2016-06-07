% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!

%% Use debuging after error occured if the next line is uncommented
% dbstop if error

%% Paths to toolbox, input and output folder
% Give the complete path to the toolbox folder (depends on your operating
% system (Mac/Unix/Linux or Windows)
currentFullInstallationFilePath = 'D:\spisop_toolbox_beta2.3';% 

%% Parallelization !requires Parallel Computing Toolbox (distcomp), otherwise there are errors!
% here one can choose if parallel computing is enabled (yes) or not (no)
% in case of enabled the number of workeres are specified
numParallelWorkers = '1';% number of workers, use at maximum one less than the number of real cpu cores in your computer default 1

%% the actual initialization happens herer
addpath(currentFullInstallationFilePath);

%% -------- run until here to initilize -------- % 


%% tutorial analysis 'test_short'

% Specify the name of your dataset below (datasetAndInputFolderName) and
% create a folder in the /input/ folder with the same name. This is
% essential! The input folder will contain text files that provide all
% important parameters and settings.
% An output folder with the same name will be created automatically if not 
% existent. These are the parameters for the tutorial dataset 'test_short' 
% with three short night datasets:
datasetAndInputFolderName = 'test_short';
outputFilesPrefixString = 'run1';
tempCoreParamFile = 'CoreParam_short.txt';

%% Hypnogram values for sleep table
spisop('hypvals',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'HypValsParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Hypnogram export function (expample first 20 min after sleep onset)
spisop('hypvals',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'HypValsParam_first20min_after_sleeponset_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Power (density) estimation
spisop('pow',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'PowParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Frequency peaks for sleep spindles (find slow and fast spindles power peaks in spectra)
spisop('freqpeaks',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SpindlesFreqPeaksParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Sleep spindle detection based on frequency band defined by center frequencies
spisop('spd',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SpindlesParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Slow oscillation detection 
spisop('sod',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SOParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Confound artifact detection, still in trial phase, not recommended yet for use
spisop('confounds',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'ConfoundsParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% scoring browser 
tempCoreParamFile_sleep_score = 'CoreParam_sleep_scoring_simple_short.txt';
spisop('browser',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile_sleep_score,'BrowserParam_short.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])


%% test analysis 'test_EMSA'

%% set prefix sting/text
%the test analyisis for dataset 'short' with three short night datasets:
datasetAndInputFolderName = 'test_EMSA';
outputFilesPrefixString = 'run1';
tempCoreParamFile = 'CoreParam_EMSA.txt';

%% Hypnogram values for sleep table
spisop('hypvals',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'HypValsParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Power (density) estimation
spisop('pow',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'PowParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Frequency peaks for sleep spindles (find slow and fast spindles power peaks in spectra)
spisop('freqpeaks',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SpindlesFreqPeaksParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Sleep spindle detection based on frequency band defined by center frequencies
spisop('spd',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SpindlesParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Slow oscillation detection 
spisop('sod',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SOParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Frequency peaks for slow oscillations(find slow oscillatoin power peaks in spectra)
spisop('freqpeaks',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SOFreqPeaksParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% find equivalent non-events (for spindles)
spisop('nonEvents',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'NonEventsParam_EMSA.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% find co-occurrences of events (spindles and slow oscillations)
spisop('eventCooccurrence',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'EventsCooccurranceParam.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])


%% test analysis 'EMSA' with integrated functions
% %example for integrated use within matlab
% outputFilesPrefixString = 'test_integrated_';
% coreParameterFileName = 'CoreParam_EMSA.txt';
% temp_listOfCoreParameters = read_mixed_csv([pathInputFolder filesep coreParameterFileName],',');
% parameterFileName = 'SpindlesParam_EMSA.txt';
% listOfParameters = read_mixed_csv([pathInputFolder filesep parameterFileName],',');
% temp_channelsOfInterestFileName = 'ChannelsOfInterest_EMSA.txt';
% channelsOfInterest = read_mixed_csv([pathInputFolder filesep temp_channelsOfInterestFileName],',');
% temp_channelsOfInterestFileName = ['temp_' temp_channelsOfInterestFileName];
% csvwrite([pathInputFolder filesep temp_channelsOfInterestFileName],channelsOfInterest);
% 
% % for factorSD within 1 to 2 in steps of 0.1 do a spindle detection 
% for iFactor = 1:0.1:2 
%     temp_str_factorSD = num2str(iFactor);
%     temp_ouputFilesPrefixString = [outputFilesPrefixString 'factorSD_' temp_str_factorSD '_']; 
%     temp_listOfParameters = listOfParameters;
%     temp_listOfParameters = setParam('factorSD',temp_listOfParameters,temp_str_factorSD);
%     %only use a subset of datasets
%     temp_listOfParameters = setParam('DataSetsWhich',temp_listOfParameters,'subset');
%     temp_str_DataSetsNumbers = num2str([1 2]);
%     temp_listOfParameters = setParam('DataSetsNumbers',temp_listOfParameters,temp_str_DataSetsNumbers);
%     [res_filters, res_channels, res_events, res_peaks, res_troughs] = spisop_spd_l1(pathInputFolder, pathOutputFolder, temp_ouputFilesPrefixString, temp_listOfCoreParameters, temp_listOfParameters);
%     %... here you can check the result datasets that have also been writen
%     %to the output folder
% end

%% analysis 'test_lindev' without artifact epochs

%% set prefix sting/text
%the test analyisis for dataset 'short' with three short night datasets:
datasetAndInputFolderName = 'test_lindev';
outputFilesPrefixString = 'run1';
tempCoreParamFile = 'CoreParam_lindev.txt';
tempCoreParamFile2 = 'CoreParam_sleep_scoring_simple_lindev.txt';

%% Hypnogram values for sleep table
spisop('hypvals',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'HypValsParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Preprocess the data and save in another format, applying linear deviation and filtering
spisop('preprocout',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'PreProcOutParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Preprocess the data and save in another format, applying linear deviation and filtering (for sleep scoring outside of toolbox)
spisop('preprocout',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile2,'PreProcOutParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString 'sleep_score_'],['parallelComputation=' numParallelWorkers])
%% Slow oscillation detection 
spisop('sod',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SOParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% Sleep spindle detection based on frequency band defined by center frequencies
spisop('spd',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'SpindlesParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% compare hypnograms
spisop('hypcomp',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'HypCompParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])
%% sleep scoring with browser from the raw data 
spisop('browser',currentFullInstallationFilePath,datasetAndInputFolderName,datasetAndInputFolderName,tempCoreParamFile,'BrowserParam_lindev.txt',['prefixExtentionName=' outputFilesPrefixString],['parallelComputation=' numParallelWorkers])


%% finalization of parallel computing (if enabled)
%at the end of computation or analysis please finalize the parallel computing
finalize_parallel_computing()


