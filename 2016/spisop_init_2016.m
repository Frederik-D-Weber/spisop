% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!

%% Use debuging after error occured if the next line is uncommented
% dbstop if error

%% Paths to toolbox, input and output folder

% Give the complete path to the toolbox folder (depends on your operating
% system (Mac/Unix/Linux or Windows)
pathPrefix = 'D:\spisop_toolbox_beta2.3';% 
%% Specify the dataset naming in the
% datasetAndInputFolderName, it gives a name you want for the datasets/study/etc 
% that you want to  analyze
% !!! There needs to be another folder in the input folder with the same
% name !!!
% an output folder with the same name is created automatically if not existent
%datasetAndInputFolderName = 'test_lindev';
datasetAndInputFolderName = 'test_short';
%datasetAndInputFolderName = 'PS_128';

%% Parallelization !requires Parallel Computing Toolbox (distcomp), otherwise there are errors!
% here one can choose if parallel computing is enabled (yes) or not (no)
% in case of enabled the number of workeres are specified
parallelComputation = 'no'; % either 'yes' or 'no' default 'yes'
numParallelWorkers = 2;% number of workers, use at maximum one less than the number of real cpu cores in your computer default 1

%% the actual initialization happens herer
addpath(pathPrefix);
[pathInputFolder, pathOutputFolder] = spisop_init_helper2016(pathPrefix,datasetAndInputFolderName,datasetAndInputFolderName,parallelComputation,numParallelWorkers);

%% -------- run until here to initilize -------- % 





%% analysis 'PS_128' of Annika Hanert without artifact epochs
outputFilesPrefixString = [datasetAndInputFolderName '_' 'run2' '_'];
tempCoreParamFile = 'CoreParam_PS128.txt';

spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SOParam_PS128.txt');



%% analysis 'test_lindev' without artifact epochs

%% set prefix sting/text
%the test analyisis for dataset 'short' with three short night datasets:
outputFilesPrefixString = [datasetAndInputFolderName '_' 'run1' '_'];
tempCoreParamFile = 'CoreParam_lindev.txt';
%tempCoreParamFile2 = 'CoreParam_sleep_scoring_lindev.txt';
tempCoreParamFile2 = 'CoreParam_sleep_scoring_simple_lindev.txt';

%% Hypnogram values for sleep table
spisop_hypvals(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'HypValsParam_lindev.txt');
spisop_preprocout(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'PreProcOutParam_lindev.txt');
spisop_preprocout(pathInputFolder, pathOutputFolder, [outputFilesPrefixString 'sleep_score_'], tempCoreParamFile2, 'PreProcOutParam_lindev.txt');
spisop_browser(pathInputFolder, pathOutputFolder, [outputFilesPrefixString 'sleep_score_'], tempCoreParamFile2, 'BrowserParam_lindev.txt');
spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SOParam_lindev.txt');
spisop_spd(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SpindlesParam_lindev.txt');

spisop_hypcomp(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'HypCompParam_lindev.txt');


%%
%spisop('browser',pathPrefix,parallelComputation,numParallelWorkers,datasetAndInputFolderName,datasetAndInputFolderName,'testme','CoreParam_sleep_scoring_simple_lindev.txt','BrowserParam_lindev.txt')
%spisop.exe hypvals D:\spisop_toolbox_beta2.3 test_lindev test_lindev_out CoreParam_sleep_scoring_simple_lindev.txt HypValsParam_lindev.txt
%spisop.exe browser D:\spisop_toolbox_beta2.3 test_lindev test_lindev_out CoreParam_sleep_scoring_simple_lindev.txt BrowserParam_lindev.txt
%spisop(pathPrefix,datasetAndInputFolderName,datasetAndInputFolderName,'CoreParam_sleep_scoring_simple_lindev.txt','browser','BrowserParam_lindev.txt','prefixExtentionName=run1','parallelComputation=1')

%spisop('manipulateParam',pathPrefix,datasetAndInputFolderName,datasetAndInputFolderName,'BrowserParam_lindev_test.txt','BrowserParam_lindev_test_out.txt','DataSetsNumbers','1 2 3 4:5')


%% test analysis 'test_short'

%% set prefix sting/text
%the test analyisis for dataset 'short' with three short night datasets:
outputFilesPrefixString = [datasetAndInputFolderName '_run1_'];
tempCoreParamFile = 'CoreParam_short.txt';

%% Hypnogram values for sleep table
spisop_hypvals(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'HypValsParam_short.txt');
%% Hypnogram export function (expample first 20 min after sleep onset)
spisop_hypvals(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'HypValsParam_first20min_after_sleeponset_short.txt');
%% Power (density) estimation
spisop_pow(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'PowParam_short.txt');
%% Frequency peaks for sleep spindles (find slow and fast spindles power peaks in spectra)
spisop_freqpeaks(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SpindlesFreqPeaksParam_short.txt');
%% sleep spindle detection based on frequency band defined by center frequencies
spisop_spd(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SpindlesParam_short.txt');
%% slow oscillation detection 
spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'SOParam_short.txt');
%% confound artifact detection, still in trial phase, not recommended yet for use
%spisop_confounds(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_short.txt', 'ConfoundsParam_short.txt');
spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_short.txt', 'SOParam_short.txt');

%% scoring browser 
tempCoreParamFile_sleep_score = 'CoreParam_sleep_scoring_simple_short.txt';
spisop_browser(pathInputFolder, pathOutputFolder, [outputFilesPrefixString 'sleep_score_'], tempCoreParamFile_sleep_score, 'BrowserParam_short.txt');


outputFilesPrefixString = [datasetAndInputFolderName '_run1_'];
tempCoreParamFile = 'CoreParam_basel.txt';
spisop_hypcomp(pathInputFolder, pathOutputFolder, outputFilesPrefixString, tempCoreParamFile, 'HypCompParam_basel.txt');



%% test analysis 'test_EMSA'

%the test analyisis for dataset 'EMSA' with two full night datasets:
outputFilesPrefixString = [datasetAndInputFolderName '_'];
%% Hypnogram values for sleep table
spisop_hypvals(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'HypValsParam_EMSA.txt');
%% Power (density) estimation
spisop_pow(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'PowParam_EMSA.txt');
%% Frequency peaks for sleep spindles (find slow and fast spindles power peaks in spectra)
spisop_freqpeaks(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'SpindlesFreqPeaksParam_EMSA.txt');
%% sleep spindle detection based on frequency band defined by center frequencies
spisop_spd(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'SpindlesParam_EMSA.txt');
%% slow oscillation detection 
spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'SOParam_EMSA.txt');
%% Frequency peaks for slow oscillations(find slow oscillatoin power peaks in spectra)
spisop_freqpeaks(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'SOFreqPeaksParam_EMSA.txt');
%% find equivalent non-events (for spindles)
spisop_nonEvents(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'CoreParam_EMSA.txt', 'NonEventsParam_EMSA.txt');
%% find co-occurrences of events (spindles and slow oscillations)
%profile on
spisop_eventCooccurrence(pathInputFolder, pathOutputFolder, outputFilesPrefixString, 'EventsCooccurranceParam.txt');
%profile viewer
%% test analysis 'EMSA' with integrated functions


%example for integrated use within matlab
outputFilesPrefixString = 'test_integrated_';
coreParameterFileName = 'CoreParam_EMSA.txt';
temp_listOfCoreParameters = read_mixed_csv([pathInputFolder filesep coreParameterFileName],',');
parameterFileName = 'SpindlesParam_EMSA.txt';
listOfParameters = read_mixed_csv([pathInputFolder filesep parameterFileName],',');
temp_channelsOfInterestFileName = 'ChannelsOfInterest_EMSA.txt';
channelsOfInterest = read_mixed_csv([pathInputFolder filesep temp_channelsOfInterestFileName],',');
temp_channelsOfInterestFileName = ['temp_' temp_channelsOfInterestFileName];
csvwrite([pathInputFolder filesep temp_channelsOfInterestFileName],channelsOfInterest);

% for factorSD within 1 to 2 in steps of 0.1 do a spindle detection 
for iFactor = 1:0.1:2 
    temp_str_factorSD = num2str(iFactor);
    temp_ouputFilesPrefixString = [outputFilesPrefixString 'factorSD_' temp_str_factorSD '_']; 
    temp_listOfParameters = listOfParameters;
    temp_listOfParameters = setParam('factorSD',temp_listOfParameters,temp_str_factorSD);
    %only use a subset of datasets
    temp_listOfParameters = setParam('DataSetsWhich',temp_listOfParameters,'subset');
    temp_str_DataSetsNumbers = num2str([1 2]);
    temp_listOfParameters = setParam('DataSetsNumbers',temp_listOfParameters,temp_str_DataSetsNumbers);
    [res_filters, res_channels, res_events, res_peaks, res_troughs] = spisop_spd_l1(pathInputFolder, pathOutputFolder, temp_ouputFilesPrefixString, temp_listOfCoreParameters, temp_listOfParameters);
    %... here you can check the result datasets that have also been writen
    %to the output folder
end

%% finalization of parallel computing (if enabled)
%at the end of computation or analysis please finalize the parallel computing
finalize_parallel_computing()


