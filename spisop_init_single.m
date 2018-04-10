% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!

%% Use debuging after error occured if the next line is uncommented
dbstop if error

%% Paths to toolbox, input and output folder

% Give the complete path to the toolbox folder (depends on your operating
% system (Mac/Unix/Linux or Windows)
pathPrefix = 'D:\spisop_toolbox_beta2.3';% 

datasetAndInputFolderName = 'davos'
%% Parallelization !requires Parallel Computing Toolbox (distcomp), otherwise there are errors!
% here one can choose if parallel computing is enabled (yes) or not (no)
% in case of enabled the number of workeres are specified
parallelComputation = 'no'; % either 'yes' or 'no' default 'yes'
numParallelWorkers = 2;% number of workers, use at maximum one less than the number of real cpu cores in your computer default 1

%% the actual initialization happens herer
addpath(pathPrefix);
[pathInputFolder, pathOutputFolder] = spisop_init_helper(pathPrefix,datasetAndInputFolderName,datasetAndInputFolderName,parallelComputation,numParallelWorkers);

%% -------- run until here to initilize -------- % 


outputFilesPrefixString = [datasetAndInputFolderName '_' 'run1' '_'];

spisop_browser('', pathOutputFolder, [outputFilesPrefixString 'sleep_score_test'], [pathPrefix filesep 'standard' filesep 'single_exec' filesep 'CoreParam.txt'], [pathPrefix filesep 'standard' filesep 'single_exec' filesep 'BrowserParam.txt']);


%% finalization of parallel computing (if enabled)
%at the end of computation or analysis please finalize the parallel computing
finalize_parallel_computing()


