% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!

%% Use debuging after error occured if the next line is uncommented
% dbstop if error
clear
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
datasetAndInputFolderName = 'test_short';

%% Parallelization !requires Parallel Computing Toolbox (distcomp), otherwise there are errors!
% here one can choose if parallel computing is enabled (yes) or not (no)
% in case of enabled the number of workeres are specified
parallelComputation = 'no'; % either 'yes' or 'no' default 'yes'
numParallelWorkers = 2;% number of workers, use at maximum one less than the number of real cpu cores in your computer default 1

%% the actual initialization happens here
addpath(pathPrefix);
[pathInputFolder, pathOutputFolder] = spisop_init_helper2016(pathPrefix,datasetAndInputFolderName,datasetAndInputFolderName,parallelComputation,numParallelWorkers);

cd(pathPrefix)
mcc -mv -I ./subfunctions -I ./mainfunctions -a ./subfunctions/*.m -a ./dummy/dummy.eeg -a ./dummy/dummy.vhdr -a ./dummy/dummy.txt -a ./dummy/dummy.vmrk  spisop2016.m