function [ out ] = spisop(functionName,currentFullInstallationFilePath,inputFolderName,outputFolderName,coreParamFile,functionParamFile,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NStandardParam = 6;
NMaxParam = 8;

prefixExtentionName = '';
parallelComputation = 'no';
numParallelWorkers = 2;



if (nargin < NStandardParam) || (nargin > NMaxParam)
    if (nargin < 1)

    hlp = helpdlg(sprintf(['SpiSOP is a command line tool!\n\n It does NOT run by double clicking it.\nYou need to give parameters in a command line... \n\n' ...
        getUsageString()]),'It does not work like this!');
    end
    uiwait(hlp);

    fprintf(['Wrong number of input arguments, should be between ' num2str(NStandardParam) ' and ' num2str(NMaxParam) ' but is ' num2str(nargin) '\n'])
    printUsage();
    if (~isdeployed)
        error('try again!');
    else
        exit;
    end
end

if nargin > 6 && ~(strcmpi(functionName,'manipulateparam')) || (nargin > 5 && (strcmpi(functionName,'eventcooccurrence')))
    for iArg = 1:length(varargin)
        currArg = lower(varargin{iArg});
        currArgPair = strsplit(currArg,'=');
        if (length(currArgPair) ~= 2)
            error(['the argmument number ' num2str(iArg) ' was not a parameter=value pair, please rewrite it']);
        end
        if (strcmp(currArgPair(2),''))
            error(['the argmument number ' num2str(iArg) ' has an empty value (parameter=VALUE), please rewrite it']);
        end
        switch currArgPair{1}
            case 'parallelcomputation'
                numParallelWorkers = str2num(currArgPair{2});
                if numParallelWorkers > 0
                    parallelComputation = 'yes';
                else
                    parallelComputation = 'no';
                end
            case 'prefixextentionname'
                prefixExtentionName = currArgPair{2};
            otherwise
               error(['the argmument valuepair' currArg 'is unknown!']);
        end
    end
end

%[currentFullInstallationFilePath,dummy_name,temp_ext]  = fileparts( mfilename('fullpath'));

% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!

%% Use debuging after error occured if the next line is uncommented
% dbstop if error

%% Paths to toolbox, input and output folder

% Give the complete path to the toolbox folder (depends on your operating
% system (Mac/Unix/Linux or Windows)
pathPrefix = currentFullInstallationFilePath;% 
%% Specify the dataset naming in the
% datasetAndInputFolderName, it gives a name you want for the datasets/study/etc 
% that you want to  analyze
% !!! There needs to be another folder in the input folder with the same
% name !!!
% an output folder with the same name is created automatically if not existent
%datasetAndInputFolderName = 'test_lindev';

%% Parallelization !requires Parallel Computing Toolbox (distcomp), otherwise there are errors!
% here one can choose if parallel computing is enabled (yes) or not (no)
% in case of enabled the number of workeres are specified
%parallelComputation = 'no'; % either 'yes' or 'no' default 'yes'
% if isdeployed
%     numParallelWorkers = str2num(numParallelWorkers);% number of workers, use at maximum one less than the number of real cpu cores in your computer default 1
% end
%% the actual initialization happens here
if (~isdeployed)
    addpath(pathPrefix);
end
%spisop_init_helper(datasetAndInputFolderName,datasetAndInputFolderName,parallelComputation,numParallelWorkers);
[pathInputFolder, pathOutputFolder] = spisop_init_helper(pathPrefix,inputFolderName,outputFolderName,parallelComputation,numParallelWorkers);

try
    if (~isdeployed)
        %...\external_code\automatic_sleep_scoring\z3score
        %Download cfslib-MATLAB from https://github.com/amiyapatanaik/cfslib-MATLAB
        addpath([pathPrefix  filesep 'external_code' filesep 'automatic_sleep_scoring' filesep 'z3score' filesep 'z3score-api' filesep 'cfslib-MATLAB']);
        addpath([pathPrefix  filesep 'external_code' filesep 'automatic_sleep_scoring' filesep 'z3score' filesep 'z3score-api' filesep 'cfslib-MATLAB' filesep 'utilities']);
        addpath([pathPrefix  filesep 'external_code' filesep 'automatic_sleep_scoring' filesep 'z3score' filesep 'z3score-api' filesep 'cfslib-MATLAB' filesep 'utilities' filesep 'jsonlab']);
        addpath([pathPrefix  filesep 'external_code' filesep 'automatic_sleep_scoring' filesep 'z3score' filesep 'z3score-api' filesep 'cfslib-MATLAB' filesep 'utilities' filesep 'encoder']);
    end
    
    if (~isdeployed)
        %...\external_code\REMs_detector\marek_adamczyk\REMdetector
        addpath([pathPrefix  filesep 'external_code' filesep 'REMs_detector' filesep 'marek_adamczyk' filesep 'REMdetector']);
    end
    
catch
    warning('Z3Score: could not be initialized');
end

%% -------- run until here to initilize -------- % 


%% set prefix sting/text
%the test analyisis for dataset 'short' with three short night datasets:
outputFilesPrefixString = [outputFolderName '_' prefixExtentionName '_'];


switch lower(functionName) 
    case 'hypvals'
        spisop_hypvals(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'preprocout'
        spisop_preprocout(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'browser'
        spisop_browser(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'sod'
        spisop_sod(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'spd'
        spisop_spd(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'pow'
        spisop_pow(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'freqpeaks'
        spisop_freqpeaks(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'confounds'
        spisop_confounds(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'remsmaad'
        spisop_remsMaAd(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'nonevents'
        spisop_nonEvents(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile);
    case 'hypcomp'
        spisop_hypcomp(pathInputFolder, pathOutputFolder, outputFilesPrefixString, coreParamFile, functionParamFile)
    case 'eventcooccurrence' 
        fprintf('ignoring the core paremater file (coreParamFile) since it is not needed d\n');
        spisop_eventCooccurrence(pathInputFolder, pathOutputFolder, outputFilesPrefixString, functionParamFile);
    case 'manipulateparam' 
        fprintf(['ignoring the output folder name (outputFolderName) and stay in input folder with the name(' inputFolderName ') since output folder is not needed\n']);
        inParamFilename = coreParamFile;
        outParamFilename = functionParamFile;
        param = varargin{1};
        value_str = varargin{2};
        %cd(pathInputFolder);
        spisop_manipulateParam(pathInputFolder, inParamFilename, outParamFilename, param, value_str);
    otherwise
        fprintf(['function with name ' functionName ' is unknown'])
        printUsage();
        if (~isdeployed)
            error('try again!');
        else
            exit;
        end
end


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

%% finalization of parallel computing (if enabled)
%at the end of computation or analysis please finalize the parallel computing
%if strcmp(parallelComputation,'yes')
    finalize_parallel_computing()
%end
out = 'finished';

end

function fstring = getUsageString()
fstring = ['USAGE: spisop(.exe) functionName currentFullInstallationFilePath inputFolderName outputFolderName coreParamFile functionParamFile\n'...
        'REQUIRED PARAMETERS: ' '\n'...
        ' functionName is the name of the function to call e.g. hypvals spd sod ...see below' '\n'...
        ' currentFullInstallationFilePath is the path to the toolbox directory e.g. D:\\spisop_toolbox_beta2.3' '\n'...
        ' inputFolderName is a name to the input folder with the parameter files in them' '\n'...
        ' outputFolderName is a name of the folder (to be created) in the general \"output\" folder to place the output files' '\n'...
        ' coreParamFile is the full name of the file containing the core parameters stored in the input folder with inputFolderName' '\n'...
        ' functionParamFile is the full name of the file containing the function parameters stored in the input folder with inputFolderName' '\n'...
        ' FUNCTIONS for functionName are:' '\n'...
        '  hypvals' '\n'...
        '  preprocout' '\n'...
        '  browser' '\n'...
        '  sod' '\n'...
        '  pow' '\n'...
        '  freqpeaks' '\n'...
        '  confounds' '\n'...
        '  remsMaAd' '\n'...
        '  nonEvents' '\n'...
        '  hypcomp' '\n'...
        '  eventCooccurrence' '\n'...
        '  manipulateparam' '\n'...
        '    here the parameters are: spisop(.exe) manipulateparam currentFullInstallationFilePath inputFolderName outputFolderName inParamFilename outParamFilename param value_str\n'...   
        '    where outputFolderName will be ingnored, inParamFilename is paramter file name and outParamFilename the parameter file to write out to, param is the parameter name to manimpulate, value_str gives the value to overwrite the old paramter value\n'...   
        '' '\n'...
        'OPTIONAL PARAMETERS as \"[parameter]=[value]\" pairs: ' '\n'...
        ' prefixExtentionName=[value] give a prefix extention string attached to outputfiles to the standard prefix' '\n'...
        ' parallelComputation=[integer >0 and < 13] give number of workers for parallel processing to initiate' '\n'...
        ];
end

function printUsage()
usagestring = getUsageString();
    disp(sprintf([ '\n' usagestring ]));
end