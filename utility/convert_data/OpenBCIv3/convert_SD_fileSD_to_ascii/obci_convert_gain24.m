function [ out ] = obci_convert_gain24(format,currentFullInstallationFilePath,fileFilterSettings,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NStandardParam = 3;
NMaxParam = 5;

% convert_channelwise_brainvision_int16(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage)

DoReReference = 'no';
fileRerefSettings = 'dummy';
ApplyLinearDeviationMontage = 'no';
linearDeviationMontagePath = 'dummy';

if (nargin < NStandardParam) || (nargin > NMaxParam)
    fprintf(['Wrong number of input arguments, should be between ' num2str(NStandardParam) ' and ' num2str(NMaxParam) ' but is ' nargin '\n'])
    printUsage();
end

if nargin > NStandardParam 
    for iArg = 1:length(varargin)
        currArg = lower(varargin{iArg});
        currArgPair = strsplit(currArg,'=');
        if (length(currArgPair) ~= 2)
            error(['the argmument number ' num2str(iArg) ' was not a parameter=value pair, please rewrite it']);
        end
        if (strcmp(currArgPair(2),''))
            error(['the argmument number ' num2str(iArg) ' has an empty value (parameter=VALUE), please rewrite it']);
        end
        switch lower(currArgPair{1})
            case 'filererefsettings'
                DoReReference = 'yes';
                fileRerefSettings = currArgPair{2};
            case 'lineardeviationmontagepath'
                ApplyLinearDeviationMontage = 'yes';
                linearDeviationMontagePath = currArgPair{2};
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

switch lower(format) 
    case 'brainvision_int16'
        convert_channelwise_brainvision_int16(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);
    case 'brainvision_int32'
        convert_channelwise_brainvision_int32(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);
    case 'brainvision_float32'
        convert_channelwise_brainvision_float32(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);
    case 'edf_autoscale'
        convert_channelwise_edf_autoscale(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);
    case 'edf_0p1uvacc_cutoff'
        convert_channelwise_edf_0p1uVacc_cutoff(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);
    case 'edf_0p01uvacc_cutoff'
        convert_channelwise_edf_0p01uVacc_cutoff(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);  
    case 'ascii'
        convert_channelwise_ascii(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage);  
    otherwise
        fprintf(['format with name ' format ' is unknown'])
        printUsage();
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
%    finalize_parallel_computing()
%end
out = 'finished';

end

function printUsage()
sprintf(['USAGE: obci_convert_gain24(.exe) format currentFullInstallationFilePath fileFilterSettings\n'...
        'REQUIRED PARAMETERS: ' '\n'...
        format,currentFullInstallationFilePath,fileFilterSettings
        ' format is the name of the format to convert to after filtering e.g. brainvision_float32 ...see below' '\n'...
        ' currentFullInstallationFilePath is the path to the toolbox directory e.g. D:\\spisop_toolbox_beta2.3' '\n'...
        ' fileFilterSettings is the full name of the file containing the filter settings located at the folder named settings' '\n'...
        ' FORMATs are:' '\n'...
        '  brainvision_int16' '\n'...
        '  brainvision_int32' '\n'...
        '  brainvision_float32' '\n'...
        '  edf_autoscale' '\n'...
        '  edf_0p1uVacc_cutoff' '\n'...
        '  edf_0p01uVacc_cutoff' '\n'...
        '  ascii' '\n'...
        'OPTIONAL PARAMETERS as \"[parameter]=[value]\" pairs: ' '\n'...
        ' fileRerefSettings=[filename] is the full name of the file containing the rereference settings located at the folder named settings' '\n'...
        ' linearDeviationMontagePath=[filename] is the full name of the file containing the linear deviation matrix located at the folder named settings' '\n'...
        ]);
    error('try again!')
end