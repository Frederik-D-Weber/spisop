% Initiation of SpiSOP (abrev. for Spindles Slow Oscillations and Power) toolbox
% Copyright Frederik D. Weber, see README.txt and COPYING.txt for more information
% if you do not agree with the licencing and use stated there, you are not
% allowed to use this Software!
clear
pathPrefix = 'D:\spisop_toolbox_beta2.3';

%% settings
%fileFilterSettings = 'split_filter_home_jingyi_setup1.txt';
%fileFilterSettings = 'split_reref_filter_home_jingyi_setup1.txt';
%fileFilterSettings = 'split_reref_lindev_filter_home_jingyi_setup1.txt';
%fileFilterSettings = 'split_reref_lindev_filter_home_jingyi_setup1_spind_slow.txt';
%fileFilterSettings = 'split_pure_lindev_filter_home_jingyi_setup1_spind_slow.txt';


%fileRerefSettings = 'reref.txt';
%fileRerefSettings = 'reref_mult.txt';

%DoReReference = 'no';

%linearDeviationMontagePath = 'linear_deviation_sleep_scoring.txt';
%ApplyLinearDeviationMontage = 'yes';


%%

serverPathPrefix = [pathPrefix filesep 'utility/convert_data/OpenBCIv3/convert_SD_file/SD_to_ascii_then_to_EEG'];
folderPath = [serverPathPrefix ''];
settingsFolderPath =  [serverPathPrefix filesep 'settings'];
cd(folderPath);

if (~isdeployed)
    addpath([pathPrefix filesep 'fieldtrip_fw']);
    ft_defaults;
end

%obci_convert_gain24('edf_autoscale',pathPrefix,'split_pure_lindev_filter_home_jingyi_setup1_spind_slow.txt','linearDeviationMontagePath=linear_deviation_sleep_scoring.txt')
%obci_convert_gain24('ascii',pathPrefix,'split_pure_lindev_filter_home_jingyi_setup1_spind_slow.txt','linearDeviationMontagePath=linear_deviation_sleep_scoring.txt')

cd(serverPathPrefix)
mcc -mv obci_convert_gain24.m