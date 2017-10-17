function convert_channelwise_edf_0p1uVacc_cutoff(currentFullInstallationFilePath,fileFilterSettings,fileRerefSettings,DoReReference,linearDeviationMontagePath,ApplyLinearDeviationMontage)
%dbstop if error
%clear

%load working direktory

%% Paths to toolbox, input and output folder

% Give the complete path to the toolbox folder (depends on your operating
% system (Mac/Unix/Linux or Windows)
pathPrefix = currentFullInstallationFilePath;% 

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

% load fieldtrip
%addpath([serverPathPrefix 'fieldtrip-20150506']);
if (~isdeployed)
    addpath([pathPrefix filesep 'fieldtrip_fw']);
    ft_defaults;
end


DelimiterLinearDeviationMontage = ',';

%sd_dataFilesNames = dir('*.TXT');
%sd_dataFilesNames = {sd_dataFilesNames.name}';


if (ispc)
    system('convert_gain24.bat');
elseif (isunix)
    system('bash convert_gain24.sh');
end

dataFilesNames = dir('*.TXT.new.csv');
dataFilesNames = {dataFilesNames.name}';

channel_settings_table = readtable([settingsFolderPath filesep fileFilterSettings],'Delimiter',',');
splitFactors = unique(channel_settings_table.split_factor);

for iData = 1:(length(dataFilesNames))
    % iData = 1;
    fprintf('dataset %i: read raw data\n',iData);
    dataFileNamePath = [folderPath filesep dataFilesNames{iData}];
    table = readtable(dataFileNamePath,'FileType','text','Delimiter',';');
    table = table2dataset(table);
    %table = dataset('file',dataFileNamePath,'Delimiter',';');
    table.valid = (mod(((table.index-1)-table.ringpar),256) == 0);
    
    all(table.valid)
    
    channelColumns = {};
    file_chanelColumns = 4:11;
    
    for iSplitFactor = 1:length(splitFactors)
        %iSplitFactor = 1;
        
        fprintf('dataset %i, split %i: process raw data\n',iData,iSplitFactor);
        curr_splitFactorsRows = find(strcmp(channel_settings_table.split_factor,splitFactors(iSplitFactor)));
        curr_channel_settings_table_analysis = channel_settings_table(curr_splitFactorsRows,:);
        curr_splitFactorsRows_raw_index = curr_channel_settings_table_analysis.channel_number > 0;
        curr_channel_settings_table_raw = curr_channel_settings_table_analysis(curr_splitFactorsRows_raw_index,:);
        curr_chanOrder_raw = curr_channel_settings_table_raw.channel_order;
        [dummy_chanOrder curr_chanIndexOrder_raw] = sort(curr_chanOrder_raw);
        curr_channel_settings_table_raw = curr_channel_settings_table_raw(curr_chanIndexOrder_raw,:);
        
        curr_channelColumns_raw = file_chanelColumns(curr_channel_settings_table_raw.channel_number(curr_splitFactorsRows(curr_splitFactorsRows_raw_index)));
        curr_splitFactorsNewChannelLabels_raw = curr_channel_settings_table_raw.new_channel_label;
                
        %mV in
        %if ~isa(table(2,curr_channelColumns(1)),'double')
        dat = double(table(:,curr_channelColumns_raw));%EEG.Data;%
        %else
        %    dat = table(:,curr_channelColumns);%EEG.Data;%
        %end
        % ! potential is in milliVolt(mV) so divide by 1000 for microV,
        % therefore no adjustments
        dat = dat';
        
        %dat = dat(curr_chanIndexOrder_raw,:);
        
        %varNames = get(table,'VarNames');
        %chanNames = varNames(curr_channelColumns)';
        
        curr_splitFactors_signal_correction_factor_raw = curr_channel_settings_table_raw.signal_correction_factor;
        for iChannel = 1:size(dat,1)
            dat(iChannel,:) = dat(iChannel,:) * curr_splitFactors_signal_correction_factor_raw(iChannel);
        end
        
        hdr.Fs = 250;%EEG.SamplingRate;
        hdr.nChans = length(curr_splitFactorsNewChannelLabels_raw);%EEG.ChannelNumber;
        hdr.nSamples = size(dat,2);%EEG.Points;
        hdr.nSamplesPre = 0;
        hdr.nTrials = 1;
        hdr.label = curr_splitFactorsNewChannelLabels_raw;%EEG.ChannelTitles;
        hdr.chantype = repmat({'eeg'},hdr.nChans,1);%Nx1 cell-array with the channel type, see FT_CHANTYPE
        hdr.chanunit = repmat({'uV'},hdr.nChans,1);%EEG.ChannelUnits;%Nx1 cell-array with the physical units, see FT_CHANUNIT
        
        dat(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(dat,1),1);
        
        % ! potential is assumed to be in microV (µV) not milliVolt(mV)
        dataformat = 'brainvision_eeg';
        hdr.brainvision_outformat = 'float32';%float32 int16 int32;
        
        
        data_file_name_unfiltered = [dataFilesNames{iData} '_' char(splitFactors(iSplitFactor)) '_' hdr.brainvision_outformat];
        ft_write_data([data_file_name_unfiltered '.eeg'], dat,'dataformat',dataformat,'header',hdr);
        
        dat = [];

        
        
        data_file_name_current = data_file_name_unfiltered;
        
        
        cfg = [];
        cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
        
        cfg.dataset = [data_file_name_current '.eeg'];
        cfg.headerfile = [data_file_name_current '.vhdr'];
        cfg.channel = 'all';
        data = ft_preprocessing(cfg);
        hdr = data.hdr;
        if strcmp(DoReReference,'yes') || strcmp(ApplyLinearDeviationMontage,'yes')
            
            if strcmp(DoReReference,'yes')
                
                
                table_reref = readtable([settingsFolderPath filesep fileRerefSettings],'FileType','text','Delimiter',',');
                for iReref = size(table_reref,1)
                    %iReref = 1;
                    %referenceChannels = {'A1', 'A2'};
                %toBeReferencedChannels = {'C*', 'F*'};
                %newImplicitRefChannelLabel = 'Cz';
                referenceChannels = strsplit(table_reref.referenceChannels{iReref});
                toBeReferencedChannels = strsplit(table_reref.toBeReferencedChannels{iReref});
                
                referenceChannels = referenceChannels(~(cellfun(@isempty,referenceChannels)));
                toBeReferencedChannels = toBeReferencedChannels(~(cellfun(@isempty,toBeReferencedChannels)));
                
                newImplicitRefChannelLabel = table_reref.newImplicitRefChannelLabel{iReref};
                
                
                fprintf('dataset %i, split %i: rereferencing data\n',iData,iSplitFactor);
                cfg = [];
                cfg.reref       = DoReReference;
                cfg.channel = ft_channelselection([referenceChannels toBeReferencedChannels],data.label);
                %cfg.channel = cellstr([referenceChannels toBeReferencedChannels]');
                cfg.implicitref = newImplicitRefChannelLabel;% the implicit (non-recorded) reference channel is added to the data representation
                cfg.refchannel     = referenceChannels;
                data_reref        = ft_preprocessing(cfg,data);
                
                cfg = [];
                notChan = strcat('-',data_reref.label);
                cfg.channel = ['all'; notChan(:)];
                data = ft_selectdata(cfg, data);

                cfg = [];
                data = ft_appenddata(cfg, data, data_reref);
                data_reref = [];
                
                end
%                 reorderIndex = [];
%                 for iLabel = 1:length(data.label)
%                     %iLabel = 2
%                     if (iLabel <= length(curr_channel_settings_table_raw.new_channel_label))
%                         %curr_channel_settings_table.new_channel_label(iLabel)
%                         reorderIndex =  [reorderIndex; find(strcmp(curr_channel_settings_table_raw.new_channel_label(iLabel),data.label),1,'first')];
%                     else
%                        reorderIndex =  [reorderIndex; iLabel]; 
%                     end
%                 end
                curr_channel_settings_table_reref = curr_channel_settings_table_analysis(ismember(curr_channel_settings_table_analysis.new_channel_label,data.label),:);
                if ~(size(curr_channel_settings_table_reref,1) == length(data.label))
                        error('not all channels emerging/created after rerefenced have been covered in the settings file for filtering')
                end
                
                [dummy_chan curr_indexOrder_setting] = ismember(data.label,curr_channel_settings_table_reref.new_channel_label);
                curr_chanOrder_reref = curr_channel_settings_table_reref.channel_order(curr_indexOrder_setting);
                [dummy_chanOrder curr_chanIndexOrder_reref] = sort(curr_chanOrder_reref);

                datreoder = data.trial{1};
                datreoder = datreoder(curr_chanIndexOrder_reref,:);
                data.trial{1} = datreoder;
                data.label = data.label(curr_chanIndexOrder_reref);
                
                
                hdr.Fs = data.fsample;%EEG.SamplingRate;
                hdr.nChans = length(data.label);%EEG.ChannelNumber;
                hdr.nTrials = 1;
                hdr.label = data.label;%EEG.ChannelTitles;
                hdr.chantype = repmat({'eeg'},hdr.nChans,1);%Nx1 cell-array with the channel type, see FT_CHANTYPE
                hdr.chanunit = repmat({'uV'},hdr.nChans,1);%EEG.ChannelUnits;%Nx1 cell-array with the physical units, see FT_CHANUNIT
                
                
                
                sigpositive_data = data.trial{1};
                sigpositive_data(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(sigpositive_data,1),1);
                data.trial{1} = sigpositive_data;
                
                hdr.brainvision_outformat = 'float32';%float32 int16 int32;
                data_file_name_current = [data_file_name_current '_' 'reref'];
                ft_write_data([data_file_name_current '.eeg'], data.trial{:},'dataformat',dataformat,'header',hdr);
                
            end
            
            if strcmp(ApplyLinearDeviationMontage,'yes')
                fprintf('dataset %i, split %i: apply linear deviation montage to data\n',iData,iSplitFactor);
                
                
                montageTable = dataset('File',[settingsFolderPath filesep linearDeviationMontagePath],'Delimiter',DelimiterLinearDeviationMontage,'ReadVarNames',true,'ReadObsNames',true);
                
                montage.labelorg = get(montageTable,'VarNames');
                montage.labelnew  = get(montageTable,'ObsNames')';
                montage.tra = double(montageTable);
                
                data = ft_apply_montage(data,montage,'keepunused','yes','inverse','no');
                
                curr_channel_settings_table_lindev = curr_channel_settings_table_analysis(ismember(curr_channel_settings_table_analysis.new_channel_label,data.label),:);
                if ~(size(curr_channel_settings_table_lindev,1) == length(data.label))
                        error('not all channels emerging/created after applying linear deviations have been covered in the settings file for filtering')
                end
                
     
                [dummy_chan curr_indexOrder_setting] = ismember(data.label,curr_channel_settings_table_lindev.new_channel_label);
                curr_chanOrder_lindev = curr_channel_settings_table_lindev.channel_order(curr_indexOrder_setting);
                [dummy_chanOrder curr_chanIndexOrder_lindev] = sort(curr_chanOrder_lindev);


                datreoder = data.trial{1};
                datreoder = datreoder(curr_chanIndexOrder_lindev,:);
                data.trial{1} = datreoder;
                data.label = data.label(curr_chanIndexOrder_lindev);
                
                
                hdr.Fs = data.fsample;%EEG.SamplingRate;
                hdr.nChans = length(data.label);%EEG.ChannelNumber;
                hdr.nTrials = 1;
                hdr.label = data.label;%EEG.ChannelTitles;
                hdr.chantype = repmat({'eeg'},hdr.nChans,1);%Nx1 cell-array with the channel type, see FT_CHANTYPE
                hdr.chanunit = repmat({'uV'},hdr.nChans,1);%EEG.ChannelUnits;%Nx1 cell-array with the physical units, see FT_CHANUNIT
  
                sigpositive_data = data.trial{1};
                sigpositive_data(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(sigpositive_data,1),1);
                data.trial{1} = sigpositive_data;
                
                hdr.brainvision_outformat = 'float32';%float32 int16 int32;
                data_file_name_current = [data_file_name_current '_' 'lindev'];
                ft_write_data([data_file_name_current '.eeg'], data.trial{:},'dataformat',dataformat,'header',hdr);
            end
        end
        
        data_filt = {};
        
        fprintf('dataset %i, split %i: apply filtering to data\n',iData,iSplitFactor);

        for iChannelbyOrder = 1:length(data.label)
            %iChannelbyOrder = 1
            curr_channel_name = data.label(iChannelbyOrder);
            
            cfg = [];
            cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
            
            curr_corresponding_setting_row = find(strcmp(curr_channel_name,curr_channel_settings_table_analysis.new_channel_label),1,'first');
            cfg.hpfilter = char(curr_channel_settings_table_analysis.high_pass(curr_corresponding_setting_row));
            cfg.hpfiltord = double(curr_channel_settings_table_analysis.high_pass_order(curr_corresponding_setting_row));
            cfg.hpfreq = [double(curr_channel_settings_table_analysis.high_pass_freq(curr_corresponding_setting_row))];%dummy values are overwritten by low level function
            cfg.hpfilttype = char(curr_channel_settings_table_analysis.high_pass_type(curr_corresponding_setting_row));
            cfg.hpfiltdir = 'twopass';
            
            cfg.lpfilter = char(curr_channel_settings_table_analysis.low_pass(curr_corresponding_setting_row));
            cfg.lpfiltord = double(curr_channel_settings_table_analysis.low_pass_order(curr_corresponding_setting_row));
            cfg.lpfreq = [double(curr_channel_settings_table_analysis.low_pass_freq(curr_corresponding_setting_row))];%dummy values are overwritten by low level function
            cfg.lpfilttype = char(curr_channel_settings_table_analysis.low_pass_type(curr_corresponding_setting_row));
            cfg.lpfiltdir = 'twopass';
            
            cfg.channel = [iChannelbyOrder];
            data_filt{iChannelbyOrder} = ft_preprocessing(cfg,data);
            
            
        end
        data = ft_appenddata([],data_filt{:});
        %data.hdr = data_filt{1}.hdr;
        
        sigpositive_data = data.trial{1};
        
        
        sigpositive_data(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(sigpositive_data,1),1);
        data.trial{1} = sigpositive_data;
        
        %data_file_name_current = 'test'
        dataformat = 'edf';
        hdr.edf_doautoscale = false;
        hdr.edf_accuracy = 0.1;
        hdr.edf_docutoff = true;
        data_file_name_current = [data_file_name_current '_' 'filtered' '_' '0p1uVaccurracy_and_Ycutoff'];
        ft_write_data([data_file_name_current '.edf' ], data.trial{:},'dataformat',dataformat,'header',hdr);
        
     
        
        
%         data_file_name_current1 = [data_file_name_current '_' 'filtered.1'];
%         hdr.edf_accuracy = 0.1;
%         hdr.edf_doautoscale = false;
%         hdr.edf_docutoff = true;
%         ft_write_data([data_file_name_current1 '.edf' ], data.trial{:},'dataformat','edf','header',hdr);
%         
%          data_file_name_current1 = [data_file_name_current '_' 'filtered.autoscale'];
%         hdr.edf_doautoscale = true;
%         ft_write_data([data_file_name_current1 '.edf' ], data.trial{:},'dataformat','edf','header',hdr);
%         
%         hdr.brainvision_outformat = 'int16';
%         ft_write_data([data_file_name_current '.eeg' ], data.trial{:},'dataformat','brainvision_eeg','header',hdr);
        
        %hdr.brainvision_outformat = 'float32'
        %ft_write_data([data_file_name_current '.eeg' ], data.trial{:},'dataformat','brainvision_eeg','header',hdr);


        
        %dataformat = 'edf';
        %%%ft_write_data([dataFilesNames{iData} '.edf'], dat,'dataformat',dataformat,'header',hdr)
        %%%ft_write_data([dataFilesNames{iData} '.edf'], int16(round(dat)),'dataformat',dataformat,'header',hdr)
        %ft_write_data([dataFilesNames{iData} '_filtered.edf'], int16(round(data.trial{:})),'dataformat',dataformat,'header',data.hdr)
    end
end
end