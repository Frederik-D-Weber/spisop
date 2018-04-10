% Requires cfslib MATLAB library: https://github.com/amiyapatanaik/cfslib-MATLAB

classdef AutoSleepScoringZ3score
    properties
        cfg
        hypnogram
        hypnogram_confidence
        artifacts
        data
        isOperable
        isScored
    end
    methods
        function obj = AutoSleepScoringZ3score(data,cfg)
            obj.data = data;
            if isempty(cfg)
                obj.cfg = [];
                obj.cfg.autoSleepsCoreAlgoserverURL = 'http://z3score.com/api/v1';
                obj.cfg.email = 'example@email.com';
                obj.cfg.key = '1234567890abcdefghijklmnopqrstuvwxyz#+##qe/=';
                if exist('Z3Score_license.txt') == 2
                    try
                        licfile = read_mixed_csv('Z3Score_license.txt',',');
                        obj.cfg.autoSleepsCoreAlgoserverURL = licfile{1};
                        obj.cfg.email = licfile{2};
                        obj.cfg.key = licfile{3};
                    catch err
                    end
                end
            else
                obj.cfg = cfg;
            end
            obj.cfg.epoch_duration_seconds = 30;
            
            obj.hypnogram = [];
            obj.hypnogram_confidence = [];
            obj.artifacts = [];
            
            obj.isOperable = false;
            obj.isScored = false;
            
            obj = obj.update_cfg();
            obj = obj.update_channel();
            obj = obj.check();
        end
        function obj = update_cfg(obj)
            prompt = {'Auto Sleep Score Algo Server URL', 'user email', 'API key'};
            title = 'Credentials Z3Score';
            dims = [1 35];
            definput = {obj.cfg.autoSleepsCoreAlgoserverURL, obj.cfg.email, obj.cfg.key};
            cfg_answers = inputdlg(prompt,title,dims,definput);
            if ~isempty(cfg_answers)
                obj.cfg.autoSleepsCoreAlgoserverURL = cfg_answers{1};
                obj.cfg.email = cfg_answers{2};
                obj.cfg.key = cfg_answers{3};
            else
                disp(['Z3Score: abort did not update config!'])
            end
        end
        function obj = update_channel(obj)
            if ~isfield(obj.cfg,'channel')
                obj.cfg.channel = [];
            end
            %% choose channels
            if ~isfield(obj.cfg.channel,'C3_A2_idx')
                obj.cfg.channel.C3_A2_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'C4_A1_idx')
                obj.cfg.channel.C4_A1_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'EOGl_A1_idx')
                obj.cfg.channel.EOGl_A1_idx = 1;
            end
            if ~isfield(obj.cfg.channel,'EOGr_A2_idx')
                obj.cfg.channel.EOGr_A2_idx = 1;
            end
            channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                        obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx];
            channel_list_request = {'C3:A2', 'C4:A1', 'EOGl:A1', 'EOGr:A2'};
            list_channel = obj.data.label;
            for iReqCh = 1:numel(channel_list_request)
            [selection,ok] = listdlg('PromptString',['Choose ' channel_list_request{iReqCh} ' channel'],...
                'SelectionMode','single',...
                'InitialValue',channels_idx(iReqCh),...
                'ListString',list_channel);
            if ~ok
                disp([channel_list_request{iReqCh} ' channel has not been updated'])
                %return
            end
            switch iReqCh
                case 1
                    obj.cfg.channel.C3_A2_idx = selection;
                case 2
                    obj.cfg.channel.C4_A1_idx = selection;
                case 3
                    obj.cfg.channel.EOGl_A1_idx = selection;
                case 4
                    obj.cfg.channel.EOGr_A2_idx = selection;
                otherwise
            end
            end

        end
        function obj = update_data(obj,data)
            obj.data = data;
            obj = obj.update_channel();
        end
        function message = status(obj)
            %check license
            try
                response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL '/check'],...
                    {'email',obj.cfg.email,'key',obj.cfg.key}));
            catch
                message = 'Z3Score: scoring server is unreachable';
                disp(message);
                return
            end
            
            if response.status == 0,
                message = 'Z3Score: License check failed';
                disp(message);
                disp(['Z3Score: Error message: ' response.message])
                return
            end
            
            %disp(response.message);
            
            num_epochs_requested = floor(size(obj.data.trial{1},2)/obj.data.fsample/obj.cfg.epoch_duration_seconds);

            message = [sprintf('Z3Score: Server reachabel\nAPI Call limit left (hourly): %d\n Epoch limit left (daily): %d\n ... next request %d epochs',response.call_limit, response.epoch_limit,num_epochs_requested)];
 
            
        end
        function obj = check(obj)
            obj.isOperable = false;
            %check license
            try
                response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL '/check'],...
                    {'email',obj.cfg.email,'key',obj.cfg.key}));
            catch
                disp('Z3Score: scoring server is unreachable');
                return
            end
            
            if response.status == 0,
                disp('Z3Score: License check failed');
                disp(['Z3Score: Error message: ' response.message])
                return
            end
            
            %disp(response.message);
            
            num_epochs_requested = floor(size(obj.data.trial{1},2)/obj.data.fsample/obj.cfg.epoch_duration_seconds);

            if ~(response.call_limit >= 1) && ~((response.epoch_limit-num_epochs_requested) > 0)
                msgbox('Limit reached' ,sprintf('Z3Score: API Call limit left (hourly): %d\n Epoch limit left (daily): %d\n ... but requested %d',response.call_limit, response.epoch_limit,num_epochs_requested),'modal');
                return
            end
            
            
            if isempty(obj.cfg.channel.C3_A2_idx) &&...
                isempty(obj.cfg.channel.C4_A1_idx) &&...
                isempty(obj.cfg.channel.EOGl_A1_idx) &&...
                isempty(obj.cfg.channel.EOGr_A2_idx)
                %size(obj.data.trial{1},1) ~= 4
                disp('Z3Score: data does not provide the 4 needed channels')
                return
            end
            obj.isOperable = true;
        end
        function obj = score(obj)
            if obj.isOperable
                if ~isfield(obj.cfg,'stream')
                    obj.cfg.stream = [];
                end
                if isempty(obj.cfg.stream)
                    %Find out sampling rate
                    samplingRate = obj.data.fsample;
                    %Construct raw data from selected channels
                    channels_idx = [obj.cfg.channel.C3_A2_idx obj.cfg.channel.C4_A1_idx...
                        obj.cfg.channel.EOGl_A1_idx obj.cfg.channel.EOGr_A2_idx];
                    EEGData = obj.data.trial{1}(channels_idx,:);
                    
                    num_epochs = floor(size(EEGData,2)/samplingRate/obj.cfg.epoch_duration_seconds);
                    
                    
                    
                    %Convert raw stream to a CFS and write to a file
                    disp('Z3Score: data to CFS stream');
                    %    tic;
                    %Convert raw PSG stream to CFS stream
                    obj.cfg.stream = streamCFS(EEGData, samplingRate);
                    EEGData = [];
                    %fileID = fopen('test.cfs','w');
                    %fwrite(fileID,stream,'*uint8','ieee-le');
                    %fclose(fileID);
                    %t = toc;
                    %fprintf('Time taken %.3f seconds\n',t);
                end
                
                disp('Z3Score: fetch scoring info from server');
     
                try
                    response = loadjson(urlreadpost([obj.cfg.autoSleepsCoreAlgoserverURL  '/score'], ...
                        {'email',obj.cfg.email,'key',obj.cfg.key,'file',obj.cfg.stream}));
                catch
                    disp('Z3Score: Scoring server is unreachable');
                    return
                end
                
                if response.status == 0,
                    disp('Z3Score: Error scoring data.');
                    disp(['Z3Score: Error message:' response.message])
                    return
                end
                
                %fprintf('Time taken %.3f seconds.\nAPI calls left (hourly limits): %d, Epochs left (daily limits): %d \n',t, response.calls_left, response.epochs_left);
                %Automatic sleep scores
                scores = response.message;
                obj.hypnogram = [scores(:,1) zeros(numel(scores(:,1)),1)];
                obj.hypnogram_confidence = scores(:,2);
                obj.hypnogram_confidence = linear_scaling_normalization(obj.hypnogram_confidence, 0, 10, 0, 1);
                
                %artifacts are in 5-seconds scoring, only 7 are
                %unscorable/artifact epochs
                artifacts_5sec_blocks = response.artifact;
                artifacts_5sec_blocks_index = (artifacts_5sec_blocks == 7);
                
                epoch_length_seconds = 5; % in seconds
                temp_starting_seconds_in_data = (1:floor(obj.data.sampleinfo(2)/(obj.data.fsample*epoch_length_seconds)))-1;
                epoch_begsample = obj.data.fsample*epoch_length_seconds*temp_starting_seconds_in_data+1;
                epoch_endsample = obj.data.fsample*epoch_length_seconds*(temp_starting_seconds_in_data+1);
                obj.artifacts = [epoch_begsample' epoch_endsample'];
                obj.artifacts = obj.artifacts(artifacts_5sec_blocks_index,:);
                
                %artifacts MA epochs mark
                for iEpoch = 1:size(obj.hypnogram,1)
                    obj.hypnogram(iEpoch,2) = any(artifacts_5sec_blocks_index(((iEpoch-1)*6+1):(iEpoch*6)))*3;
                end
                                
                %Save the sleep scores
                %csvwrite('test_score.csv',response.message)
                %Read expert sleep scores
                %expert = csvread('test_expert.csv');
                %expert = expert(1:num_epochs);
                
                %C = confusionmat(scores(:,1),expert);
                %fprintf('Z3Score: Auto scoring agreement with expert scorer: %.2f%%\n',sum(scores(:,1) == expert)*100/num_epochs);
                %kappa(C);
                
                disp('Z3Score: Done');
                
            else
                obj.check(obj);
            end
            function newvalue = linear_scaling_normalization(toNorm, minInterval, maxInterval, minNormInterval, maxNormInterval)
                newvalue =  ((toNorm - minInterval) / (maxInterval - minInterval)) * (maxNormInterval - minNormInterval) + minNormInterval;
            end
        end
    end
end




