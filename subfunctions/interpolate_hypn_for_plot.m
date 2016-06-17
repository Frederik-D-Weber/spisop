function [hypn_plot_interpol hypn_plot_interpol_MA] = interpolate_hypn_for_plot(hypn,epochLengthSamples,plot_MA_offset)

        plot_MA_offset = -5.5;
        hypn_plot = hypn;
        hypn_plot(hypn_plot(:,1) == 5,1) = 0.5;
        hypn_plot(hypn_plot(:,1) == 8,1) = 0;
        hypn_plot(:,1) = hypn_plot(:,1)*-1;
        hypn_plot_MA = hypn_plot(:,2) ;
        hypn_plot_MA = hypn_plot_MA*0.5;
        hypn_plot_MA(hypn_plot_MA > 1) = 1.35;
        hypn_plot = hypn_plot(:,1) ;
        hypn_plot_interpol = [];
        hypn_plot_interpol_MA = [];
        for iEp = 1:length(hypn_plot)
            temp_samples = repmat(hypn_plot(iEp),epochLengthSamples,1);
            if (hypn_plot(iEp) == -0.5) %REM
                temp_samples(1:2:end) = -0.70;
                temp_samples(2:2:end) = -0.20;
%                 for iSamp = 1:length(temp_samples)
%                     if mod(iSamp,2) == 0
%                         temp_samples(iSamp) = -0.20;
%                     else
%                         temp_samples(iSamp) = -0.70;
%                     end
%                 end
            end
            
            hypn_plot_interpol = [hypn_plot_interpol; temp_samples];
            
            temp_samples_MA = repmat(plot_MA_offset+hypn_plot_MA(iEp),epochLengthSamples,1);
            if (hypn_plot_MA(iEp) > 0) %REM
                for iSamp = 1:length(temp_samples_MA)
                    if mod(iSamp,2) == 1
                        temp_samples_MA(iSamp) = plot_MA_offset;
                    end
                end
            end
            hypn_plot_interpol_MA = [hypn_plot_interpol_MA; temp_samples_MA];
            
        end

end
