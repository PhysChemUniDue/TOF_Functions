function fh = tof_show(data_set, plotResiduals)

colors = lines();

for i = 1:numel(data_set)
    
    fh(i) = figure('WindowStyle', 'normal'); %#ok
    if isfield(data_set(i), 'fitobj') && plotResiduals
        subplot(2, 1, 1)
    end
    hold on
    plot(data_set(i).time, data_set(i).flux, ...
        'Color', [0.7 0.7 0.7], ...
        'Marker', '.', ...
        'LineStyle', 'none')
    
    set(gca,  'XLim', [0, 0.5e-3])
    
    if isfield(data_set(i), 'fitobj')
        
        num_channels = (numel(coeffnames(data_set(1).fitobj)) - 2) / 3;
        
        fit_plot = plot(data_set(i).fitobj);
        fit_plot.Color = 'k';
        fit_plot.DisplayName = 'Sum of Channels';
        
        % Plot Single Channels
        s = data_set(i).fitobj;
        for j = 1:num_channels
            % Set all amplitudes to 0
            s.(sprintf('A%.2u', j)) = 0;
        end
        
        for j = 1:num_channels
            % Restore
            s.(sprintf('A%.2u', j)) = data_set(i).fitobj.(sprintf('A%.2u', j));
            
            ch_plot = plot(s);
            ch_plot.Color = colors(j,:);
            ch_plot.DisplayName = sprintf('Channel %.2u', j);
            
            % Reset to 0
            s.(sprintf('A%.2u', j)) = 0;
        end
        
    end
    
    ylim([-max(data_set(i).flux(40:end)) * 0.05, ...
        max(data_set(i).flux(40:end)) * 1.05])
    title(sprintf('%s', data_set(i).fileName), 'Interpreter', 'none')
    xlabel('Time in s')
    ylabel('Flux (arb. units)')
    
    %% Plot Residuals
    
    if isfield(data_set(i), 'fitobj') && plotResiduals
        
        start_index = data_set(i).fitopts.startindex;
        end_index = data_set(i).fitopts.endindex;
        
        subplot(2, 1, 2)
        hold on
        stem(data_set(i).time(start_index:end_index), ...
            data_set(i).fitout.residuals, ...
            'Color', [0.7 0.7 0.7], ...
            'Marker', '.')
        plot(data_set(i).time(start_index:end_index), ...
            smooth(data_set(i).fitout.residuals, ...
            ceil((numel(data_set(i).fitout.residuals))/10), 'loess'), 'r')
        
        xlabel('Time in s')
        ylabel('Residuals (arb. units)')
    else
        end_index = numel(data_set(i).time);
    end
    
    %% Set x axis
    
    if plotResiduals
        for j = 1:2
            subplot(2, 1, j)
            xlim([0, data_set(i).time(end_index)])
        end
    else
        xlim([0, data_set(i).time(end_index)])
    end
end

end