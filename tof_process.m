function data_set = tof_process(data_set, varargin)

%% Definitions

import akpack.J2eV

R = 8.314;

% Suppress warning about overwriting confidence bounds
warning( 'off', ...
    'curvefit:cfit:subsasgn:coeffsClearingConfBounds' );

%% Input Parsing

p = inputParser();
p.FunctionName = 'tof_process';
p.addOptional('StartT', [141, 516, 122])
p.addOptional('StartV', [1176, 0, 0])
p.addOptional('NumberOfChannels', 3)
p.addOptional('DistanceToSpectrometer', 0.075)
p.addOptional('UseExternalDistance', true)
p.addOptional('MoleculeMassAMU', 28)
p.addOptional('TolT', 0.25)
p.addOptional('TolV', 0.25)
p.addOptional('TolD', 0.20)
p.addOptional('StartIndex', 1)
p.addOptional('EndIndex', [])
p.addOptional('ShowFit', false)

parse(p, varargin{:});
r = p.Results;

%% Start Values and Bounds
s.b = r.MoleculeMassAMU * 1e-3 ./ (2*R);
s.a = zeros(1, r.NumberOfChannels);
s.t = r.StartT;
s.v = r.StartV;
s.d = r.DistanceToSpectrometer;

l.a = zeros(1, r.NumberOfChannels);
l.t = s.t - s.t * r.TolT;
l.v = s.v - s.v * r.TolV;
l.b = s.b;
l.d = s.d - s.d * r.TolD;

u.a = zeros(1, r.NumberOfChannels) + inf;
u.t = s.t + s.t * r.TolT;
u.v = s.v + s.v * r.TolV;
u.b = s.b;
u.d = s.d + s.d * r.TolD;

%% Write fit options
for i=1:numel(data_set)
    if isempty(r.EndIndex)
        end_index = numel(data_set(i).time);
    end
    data_set(i).fitopts.start = s;
    data_set(i).fitopts.lower = l;
    data_set(i).fitopts.upper = u;
    data_set(i).fitopts.startindex = r.StartIndex;
    data_set(i).fitopts.endindex = end_index;
    data_set(i).fitopts.tolt = r.TolT;
    data_set(i).fitopts.tolv = r.TolV;
    data_set(i).fitopts.told = r.TolD;
end

%% Fitting

ft = make_tof_fittype(r.NumberOfChannels);

disp('Fitting ...')
for i = 1:numel(data_set)
    fprintf('#')
    
    % Get distance to mass spectrometer from data set entry if
    % 'UseExternalDistance' is not set to false
    if isfield(data_set(i), 'externalDistance') ...
            && r.UseExternalDistance ...
            && ~isnan(data_set(i).externalDistance)
        s.d = tof_distance(data_set(i).externalDistance)*1e-3;
        l.d = s.d - s.d * r.TolD;
        u.d = s.d + s.d * r.TolD;
    end
    
    end_index = data_set(i).fitopts.endindex;
    x = data_set(i).time(r.StartIndex:end_index);
    y = data_set(i).flux(r.StartIndex:end_index);
    
    [obj, gof, out] = fit(x, y, ft, ...
        'Start', [s.a, s.t, s.b, s.d, s.v], ...
        'Lower', [l.a, l.t, l.b, l.d, l.v], ...
        'Upper', [u.a, u.t, u.b, u.d, u.v], ...
        'MaxIter', 400, ...
        'TolFun', 1e-18, ...
        'TolX', 1e-32, ...
        'MaxFunEvals', 100);
    data_set(i).fitobj = obj;
    data_set(i).fitgof = gof;
    data_set(i).fitout = out;
end
fprintf('\nDone.\n\n')

%% View Fit Results

for coeff = (coeffnames(ft))'
    disp(coeff{1})
    v = zeros(1, numel(data_set));
    for i = 1:numel(data_set)
        v(i) = data_set(i).fitobj.(coeff);
        fprintf('%g, ', v(i))
    end
    
    fprintf('\nMean Value: %1.f', mean(v))
    fprintf('\n------------------------------------\n')
end

disp('Goodnes of Fits (rsquare):')
gof = zeros(1, numel(data_set));
for i = 1:numel(data_set)
    gof(i) = data_set(i).fitgof.rsquare;
    disp(gof(i))
end
fprintf('Mean: %f\n\n', mean(gof))

disp('Translational Energies:')
for i=1:numel(data_set)
    [e_trans(i,:), t_trans(i,:)] = translational_energy(data_set(i).fitobj); %#ok
    data_set(i).e_trans = e_trans(i,:);
    data_set(i).t_trans = t_trans(i,:);
    fprintf('%s:\n', data_set(i).fileName)
    fprintf('Energies: \t%0.3f meV\n', J2eV(e_trans(i,:)) * 1e3)
    fprintf('Temperatures: \t%0.3f K\n', t_trans(i,:))
    disp('------------------------------------------------------')
end
fprintf('\nMean Values:\n')
fprintf('\tEnergy: %19.3g meV\n', J2eV(mean(e_trans, 1)) * 1e3)
fprintf('\tTemperature: %0.f K\n', mean(t_trans, 1))

fprintf('\nIntegrals:\n')
for i=1:numel(data_set)
    
    fprintf('%s:\n', data_set(i).fileName)
    
    data_set(i).integral_total = integral(@(x) data_set(i).fitobj(x), ...
        0, inf, 'ArrayValued', true);
    
    fprintf('Total: %g\n', data_set(i).integral_total)
    
    s = data_set(i).fitobj;
    for j = 1:r.NumberOfChannels
        % Set all amplitudes to 0
        s.(sprintf('A%.2u', j)) = 0;
    end
    
    for j = 1:r.NumberOfChannels
        % Restore
        s.(sprintf('A%.2u', j)) = data_set(i).fitobj.(sprintf('A%.2u', j));
        data_set(i).integral_ch(j) = integral(@(x) s(x), ...
            0, inf, 'ArrayValued', true);        
        data_set(i).int_ch_frac(j) = data_set(i).integral_ch(j) ./ ...
            data_set(i).integral_total;
        
        % Get error
        errors = confint(data_set(i).fitobj);
        
        for k=1:2
            error_amplitude = errors(k,j);
            
            % Calculate integrals with error
            s.(sprintf('A%.2u', j)) = error_amplitude;
            data_set(i).integral_ch_error(k,j) = integral(@(x) s(x), ...
                0, inf, 'ArrayValued', true);
            data_set(i).int_ch_frac_error(k,j) = ...
                data_set(i).integral_ch_error(k,j) ./ ...
                data_set(i).integral_total;
        end
        
        % Reset to 0
        s.(sprintf('A%.2u', j)) = 0;
        
        fprintf('Channel %.2u: %.2f (%.2f, %.2f)\n', ...
            j, data_set(i).int_ch_frac(j), ...
            data_set(i).int_ch_frac_error(1,j), ...
            data_set(i).int_ch_frac_error(2,j))
    end
    
    fprintf('--------------------------------------------------\n')
    
end


%% View Data

if r.ShowFit
    colors = lines();
    
    for i = 1:numel(data_set)
        figure('WindowStyle', 'normal')
        subplot(2, 1, 1)
        hold on
        plot(data_set(i).time, data_set(i).flux, 'k.')
        fit_plot = plot(data_set(i).fitobj);
        fit_plot.Color = 'k';
        fit_plot.DisplayName = 'Sum of Channels';
        ylim([-0.25e18, max(data_set(i).flux(50:end))])
        title(sprintf('%s', data_set(i).fileName), 'Interpreter', 'none')
        
        % Plot Single Channels
        s = data_set(i).fitobj;
        for j = 1:r.NumberOfChannels
            % Set all amplitudes to 0
            s.(sprintf('A%.2u', j)) = 0;
        end
        
        for j = 1:r.NumberOfChannels
            % Restore
            s.(sprintf('A%.2u', j)) = data_set(i).fitobj.(sprintf('A%.2u', j));
            
            ch_plot = plot(s);
            ch_plot.Color = colors(j,:);
            ch_plot.DisplayName = sprintf('Channel %.2u', j);
            
            % Reset to 0
            s.(sprintf('A%.2u', j)) = 0;
        end
        
        % Plot Residuals
        subplot(2, 1, 2)
        hold on
        stem(data_set(i).time(r.StartIndex:r.EndIndex), ...
            data_set(i).fitout.residuals, 'k.')
        plot(data_set(i).time(r.StartIndex:r.EndIndex), ...
            smooth(data_set(i).fitout.residuals, ...
            ceil((r.EndIndex - r.StartIndex)/10), 'loess'), 'r')
        for j = 1:2
            subplot(2, 1, j)
            xlim([0, 0.4e-3])
        end
    end
end

end

%% HELPER FUNCTIONS

function ft = make_tof_fittype(numChannels)
% Generate a fit model based on the number of channels

fitTypeString = '';
for i = 1:numChannels
    if i > 1
        fitTypeString = [fitTypeString, '+']; %#ok
    end
    fitTypeString = [...
        fitTypeString, ...
        sprintf('A%.2u / x^5*exp(-b/T%.2u*( d/x - v%.2u )^2)', i, i, i)...
        ]; %#ok
end

ft = fittype(fitTypeString);
end

function [e, temperature] = translational_energy(fit_obj)

    coeffs = coeffnames(fit_obj);
    % Get Temperature fields
    temp_fields = coeffs(contains(coeffs, 'T'));
    % Get velocity fields
    vel_fields = coeffs(contains(coeffs, 'v'));
    % Get the mass
    mass = fit_obj.b * 2 * 1.3807e-23;
    
    for i=1:numel(temp_fields)
        temp(i) = fit_obj.(temp_fields{i}); %#ok
        vel(i) = fit_obj.(vel_fields{i}); %#ok
    end
    
    e = 1/2 * mass .* vel.^2 + 2 * 1.3807e-23 .* temp;
    temperature = e ./ (2 * 1.3807e-23);
    
end