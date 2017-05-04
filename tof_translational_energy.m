function [e, temperature] = tof_translational_energy(fit_obj)

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