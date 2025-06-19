function [ps_data_comp, ps_data_ang] = bf_planewave_comp(raw_data_all, fs, fnumber, angles_deg, pars)
% function ps_data_comp = bf_planewave_comp(raw_data_all, fs, fnumber, angles_deg, pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs plane wave compounding from a vertically stacked RF matrix
%
% INPUTS:
% - raw_data_all : (nSamples * nAngles) x nElements matrix
% - fs           : sampling frequency [Hz]
% - fnumber      : dynamic aperture control
% - angles_deg   : vector of transmit angles [deg] (length = nAngles)
% - pars         : struct with fields pitch, sos, and optional x
%
% OUTPUT:
% - ps_data_comp : compounded beamformed image [nSamples x nElements]
% - ps_data_ang  : Beamformed image per angle [nSamples x nElements x nAngles] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_angles = numel(angles_deg);
[num_samples_total, num_elements] = size(raw_data_all);
num_samples = num_samples_total / n_angles;

% Validation
if mod(num_samples_total, n_angles) ~= 0
    error('Total number of rows in raw_data_all is not divisible by number of angles.');
end

% Lateral positions
if isfield(pars, 'x')
    x = pars.x;
else
    x = (0:num_elements-1) * pars.pitch;
    x = x - mean(x);  
end

tx_positions = x;
c = pars.sos;

ps_data_comp = zeros(num_samples, num_elements);
ps_data_ang  = zeros(num_samples, num_elements, n_angles);

for i_a = 1:n_angles

    % Extract current RF block
    ini_idx   = (i_a-1)*num_samples + 1;
    end_idx   = i_a * num_samples;
    data      = raw_data_all(ini_idx:end_idx, :);
    angle_rad = deg2rad( angles_deg(i_a) );
    
    % ps_single = zeros(num_samples,num_elements); %ps_single = ps_data_ang

    for j = 1:num_elements  % lateral (x)
        x_point = x(j);
        for i = 1:num_samples  % axial (z)
            z_point = c * (i / fs) / 2;

            % Compute transmit and receive delays
            d_ec = z_point * cos(angle_rad) + x_point * sin(angle_rad);
            tau  = (d_ec + sqrt(z_point^2 + (x_point - tx_positions).^2)) / c;

            % Sample index calculation
            samples = round(tau * fs);
            samples = max(min(samples, num_samples), 1);

            indices = samples + num_samples * (0:num_elements - 1);
            values  = data(indices);

            % Dynamic aperture (F-number)
            a      = z_point / (2 * fnumber + 1e-10);
            liminf = x_point - a;
            limsup = x_point + a;
            apod   = (tx_positions >= liminf) & (tx_positions <= limsup);

            ps_data_ang(i, j, i_a) = sum(apod .* values);
        end
    end
    
    ps_data_comp = ps_data_comp + ps_data_ang(:,:,i_a);
end

% Average?
% ps_data_comp = ps_data_comp / n_angles;
end
