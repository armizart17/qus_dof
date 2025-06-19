function ps_data = bf_planewave(raw_phan, fs, fnumber, steering_angle_deg, pars)
% function ps_data = bf_planewave(raw_phan, fs, fnumber, steering_angle_deg, pars)

[num_samples, num_elements] = size(raw_phan);
data = raw_phan;

if isfield(pars, 'x')
    x = pars.x;
else
    x = (0:num_elements-1)*pars.pitch;
    x = x - mean(x);
end

if isfield(pars, 'sos')
    c = pars.sos;
else
    c = 1540; % [m/s]
end
tx_positions = x;

ps_data = zeros(num_samples, num_elements);

for j = 1:length(x) % x iteration
    for i = 1:num_samples % z iteration

        % x,z position
        x_point = x(j);
        z_point = c * (i / fs) / 2;

        % calculate delays
        ang  = steering_angle_deg*pi/180; % deg2rad(steering_angle_deg)
        d_ec = z_point*cos(ang) + x_point*sin(ang);
        tau  = (d_ec + sqrt(z_point^2 + (x_point - tx_positions).^2)) / c;
        
        % calculate in samples
        samples = round(tau * fs);
        samples(samples <= 0) = 1;
        samples(samples > num_samples) = 1;
        
        indices = samples + num_samples * (0:num_elements - 1);
        values  = data(indices);
        
        % Sum in aperture
        F = fnumber;
        a = z_point /( (2 * F) + 1e-10);
        liminf = x_point - a;
        limsup = x_point + a;

        apodization = (tx_positions >= liminf & tx_positions <= limsup);
        ps_data(i, j) = sum(apodization .* values);
    end
end
end