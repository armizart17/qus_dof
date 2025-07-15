[xPolar2,zPolar2,z02,xFull2,zFull2,r02] = toPolarGridRec(siz,zmax,param);                    
 


function [xPolar,zPolar,z0,xFull,zFull,r0] = toPolarGridRec(siz,zmax,param)                     
    N = param.Nelements;
    R = param.radius;
    p = param.pitch;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+p,-z0+zmax,siz(1)));
    [xPolar,zPolar] = pol2cart(th,r);

    zPolar = zPolar+z0;

    % NEW EMZ
    th_r = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
    r_r = linspace(R+p,-z0+zmax,siz(1));
    
    xFull = th_r; % [deg]
    r0 = r_r(1);
    zFull = (r_r-r0); % [m]

end
%%



%%
% Get BMode & RF @EMZ
sampleFreq = Receive(1).decimSampleRate*1e6;
rxSamples = round(rx_delays*sampleFreq);
rxSamples(rxSamples<1) = 1;
rxSamples(rxSamples>num_samples) = num_samples;

apodAperture = dyn_aperture;

centralFreq = Receive(1).demodFrequency*1e6;
tgcDbCm = 0.1 * (centralFreq * 1e-6) ^ 1;
tgcNpM = tgcDbCm / 8.686 * 100;
r = exp(tgcNpM * z*1e-3); % [m]

plotDepth = floor(P.endDepth*2*4);

[bMode, rf] = getBModeC52v(rf_channel(:, :, :), rxSamples, apodAperture, num_samples, ...
    n_pulses, n_elements, plotDepth, r);
   
t = (0:(num_samples-1))/fs; % [sec.] time domain 0:T:(N_sample-1)*T
z = sound_speed*t/2;
z2 = (0:size(bMode, 1) - 1)/sampleFreq*sound_speed/2; % [m]
   

function [bMode, rf_data] = getBModeC52v(rf_channel, rx_samples, apodAperture, n_samples, n_pulses, n_elements, plotDepth, tgcr)
    rf_data = zeros(n_samples, n_pulses);
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            rf_channel(:, e, n) = rf_channel(rx_samples(:, e, n), e, n);
        end
        % Apodization
        % rf_channel(:, :, n) = rf_channel(:, :, n) .* apodAperture(n, :);

         rf_channel(:, :, n) = rf_channel(:, :, n) .* apodAperture(:, :, n);
   
        % Summing
        rf_data(:, n) = sum(rf_channel(:, :, n), 2);
    end

    % Depth segmentation
    rf_data = rf_data(1:plotDepth, :);

    % TGC
    tgcr = tgcr(1:plotDepth);
    rf_data = bsxfun(@times, tgcr', rf_data);

    rf_data_padded = zeros(floor(size(rf_data, 1)*1.5), size(rf_data, 2));
    rf_data_padded(1:size(rf_data, 1), :) = rf_data;

    % BMode
    bMode = abs(hilbert(rf_data_padded));
    bMode = 20*log10(bMode);
    bMode = bMode(1:size(rf_data ,1), :);
    bMode = bMode-max(bMode(:));
end