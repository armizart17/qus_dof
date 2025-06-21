function [freqs, bsc] = duke_bsc_faranTheory(c0, N_freqs, beads_per_mm3, anderson_coeffs, low_freq, high_freq, diam_um)
% function [freqs, bsc] = duke_bsc_faranTheory(c0, N_freqs, beads_per_mm3, anderson_coeffs, low_freq, high_freq, diam_um)
% Adapted from DUKE by EMZ
    freqs = linspace(low_freq, high_freq, N_freqs);   freqs = freqs(:);
    k0s = 2*pi*freqs/c0;
    
    % The scatterer radii in mm:
    % a_values = (38:1:45)/2000; %um a mm, 
    a_values = (diam_um/2)*1e-3; %um a mm, (diameter to radious)
    % The radii pdf:
    n_values = ones(size(a_values));
    n_values = n_values/sum(n_values);
    N_vals = length(a_values);
    
    bsc = zeros(N_freqs, 1);
    
    for index = 1:N_vals
        bsc = bsc + n_values(index)*pi*(a_values(index))^2*ppval(anderson_coeffs, k0s*a_values(index));
    end
    
    scats_per_mm3 = beads_per_mm3/4/pi;
    bsc = scats_per_mm3*bsc;
    
    % Finally, converting to cm^-1:
    bsc = 10*bsc;
end

