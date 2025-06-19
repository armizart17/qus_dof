function [ps_data_comp, ps_data_ang] = bf_planewave_comp_obliq(raw_data_all, fs, fnumber, angles_deg, pars)
% [ps_data_comp, ps_data_ang] = bf_planewave_comp_obliq(raw_data_all, fs, fnumber, angles_deg, pars) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BF_PLANEWAVE_COMP  Plane-wave compounding with oblique imaging-plane support
% INPUTS:
%   raw_data_all : (nSamples*nAngles)×nElements RF matrix, stacked by angle
%   fs           : sampling frequency [Hz]
%   fnumber      : receive F-number for dynamic aperture
%   angles_deg   : 1×nAngles transmit steering angles [deg]
%   pars         : struct with fields
%                  • pitch        – element spacing [m]
%                  • sos          – speed of sound [m/s]
%                  • x            – (opt.) lateral element positions [m]
%                  • scanAngleDeg – (opt.) imaging-plane rotation [deg] (default=0)
%
% OUTPUTS:
%   ps_data_comp : nSamples×nElements compounded image
%   ps_data_ang  : nSamples×nElements×nAngles individual angle images
%
% EXAMPLE:
%   pars.pitch       = 0.3e-3;
%   pars.sos         = 1540;
%   pars.scanAngleDeg= 15;           % tilt imaging plane by 15°
%   [C, A] = bf_planewave_comp(RF, 40e6, 2, [-10 0 10], pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%— unpack and sanity ---------------------------------------------------
  nAngles = numel(angles_deg);
  [nTotal, nElem] = size(raw_data_all);
  nSamples = nTotal / nAngles;
  if mod(nTotal, nAngles)~=0
    error('raw_data_all has %d rows but %d angles → not divisible.', nTotal, nAngles);
  end

  %— lateral element positions -------------------------------------------
  if isfield(pars,'x')
    x_el = pars.x(:).';
  else
    x_el = (0:num_elements-1) * pars.pitch;
    x_el = x_el - mean(x_el);

  end

  %— depth vector and base grid ------------------------------------------
  c    = pars.sos;
  z    = (0:nSamples-1).' * (c/(2*fs));      % nSamples×1
  [XX,ZZ] = meshgrid(x_el, z);               % nSamples×nElem

  %— allocate outputs ----------------------------------------------------
  ps_data_comp = zeros(nSamples, nElem);
  ps_data_ang  = zeros(nSamples, nElem, nAngles);

  %— loop over each angle -----------------------------------------------
  for a = 1:nAngles
    th = deg2rad(angles_deg(a));            % steering/rotation angle
    % rotate imaging grid
    Xq = XX*cos(th) - ZZ*sin(th);
    Zq = XX*sin(th) + ZZ*cos(th);

    % slice out this angle's RF block
    i0   = (a-1)*nSamples + 1;
    i1   =  a   *nSamples;
    data = raw_data_all(i0:i1, :);

    % beamform for this angle
    ps_single = zeros(nSamples, nElem);
    for j = 1:nElem
      for i = 1:nSamples
        x_pt = Xq(i,j);
        z_pt = Zq(i,j);

        % transmit+receive path delay
        d_ec = z_pt*cos(th) + x_pt*sin(th);
        tau  = (d_ec + sqrt(z_pt^2 + (x_pt - x_el).^2)) / c;

        % sample‐index
        samp = round(tau*fs);
        samp = min(max(samp,1), nSamples);
        idx  = samp + nSamples*(0:nElem-1);
        vals = data(idx);

        % dynamic aperture
        ap_rad = z_pt/(2*fnumber + eps);
        apod   = abs(x_el - x_pt) <= ap_rad;

        ps_single(i,j) = sum(vals .* apod);
      end
    end

    ps_data_ang(:,:,a) = ps_single;
    ps_data_comp      = ps_data_comp + ps_single;
  end

  %— average to get the compounded image -------------------------------
  ps_data_comp = ps_data_comp / nAngles;
end