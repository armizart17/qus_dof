% ====================================================================== %
% Script for clinical data.
% Created on June 19th, 2024
% Based on SWIFT
% Designed for DoF testing in Raw Data (Healthy Data, first acquisition 2024-09-07)
% Description:
% Reading Channel Data -> Beamforming -> Spectral Power Law
% ====================================================================== %

%% 544 specs
% BSC_CIRS_liver = readmatrix(fullfile(baseDir,"2096-00 Full acoustic characterization data.xlsx"), ...
%     "Range", "A2:B21","Sheet",2);
%%
% setup,
clear, clc
close all

mean2d  = @(x) mean(x(:));
std2d   = @(x) std(x(:));
cv2d    = @(x) 100*std(x(:))/mean(x(:));

% General directory
baseDir     = 'D:\emirandaz\qus\data\liver\healthy'; 

% Sample and Reference directories
sampleDir   = fullfile(baseDir,'samples');
refsDir     = fullfile(baseDir,'refs','joined');

% Outcomes directory
resultsDir  = fullfile(baseDir,'results_LS');
figsDir     = fullfile(baseDir,'figures_LS');
if ~exist(resultsDir) mkdir(resultsDir); end
if ~exist(figsDir) mkdir(figsDir); end


range_bmode     = [-60 0];
range_acs       = [0 1.1];
Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB;
calc2dStats     = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100*std(x(:))/mean(x(:))};
plotBmode       = false;
plotBSCdB       = true;  % plot \Delta b in dB
fontSize        = 14;

% ROI MANUAL 
roi_already     = true;
roisDir         = fullfile(baseDir,'rois');
if ~exist(roisDir) mkdir(roisDir); end

% Read sample files
sampleFiles = dir(fullfile(sampleDir,'*.mat'));

% Specific sample
% acqDir = dir(fullfile(sampleDir,'016-03.mat')); %65*ma
% acqDir = dir(fullfile(sampleDir,'007-05.mat')); %*ma
% acqDir = dir(fullfile(sampleDir,'014-01.mat')); %*ma
% acqDir = dir(fullfile(sampleDir,'016-06.mat')); %emz
% sampleFiles     = dir(fullfile(sampleDir,'020-05.mat')); %dv



%%

for iFile = 1:length(sampleFiles)
% samName = "022-02";
samName = sampleFiles(iFile).name(1:end-4);

%% Spectral Methods estimation hyperparameters

%% SPECTRAL METHOD PARAMETERS

pars.P           = 512; % NFFT for 10wl is 256, 20wl 512
pars.bw          = [1.5 4]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
pars.blocklines  = 8;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

% blocksize       = 12;   % Axial block size in wavelengths
% blocklines      = 8;   % Num of lines, lateral block size
% overlap_pc      = 0.8;
% freq_L          = 1.5e6; 
% freq_H          = 4e6;



% Bandwidth
fixedBW         = true;
ratio           = db2mag(-30);

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting constants
dynRange = [-60,0];
range_acs = [0,1.5];
bsRange = [-15 15];


%% Loading file and variables
load(fullfile(sampleDir,samName+".mat"));
fprintf("Loading sample %s \n", samName)

RcvData = cell2mat(RcvData);
n_frame = size(RcvData,3); % Frame selector
RcvData = RcvData(:, :, n_frame); % Select frame of RF Data

% Additional variables
central_freq = Receive(1).demodFrequency*1e6; % Central frequency of pulse
fs = Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
n_pulses = P.numRays; % number of pulses
n_elements = Trans.numelements; % number of transducer elements
num_samples = Receive(1).endSample - Receive(1).startSample +1; % samples per channel
sound_speed = Resource.Parameters.speedOfSound; % [m/s]
wvl = sound_speed/central_freq; % [m] wavelength in meters
scalemm2wvl = 1/wvl;

% Initialize variables
rf_channel = zeros(num_samples , n_elements, n_pulses);
rx_apods = zeros(1, n_elements, n_pulses);
rf = zeros(num_samples, n_pulses);

%% Organize data
for n = 1:n_pulses % Iterate through pulses
    rf_channel(:, :, n) = RcvData(Receive(n).startSample:Receive(n).endSample, :);
end

% Acommodate to time delays and rf signals
focus = 20/1000;
t = (0:(num_samples-1))/fs; % [sec.] time domain 0:T:(N_sample-1)*T
[rx_delays] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl);

% Dynamic Aperture
f_num = 3;
z = sound_speed*t/2;
elem_pitch = Trans.spacingMm*1e-3;
maxAprSz = 32;
dyn_aperture = zeros(length(z), n_elements, n_pulses);
for n = 1:n_pulses
    for z_i = 1:length(z)
        a = z(z_i)/(2*f_num);
        hlfAprSz = floor(a / elem_pitch);
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        a_i = -hlfAprSz: hlfAprSz;    % aperture indices
        fulAprSz = 2*hlfAprSz + 1;
        aper_center = n;
        aper = aper_center + a_i;
        aper = aper(aper>=1);
        aper = aper(aper<=128);
        dyn_aperture(z_i, aper, n) = 1;
    end
end

%% Beamforming

% Delay-and-sum
for n = 1:n_pulses
    % Delaying
    for e = 1:n_elements
        % rx_apods(e, n);
        rf_channel(:, e, n) = ...
            interp1(t, rf_channel(:, e, n), rx_delays(:,e, n), 'linear', 0);
    end
    rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
   
     % save(fullfile(processed_data_save_path, processed_data_filename), ...
     %     'rf_channel');
    % Summing
    rf(:, n) = sum(rf_channel(:, :, n), 2);
end

%% B-Mode and coordinates

b_mode = 20*log10(abs(hilbert(rf)));
b_mode = b_mode-max(b_mode(:));

param = getparam('C5-2v');
siz = size(rf);
z = sound_speed*t/2;
zmax = z(end);
R = param.radius;
p = param.pitch;
N = param.Nelements;
L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
d = sqrt(R^2-L^2/4); % apothem
z0 = -d;

th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
r = linspace(R+p,-z0+zmax,siz(1));

% To Polar Coordinates
[xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);

%% Selecting ROI
BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));

xFull = th; % [deg]
r0 = r(1);
zFull = (r-r0)*1e2; % [cm]

if ~roi_already
    % first time
    % manual roi
    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(xFull,zFull,BmodeFull,dynRange);
    ylim([0 10])
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfAngle [°]'); ylabel('\bfDepth [cm]');
    title(samName)
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        rect = getrect;
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            disp(rect)
            break
        end
    end
    
    close,

else
% rois already saved execute full code
load(fullfile(roisDir,samName+".mat"),'rect');
    
pars.x_roi     = [rect(1), rect(1)+rect(3)];      % [º]
pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-2; % [m]

%% PACKAGE DATA FOR SPECTRAL FUNCTIONS

% SAMPLE
SAM.z       = zFull*1e-2; % to [m]
SAM.x       = xFull;
SAM.fs      = fs;
SAM.rf      = rf;
SAM.bMode   = BmodeFull;
bmode_sam   = SAM.bMode;
j_sam       = 1.0;

% REFERENCE 544
alpha_ref   = 0.53; % [dB/cm/MHz]
sos_ref     = 1539; % [m/s]
j_ref       = 1.0;

refFiles    = dir([refsDir,'\*.mat']);
% refFiles = refFiles(x,:); % for select specific "x" refFiles
numRefs     = length(refFiles); 
REF         = load( fullfile(refsDir, refFiles(1).name ) );
newrf       = nan([size(REF.rf), numRefs], 'like', REF.rf); % Use 'like' for type consistency
for i = 1:numRefs
    newrf(:,:,i) = load(fullfile(refsDir,refFiles(i).name ), 'rf').rf(:,:,1); % Directly extract rf, avoiding redundant variables
end

REF.rf  = newrf;
REF.acs = alpha_ref;
REF.c0  = sos_ref;

% Just in case
REF.x   = SAM.x;
REF.z   = SAM.z;

%% BMODE CHECK

caption  = samName; 
bmode_ref = my_RF2Bmode(REF.rf(:,:,1));
if (plotBmode)
deadBand = 0.1e-2;
figure,
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1200, 600]); % [x, y, width, height]

subplot(121), 
imagesc(SAM.x, SAM.z*1E2, bmode_sam, range_bmode), 
% axis("image"), 
hold on;
rectangle('Position', [1, 1E2, 1, 1E2].*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [°]'), ylabel('Depth [cm]');
title(caption)
ylim([deadBand*1000 10])
colormap('gray')

subplot(122), 
imagesc(REF.x, REF.z*1E2, bmode_ref, range_bmode), 
% axis("image");
rectangle('Position',  [1, 1E2, 1, 1E2].*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [°]'), ylabel('Depth [cm]');
title('REF')
colormap('gray')
ylim([deadBand*1000 10])
end


%% BW from spectrogram

freq_L      = pars.bw(1)*1e6;
freq_H      = pars.bw(2)*1e6;
blocksize   = pars.blocksize; 
overlap_pc  = pars.overlap;

x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);
dz = (zFull(2) - zFull(1))/100;

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = rf(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% Axial samples
dz = (zFull(2) - zFull(1))/100;
wz = round(blocksize*wl*(1-overlap_pc)/dz ); % Between windows
nz = 2*round(blocksize*wl/dz /2); % Window size

NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot
[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;
figure,
plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
xlim([0, fs/2e6])
hold on
xline(freq_L/1e6, 'k--')
xline(freq_H/1e6, 'k--')
hold off
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
% saveas(gcf,fullfile(figsDir,"sample"+samName+"_Spectrum.png"))
%close

%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

SR_emz = S_sam ./ S_ref;
SR = permute(SR_emz,[3,1,2]); clear SR_emz

%% GENERAL REGULARIZATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-5;
par_rpl.df_op      = 0;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

% mu_rpl_tv          = [1E3; 1E3; 10^4.2];
mu_rpl_tv          = [0.001; 0.001; 0.001];

%% DOF METHODS

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
f       = band(:); 

% In case mismatch of j_sam & j_ref, otherwise is 1
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

% Matrix (Gaussian)
% X = kron( speye(p*q), ones(size(f)) );
% Z = kron( speye(p*q), -f.^2 );
% W = kron( speye(p*q), -4*f.^j_sam );

% Matrix (Power Law)
X = kron( speye(p*q), ones(size(f)) );
Z = kron( speye(p*q), log(f) ); 
W = kron( speye(p*q), -4*f.^j_sam );

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);  

methods          = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
maps_results_dof = cell(1,length(methods));
bsc_results_dof  = cell(1,length(methods));

% PRIORS
% if ~deltaPriorFix
% delta_b_prior       = log(db2pow(cell2mat(bsc_results_all{iSam, iRef}(3))));
% delta_n_prior       = cell2mat(bsc_results_all{iSam, iRef}(4));
% end


% Loop over methods
    for iMet = 1:length(methods)
        estim_method = methods{iMet};
        clear b n a

        if strcmp(estim_method, '3-DoF')

            SR_comp = SR .* comp_freq_a;
            Y = log(SR_comp);

            u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
            [u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

            b = u_opt(1:p*q);
            n = u_opt(p*q+1:2*p*q);
            a = u_opt(2*p*q+1:3*p*q);
            a_Np2dB = Np2dB*Dy*a./dz(:);

        elseif strcmp(estim_method, '2-DoF-a')
            
            delta_alpha_prior = median(a_Np2dB(:));

            comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);                
            SR_comp     = SR .* comp_ref .* comp_freq_a;
            Y = log(SR_comp);

            u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
            [u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);
            
            b = u_opt(1:p*q);
            n = u_opt(p*q+1:2*p*q);
            a_Np2dB = delta_alpha_prior*ones(p*q, 1);

        elseif strcmp( estim_method, '2-DoF-b')
             
            delta_b_prior = log(db2pow( median(b_ratio_dB(:))   ));
            
            comp_ref    = comp_ref_b_bsc(delta_b_prior);
            SR_comp     = SR .* comp_ref .*comp_freq_a;
            Y = log(SR_comp);

            u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
            [u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);
            
            n = u_opt(1:p*q);
            a = u_opt(p*q+1:2*p*q);
            b = delta_b_prior*ones(p*q, 1);
            a_Np2dB = Np2dB*Dy*a./dz(:);

        elseif strcmp( estim_method, '2-DoF-n')

            delta_n_prior = median(n_ratio(:));

            comp_ref    = comp_ref_n_bsc(delta_n_prior, band, p, q);
            SR_comp     = SR .* comp_ref .* comp_freq_a;
            Y = log(SR_comp);

            u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
            [u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

            b = u_opt(1:p*q);
            a = u_opt(p*q+1:2*p*q);
            n = delta_n_prior*ones(p*q, 1);
            a_Np2dB = Np2dB*Dy*a./dz(:);

        end

        % Compute final parameters
        b_ratio     = reshape(exp(b), p, q);
        b_ratio_dB  = 10*log10(b_ratio);
        alpha_ratio = reshape(a_Np2dB, p, q);
        n_ratio     = reshape(n, p, q); 
        acs_sam     = alpha_ratio + alpha_ref;

        % Compute statistics
        [m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam), calc2dStats{2}(acs_sam), calc2dStats{3}(acs_sam));
        [m_b, s_b, cv_b] = deal(calc2dStats{1}(b_ratio_dB), calc2dStats{2}(b_ratio_dB), calc2dStats{3}(b_ratio_dB));
        [m_n, s_n, cv_n] = deal(calc2dStats{1}(n_ratio), calc2dStats{2}(n_ratio), calc2dStats{3}(n_ratio));

        % Save maps results
        % maps_results_dof{iMet} = {acs_sam, b_ratio_dB, n_ratio};

        maps_results_dof{iMet}.alpha = acs_sam;
        maps_results_dof{iMet}.b_dB  = b_ratio_dB;
        maps_results_dof{iMet}.n     = n_ratio;

        % Compute and save BSC results
        b_est = median(b_ratio(:));
        n_est = median(n_ratio(:));
        
        % Delta
        % bsc_est_gauss = b_est .* exp(-n_est .* band.^2);
        bsc_est_powlaw = b_est .* band.^n_est;

        % Reference fully known
        % OPTION A
        % b_sam     = b_ratio * b_ref;
        % n_sam     = n_ratio + n_ref;
        % bsc_sam   = b_sam *(freq_bsc.^n_sam);
        % med_bsc   = median(bsc_sam(:));
        % Delta_b*b_ref*f^(delta_n_prior +n_ref)

        bsc_results_dof{iMet}.bsc_powlaw = bsc_est_powlaw;

        fprintf('-----%s---\n', estim_method);
        fprintf('α_s        : %.4f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
            if plotBSCdB 
        fprintf('Δb [dB]    : %.4f ± %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
            else
        fprintf('Δb         : %.4f ± %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
            end    
        fprintf('Δn         : %.4f ± %.4f, %%CV = %.4f\n', round(m_n, 4), round(s_n, 4), round(cv_n, 4));
        fprintf('--------\n');

    end

%% DISPLAY OVERLAY SIMPLE RECT

methods      = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};
range_depth  = [0 10];
range_acs    = [0.1 1.2];
transparency = 0.65;

acs_sam = alpha_ratio + alpha_ref;

units           = 1E2;
bmodeFull       = bmode_sam;
% colorImg        = bigImg(acs_sam, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
x_img           = spectralData_sam.lateral;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x;
zFull           = SAM.z*units;

[X,Z]   = meshgrid(xFull,zFull);
roi     = X >= x_img(1) & X <= x_img(end) & Z >= z_img(1) & Z <= z_img(end);

%%%%%%%%%%%%%%%%%%%%%%%% alpha %%%%%%%%%%%%%%%%%%%%%%%%
indices_alpha = [1, 3, 4];  % Corresponden a 3-DoF, 2-DoF-b, 2-DoF-n

figure;
tiledlayout(1,3, 'TileSpacing', 'tight')
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 600]); % [x, y, width, height]
for i = 1:length(indices_alpha)
    idx = indices_alpha(i);
    nexttile;

    [~,~,hColor] = imOverlayInterp(bmodeFull, maps_results_dof{idx}.alpha, range_bmode, range_acs, ...
                                   transparency, x_img, z_img, roi, xFull, zFull);
    axis normal;
    hold on;
    contour(xFull, zFull, roi, 2, 'w--')
    hold off;
    ylim(range_depth)

    xlabel('Lateral [°]'), ylabel('Depth [cm]');
   
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title(sprintf('$%s: \\alpha_s = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.alpha(:), 'omitnan'), ...
        std(maps_results_dof{idx}.alpha(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.alpha))), ...
        'Interpreter', 'latex');
    set(gca, 'fontsize', fontSize);
end
exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_a_rect.png"), ...
   'Resolution','300')

% %%%%%%%%%%%%%%%%%%%%%%%% b %%%%%%%%%%%%%%%%%%%%%%%%
indices_b = [1, 2, 4];  

figure;
tiledlayout(1,3, 'TileSpacing', 'tight')
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 600]); % [x, y, width, height]
for i = 1:length(indices_b)
    idx = indices_b(i);
    nexttile;

    [~,~,hColor] = imOverlayInterp(bmodeFull, maps_results_dof{idx}.b_dB, range_bmode, [], ...
                                   transparency, x_img, z_img, roi, xFull, zFull);
    axis normal;
    hold on;
    contour(xFull, zFull, roi, 2, 'w--')
    hold off;
    ylim(range_depth)

    xlabel('Lateral [°]'), ylabel('Depth [cm]');
   
    hColor.Label.String = 'dB';
    title(sprintf('$%s: \\Delta b = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.b_dB(:), 'omitnan'), ...
        std(maps_results_dof{idx}.b_dB(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.b_dB))), ...
        'Interpreter', 'latex');
    set(gca, 'fontsize', fontSize);
end
exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_b_rect.png"), ...
   'Resolution','300')

% %%%%%%%%%%%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%%%%%%%%

indices_n = [1, 2, 3];  

figure;
tiledlayout(1,3, 'TileSpacing', 'tight')
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 600]); % [x, y, width, height]
for i = 1:length(indices_n)
    idx = indices_n(i);
    nexttile;

    [~,~,hColor] = imOverlayInterp(bmodeFull, maps_results_dof{idx}.n, range_bmode, [], ...
                                   transparency, x_img, z_img, roi, xFull, zFull);
    axis normal;
    hold on;
    contour(xFull, zFull, roi, 2, 'w--')
    hold off;
    ylim(range_depth)

    xlabel('Lateral [°]'), ylabel('Depth [cm]');
   
    hColor.Label.String = 'a.u.';
    title(sprintf('$%s: \\Delta n = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.n(:), 'omitnan'), ...
        std(maps_results_dof{idx}.n(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.n))), ...
        'Interpreter', 'latex');
    set(gca, 'fontsize', fontSize);
end
exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_n_rect.png"), ...
    'Resolution','300')

%%  DISPLAY OVERLAY SIMPLE POLAR

% Plot in cartesian cords
[TH_acs,R_acs]          = meshgrid(-x_img*pi/180 + pi/2,z_img/100 + r0);
[xPolarACS,zPolarACS]   = pol2cart(TH_acs,R_acs);
zPolarACS               = zPolarACS + z0Polar;

%%%%%%%%%%%%%%%%%%%%%%% BMODE %%%%%%%%%%%%%%%%%%%%%%%
figure('Units', 'pixels', 'Position', [50, 100, 400, 300]); % [x, y, width, height]
idx = indices_alpha(i);

[ax1,~] = imOverlayPolar(BmodeFull,ones(size(maps_results_dof{1}.alpha)),range_bmode,range_acs,0, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')

title(ax1, samName, 'Interpreter', 'none');

xlim([-9 9])
ylim([0 15])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off

set(gca, 'fontsize', fontSize);

% exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_pol"+num2str(idx)+".png"), ...
% 'Resolution','300')
%% QUS PARAMETERS

%%%%%%%%%%%%%%%%%%%%%%%% alpha %%%%%%%%%%%%%%%%%%%%%%%%
indices_alpha = [1, 3, 4];  % Corresponden a 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

for i = 1:length(indices_alpha)
    figure('Units', 'pixels', 'Position', [50, 100, 400, 300]); % [x, y, width, height]
    idx = indices_alpha(i);

    % [ax1,~] = imOverlayPolar(BmodeFull,ones(size(maps_results_dof{1}.alpha)),range_bmode,range_acs,0, ...
    %     xPolar,zPolar,xPolarACS,zPolarACS);
    [ax1,~] = imOverlayPolar(BmodeFull,maps_results_dof{idx}.alpha,range_bmode,range_acs,transparency, ...
        xPolar,zPolar,xPolarACS,zPolarACS);
    % yticks(ax1,[4 8 12 16])
    xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')

    title(ax1, sprintf('$%s: \\alpha_s = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.alpha(:), 'omitnan'), ...
        std(maps_results_dof{idx}.alpha(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.alpha))), ...
        'Interpreter', 'latex');

    xlim([-9 9])
    ylim([0 15])
    hold on
    contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
    hold off
    % colorbar
   
    % hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

    set(gca, 'fontsize', fontSize);
    
    exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_a_pol"+num2str(idx)+".png"), ...
    'Resolution','300')
end

%% %%%%%%%%%%%%%%%%%%%%%% b %%%%%%%%%%%%%%%%%%%%%%%%
indices_b = [1, 2, 4];  % Corresponden a 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

for i = 1:length(indices_b)
    figure('Units', 'pixels', 'Position', [50, 100, 400, 300]); % [x, y, width, height]
    idx = indices_b(i);

    % [ax1,~] = imOverlayPolar(BmodeFull,ones(size(maps_results_dof{1}.alpha)),range_bmode,range_acs,0, ...
    %     xPolar,zPolar,xPolarACS,zPolarACS);
    [ax1,~] = imOverlayPolar(BmodeFull,maps_results_dof{idx}.b_dB,range_bmode,[ ],transparency, ...
        xPolar,zPolar,xPolarACS,zPolarACS);
    % yticks(ax1,[4 8 12 16])
    xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
    
    title(ax1, sprintf('$%s: \\Delta b = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.b_dB(:), 'omitnan'), ...
        std(maps_results_dof{idx}.b_dB(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.b_dB))), ...
        'Interpreter', 'latex');
    
    xlim([-9 9])
    ylim([0 15])
    hold on
    contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
    hold off
    % colorbar 
    % hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

    set(gca, 'fontsize', fontSize);

     exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_b_pol"+num2str(idx)+".png"), ...
    'Resolution','300')
end

%% %%%%%%%%%%%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%%%%%%%%
indices_b = [1, 2, 3];  % Corresponden a 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

for i = 1:length(indices_b)
    figure('Units', 'pixels', 'Position', [50, 100, 400, 300]); % [x, y, width, height]
    idx = indices_b(i);

    % [ax1,~] = imOverlayPolar(BmodeFull,ones(size(maps_results_dof{1}.alpha)),range_bmode,range_acs,0, ...
    %     xPolar,zPolar,xPolarACS,zPolarACS);
    [ax1,~] = imOverlayPolar(BmodeFull,maps_results_dof{idx}.n,range_bmode,[ ],transparency, ...
        xPolar,zPolar,xPolarACS,zPolarACS);
    % yticks(ax1,[4 8 12 16])
    xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')

    title(ax1, sprintf('$%s: \\Delta n = %.2f \\pm %.2f, \\%%CV = %.2f$', ...
        method_labels{idx}, ...
        mean(maps_results_dof{idx}.n(:), 'omitnan'), ...
        std(maps_results_dof{idx}.n(:), 'omitnan'), ...
        abs(cv2d(maps_results_dof{idx}.n))), ...
        'Interpreter', 'latex');

    xlim([-9 9])
    ylim([0 15])
    hold on
    contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
    hold off
    % colorbar 
    % hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';

    set(gca, 'fontsize', fontSize);

    exportgraphics(gcf,fullfile(figsDir,"sam_"+samName+"_n_pol"+num2str(idx)+".png"), ...
    'Resolution','300')
end

end
%% SAVING DATA
if roi_already
    % Saving ACS maps and ROI
    save(fullfile(resultsDir,samName+".mat"), ...
        'maps_results_dof','bsc_results_dof','rect')
else
    % first time working with roi
    save(fullfile(roisDir,samName+".mat"), ...
    'rect')

end

close all
pause (0.5);
end
%% Auxiliary functions
% Get delays
function [t_delay] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)

t_delay = zeros(length(t), n_elements, n_pulses);
% (x, z) [m] Obtain positions of center of every element
element_pos_x = Trans.ElementPos(:, 1)*wvl;
element_pos_z = Trans.ElementPos(:, 3)*wvl;
phi = Trans.ElementPos(:, 4);

for n = 1:n_pulses
    for e = 1:n_elements
        focus = sound_speed*t(:)/2;
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + ...
            (xfocus - element_pos_x(e)).^2))/sound_speed;
    end
end

end