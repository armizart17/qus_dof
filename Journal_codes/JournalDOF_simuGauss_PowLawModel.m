% General Description: 
% *** Script for DoF test in Gauss Simulation using Regularized Power Law TV (with Gaussian Model)
% *** All results (metrics Table and figures can be saved in "dirOutcomes"
% Data Description:
% *** Data simulated in k-wave and coloured power spectra in post-processing
% *** and is read from directory "dirData" and folders "folderDataSam" "folderDataRef"
% *** Data is saved in NAS2 'Q:\emiranda\data4Prociencia\simus'
% Requirements
% *** utilsPowerLaw, utilsRPLTV, utilsBSC, utilsUS, utilsMetrics, init.m
% Author:
% *** EAMZ based on LIM codes
% Date: W3M06Y25
%
%%%%% FIT RPM %%%%%
% BW = 3.75 - 8.25 MHz, Fc = 6 MHz
% wavelength (wl) = 0.256667 mm
% Exec Time: 47.5765 
% -----RPM PowLaw (b.(f.^n))-----
% Δb           = 4.605340
% Δb [dB]      = 6.632617
% Δn           = -1.280963
% ---------
% -----RPM Gauss (g.exp(-s.f^2))-----
% Δg           = 0.967725
% Δg      [dB] = -0.142481
% Δs           = 0.018690
% ---------
% BW = 3.7 - 8.2 MHz, Fc = 5.95 MHz
% wavelength (wl) = 0.258824 mm
% Exec Time: 41.5274 
% -----RPM PowLaw (b.(f.^n))-----
% Δb           = 4.381608
% b_s/b_r [dB] = 6.416335
% Δn           = -1.255213
% ---------
% -----RPM Gauss (g.exp(-s.f^2))-----
% Δg           = 0.964408
% Δg      [dB] = -0.157391
% Δs           = 0.018717
% ---------
%% INITIALIZATION
init
warning('off');
%% PARAMETERS

% Important constants
Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB;
range_bmode     = [-80 0];

% Options
saveOutcomes    = true;  % save data

methodsRegu     = true;  % regularization activater RPL
plotBmode       = false; % plot bmode sample and reference ROI
plotBmodeOverl  = false;  % plot bmode overlay with colorimage
plotBSCdB       = true;  % plot \Delta b in dB
plotMaps        = false; % plot maps using imagesc
bsc_gt_by_rpm   = true;  % calculate GT by Reference Phantom Method

% Directory outcomes (can be changed)
if saveOutcomes % (FIGURES & METRICS)
    dirOutcomes = './out/dofJournal/simuGauss/powLawModel';
    if (~exist(dirOutcomes)); mkdir (dirOutcomes); end
end

% Data Directories (can be changed)
dirData         = 'Q:\emiranda\data4Prociencia\simus';
folderDataSam   = 'gaussSimu';
folderDataRef   = 'gaussSimu';

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

% This label for visualization and standarizing the name
label_methods = {'3-DoF', '2-DoF_{b,n}', '2-DoF_{n,a}', '2-DoF_{b,a}'}; 

% First row for headers, second for data
bsc_results     = cell(2, length(methods)); 
maps_results    = cell(2, length(methods));

% Store headers
bsc_results(1, :)   = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :)  = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};

%% LOAD DATA

% SAM DATA SPECS
rf_sam_name = 'rf_sd2pcSCALE4_bsc2_gauss_mask7_att0p5'; 
alpha_sam   = 0.5;
j_sam       = 1.1;

SAM             = load(fullfile(dirData, folderDataSam, rf_sam_name));
SAM.acs         = alpha_sam; % [dB/cm/MHz] 
SAM.alpha_power = j_sam;


% REF DATA SPECS 
rf_ref_name = 'rf_sd2pcSCALE4_bsc4_att0p4';  % AC
alpha_ref   = 0.4;
j_ref       = 1.1;

REF             = load(fullfile(dirData, folderDataRef, rf_ref_name)); 
REF.acs         = alpha_ref; % [dB/cm/MHz]
REF.alpha_power = j_ref; 

% FROM RPM (GROUNTRUTH) ESTIMATE DELTA_S_PRIOR
delta_alpha_prior = alpha_sam - alpha_ref; % [dB/cm/MHz]
% BW 3.75 - 8.25 MHz
delta_b_prior     = log(db2pow(6.416335)); 
delta_n_prior     = -1.280963; % n_sam - n_ref

% GT
delta_b_theo = 4.605340;
delta_n_theo = delta_n_prior;

% B-MODE CHECK
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL METHOD PARAMETERS
pars.P           = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
% pars.bw          = [4 7]; % [MHz] % just test
% pars.bw          = [3.75 8.25]; % [MHz] % just test
pars.bw          = [3.7 8.2]; % [MHz] % last
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths
pars.z_roi       = [5 42.5]*1E-3; % all
pars.x_roi       = [-18 18]*1E-3;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

if (plotBmode)
figure,

subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image"), hold on;
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('SAM')
colormap('gray')

subplot(122), 
imagesc(REF.x*1E3, REF.z*1E3, bmode_ref, range_bmode), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title('REF')
colormap('gray')
end

%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

% spectralData_ref = calc_powerSpectra(REF, pars);
spectralData_ref = calc_powerSpectra_vSimple(REF, pars); % @
S_ref = spectralData_ref.powerSpectra;

SR_emz = S_sam ./ S_ref;

SR = permute(SR_emz,[3,1,2]); clear SR_emz
 
%% BSC RPM (GT)
if (bsc_gt_by_rpm)
BSC     = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
% bsc_rpm = BSC.BSCcurve_Uni(:,1); % mean
bsc_rpm = BSC.BSCcurve_Uni(:,2); % median
freq    = BSC.band;

%%%%%%%%%%%%%%%%% POWER LAW %%%%%%%%%%%%%%%%%
% Perform linear regression  ln(bsc) = d_n . ln(f) + ln(d_b) 
coeffs   = polyfit(log(freq), log(bsc_rpm), 1); % Fit y = mx + c
d_n      = coeffs(1); % Slope = d_n
ln_db    = coeffs(2); % Intercept = ln(d_b) 
d_b      = exp(ln_db); % 

% Display results
fprintf('-----RPM PowLaw (b.(f.^n))-----\n')
fprintf('Δb           = %f\n', d_b);
fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b));
fprintf('Δn           = %f\n', d_n);
fprintf('---------\n')

bsc_rpm_powlaw = d_b*(freq.^d_n);
%%%%%%%%%%%%%%%%% POWER LAW %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% GAUSSIAN %%%%%%%%%%%%%%%%%
% Perform linear regression  bsc = d_s . f^2 + ln(d_g) 
coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
ln_dg    = coeffs(2); % Intercept = ln(d_g) 
d_g      = exp(ln_dg); % 1.0917 (mean), 0.9079(median)
% Display results
fprintf('-----RPM Gauss (g.exp(-s.f^2))-----\n')
fprintf('Δg           = %f\n', d_g);
fprintf('Δg      [dB] = %f\n', 10*log10(d_g));
fprintf('Δs           = %f\n', d_s);
fprintf('---------\n')
bsc_rpm_gauss = d_g*exp(-d_s* freq.^2);
%%%%%%%%%%%%%%%%% GAUSSIAN %%%%%%%%%%%%%%%%%
end
%% GENERAL REGULARIZTATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

%% FOR BUCLE
for iMet = 1:length(methods)

estim_method = methods{iMet};

%% COMPENSATE GAUSS ATTEMPT 2-DoF-a
if strcmp(estim_method, '2-DoF-a')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; NaN]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .* comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_a = X.g + Z.s 
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
% Z = kron( speye(p*q), -f.^2 ); % EMZ Gaussian
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  


% initialization for RPL-based methods
u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
s = u_opt(p*q+1:2*p*q);

% Prior "a" known
a_Np2dB = delta_alpha_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

%% COMPENSATE 2-DoF-n
elseif strcmp( estim_method, '2-DoF-n')
    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 10^4]; % [mu_b, mu_n, mu_a]
    % if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_n_bsc(delta_n_prior, band, p, q);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .* comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_s = X.g + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
% W = kron( speye(p*q), -4*f );
W = kron( speye(p*q), -4*f.^j_sam );

% initialization for RPL-based methods
u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

% Prior "n"
s = +delta_n_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);  

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE 2-DoF-b
elseif strcmp(estim_method, '2-DoF-b')
    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 10^4]; % [mu_b, mu_n, mu_a]
    % if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_g_bsc(delta_b_prior);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_g = Z.s + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
% Z = kron( speye(p*q), -f.^2 ); % EMZ GAUSSIAN
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 
% W = kron( speye(p*q), -4*f );
W = kron( speye(p*q), -4*f.^j_sam );

% initialization for RPL-based methods
u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);

if par_rpl.df_op == 1
dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));
else
dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
Dy = sparse(kron(speye(q),dy)); %diag(ones(p*q -1,1),1) - diag(ones(p*q,1));    
end

% \Deltas
s = u_opt(1:p*q);
a = u_opt(p*q+1:2*p*q);

% Prior "g"
g = delta_b_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);       

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 3-DoF
elseif strcmp(estim_method, '3-DoF')
    
    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 10^3.8]; % [mu_b, mu_n, mu_a] * old
    % if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 10^3.85]; % [mu_b, mu_n, mu_a] *new tuffc
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);

comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR.*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y = X.g + Z.s + W.a
Y = log(SR_comp);

% matrices for RPL-based algorithms
X = kron( speye(p*q), ones(size(f)) );
% Z = kron( speye(p*q), -f.^2 ); % EMZ Gaussian
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  
W = kron( speye(p*q), -4*f.^j_sam );
% W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

% \Deltas
g = u_opt(1:p*q);
s = u_opt(p*q+1:2*p*q);
a = u_opt(2*p*q+1:3*p*q);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);

a_Np2dB = Np2dB*Dy*a./dz(:);

end

%% QUS PARAMETERS 

g_ratio     = reshape(exp(g), p, q);
g_ratio_dB  = 10*log10(g_ratio);
alpha_ratio = reshape(a_Np2dB, p, q);
s_ratio     = reshape(s, p, q); 

% METRICS 3 DOF
acs_sam   = alpha_ratio + alpha_ref;

mean2d = @(x) mean(x(:));
std2d = @(x) std(x(:));
cv2d = @(x) 100*std(x(:))/mean(x(:));

calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

[m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio), calc2dStats{2}(g_ratio), calc2dStats{3}(g_ratio));
if plotBSCdB 
    [m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio_dB), calc2dStats{2}(g_ratio_dB), calc2dStats{3}(g_ratio_dB));
end

[m_s, s_s, cv_s] = deal(calc2dStats{1}(s_ratio), calc2dStats{2}(s_ratio), calc2dStats{3}(s_ratio));
[m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam), calc2dStats{2}(acs_sam), calc2dStats{3}(acs_sam));

fprintf('-----%s---\n', estim_method);
fprintf('α_s        : %.3f ± %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
    if plotBSCdB 
fprintf('Δb [dB]    : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    else
fprintf('Δb         : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    end

fprintf('Δn         : %.4f ± %.4f, %%CV = %.4f\n', round(m_s, 4), round(s_s, 4), round(cv_s, 4));
fprintf('--------\n');
 
%% IMAGESC PLOTS
Xaxis = spectralData_ref.lateral;
Zaxis = spectralData_ref.depth;
cm = 1e2;

axis_n = [0 1.2];
axis_a = [0 3];
axis_b = [-60 0]; % dB
fontSize = 16;

if plotMaps
figure, 
set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
sgtitle(label_methods{iMet}, 'FontSize', fontSize+2, 'FontWeight', 'bold');

subplot(1,3,1)
imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
axis("image");
xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
% title('\Delta \alpha ');
% title(['ACS: ', num2str(round(m_a, 3)), ' \pm ', num2str(round(s_a, 2)), ', CV = ', num2str(round(cv_a, 3))])
title({'$\alpha_s$:', ...
       [num2str(round(m_a, 3)), ' $\pm$ ', num2str(round(s_a, 3)), ', CV = ', num2str(round(cv_a, 3))]}, ...
      'Interpreter', 'latex');

h2 = colorbar; 
ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}','FontSize', fontSize);
set(gca,'fontsize',fontSize)

subplot(1,3,2)
imagesc(Xaxis*cm, Zaxis*cm, g_ratio)
h2 = colorbar; 
if plotBSCdB 
   imagesc(Xaxis*cm, Zaxis*cm, g_ratio_dB)
   h2 = colorbar;
   ylabel(h2,'dB','FontSize', fontSize);
end
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['$\frac{g_s}{g_r}: ', num2str(round(m_g, 3)), ' \pm ', num2str(round(s_g, 2)), ', CV = ', num2str(round(cv_g, 3)), '$'], ...
%       'Interpreter', 'latex')
title({'$\Delta b$:', ...
       [num2str(round(m_g, 3)), ' $\pm$ ', num2str(round(s_g, 3)), ', CV = ', num2str(round(cv_g, 3))]}, ...
      'Interpreter', 'latex');

set(gca,'fontsize',fontSize)

subplot(1,3,3)
imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
axis("image");
xlabel('Lateral [cm]'), colormap turbo;
% title(['$\Delta s$: ', num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))], ...
%       'Interpreter', 'latex');
title({'$\Delta n$:', ...
       [num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))]}, ...
      'Interpreter', 'latex');
set(gca,'fontsize',fontSize)

end
%% FIGURE INTERP OVERLAY BMODE, DELTA SNR, ACS, DELTA BSC, DELTA N
if plotBmodeOverl
fontSize = 16;

figure,
set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;

tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
sgtitle(label_methods{iMet}, 'FontSize', fontSize+2, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%% alpha_s (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%
acs_sam = alpha_ratio + alpha_ref;

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = acs_sam;
% range_bmode     = [-100 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);   
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('$\alpha_s$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = g_ratio;
% range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));

    if plotBSCdB 
       colorImg = g_ratio_dB;
    end

t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='dB';
        end
    title('$\Delta b$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = s_ratio;
% range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.lateral*units;
z_img           = spectralData_sam.depth*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;
[X, Z] = meshgrid(xFull, zFull);
roi = and(X >= x_img(1), X <= x_img(end)) & ...
      and(Z >= z_img(1), Z <= z_img(end));


t = nexttile;
    [~,hB,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, roi, xFull, zFull);
    hold on;
    contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2)
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'a.u.';
    title('$\Delta n$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% BSC RECONSTRUCTION DATA AC (change later)
freq = spectralData_sam.band;

g_est = median(g_ratio(:));
s_est = median(s_ratio(:));

% bsc_est_gauss = g_est.*exp(-s_est.*freq.^2);
bsc_est_powlaw = g_est.*(freq.^s_est);

% SAVE ALL BSC EST POW LAW
bsc_results{2, iMet} = bsc_est_powlaw;

% SAVE ALL MAPS
maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 

end

%% BOX PLOT distribution a, delta b and delta n

font_size  = 34;  
numMethods = size(maps_results, 2); % Number of methods (iMet values)

% Extract Data
acs_data        = cell(1, numMethods);
b_ratio_data    = cell(1, numMethods);
n_ratio_data    = cell(1, numMethods);

for iMet = 1:numMethods
    acs_data{iMet}     = maps_results{2, iMet}(:);  % Flatten to column
    b_ratio_data{iMet} = maps_results{3, iMet}(:);
    n_ratio_data{iMet} = maps_results{4, iMet}(:);
end

% Convert to matrix for plotting (ensuring correct format)
acs_mat     = padconcatenation(acs_data, NaN, 1); % Pad with NaN for different lengths
b_ratio_mat = padconcatenation(b_ratio_data, NaN, 1);
n_ratio_mat = padconcatenation(n_ratio_data, NaN, 1);

% method_labels = string({maps_results{1, :}}); % Convert first row to string array
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

% Exclude the second column in plot a
acs_mat_filtered = acs_mat(:, [1, 3, 4]);
method_labels_a  = method_labels([1, 3, 4]);

% Box Plot a
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]); % [x, y, width, height] in pixels
box on;
boxplot(acs_mat_filtered, 'Labels', method_labels_a);
% axis("image")
yline(alpha_sam, 'k--')
ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
if (methodsRegu)  ylim([0.05 1.01])
else              ylim([-20 20])
end
xt = get(ax, 'XTick');
for i = 1:length(method_labels_a)
    text(xt(i), ax.YLim(1)*0.005, ['$' method_labels_a{i} '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end

title('\bf\alpha')
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
set(gca, 'FontSize', font_size);

%% Exclude the third column in plot b
b_ratio_mat_filtered = b_ratio_mat(:, [1, 2, 4]);
method_labels_b      = method_labels([1, 2, 4]);

% Box Plot b
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]); % [x, y, width, height] in pixels
box on;
boxplot(b_ratio_mat_filtered, 'Labels', method_labels_b);
% axis("image")
yline(10*log10(delta_b_theo), 'k--')
ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
xt = get(ax, 'XTick');
if (methodsRegu) ylim([-4 7.5])
else              ylim([-60 20])
end
for i = 1:length(method_labels_b)
    text(xt(i), 1.15*ax.YLim(1), ['$' method_labels_b{i} '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end

title('\bf\Deltab');
ylabel('\Deltab [dB]');
set(gca, 'FontSize', font_size);
%
% Exclude the fourth column in plot n
n_ratio_mat_filtered = n_ratio_mat(:, [1, 2, 3]);
method_labels_n      = method_labels([1, 2, 3]);

%% Box Plot n
figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]); % [x, y, width, height] in pixels
box on;
boxplot(n_ratio_mat_filtered, 'Labels', method_labels_n);
% axis("image")
yline(delta_n_theo, 'k--')
ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
xt = get(ax, 'XTick');
if (methodsRegu) % ylim([-0.01 0.05])
else              ylim([-15 15])
end
for i = 1:length(method_labels_n)
    
    % text(xt(i), ax.YLim(1)-0.6*diff(ax.YLim), ['$' method_labels_n{i} '$'], ...
    %     'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
text(xt(i), 1.15*ax.YLim(1), ['$' method_labels_n{i} '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end

title('\bf\Deltan');
ylabel('\Deltan [a.u.]');
set(gca, 'FontSize', font_size);

function M = padconcatenation(C, padval, dim)
    % C: Cell array to concatenate
    % padval: Value used to pad (e.g., NaN)
    % dim: Dimension along which to concatenate (1 = rows, 2 = columns)
    max_length = max(cellfun(@numel, C));
    M = cellfun(@(x) padarray(x, [max_length - numel(x), 0], padval, 'post'), C, 'UniformOutput', false);
    M = cell2mat(M);
end

%%  METRICS BSC

% BSC GT
BSC_gt          = bsc_rpm;

diff_fit_dB     = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );

clear m_3dof m_2dofa m_2dofb m_2dofn MetricsBSC
m_3dof          = get_metrics_homo_gt(bsc_results{2, 1}, true(size(bsc_results{2, 1})), BSC_gt, '3-DoF');
m_3dof.diff_dB  = diff_fit_dB(bsc_results{2, 1}, BSC_gt);
m_3dof.param    = 'BSC';

m_2dofa         = get_metrics_homo_gt(bsc_results{2, 2}, true(size(bsc_results{2, 2})), BSC_gt, '2-DoF-a');
m_2dofa.diff_dB = diff_fit_dB(bsc_results{2, 2}, BSC_gt);
m_2dofa.param   = 'BSC';

m_2dofb         = get_metrics_homo_gt(bsc_results{2, 3}, true(size(bsc_results{2, 3})), BSC_gt, '2-DoF-b');
m_2dofb.diff_dB = diff_fit_dB(bsc_results{2, 3}, BSC_gt);
m_2dofb.param   = 'BSC';

m_2dofn         = get_metrics_homo_gt(bsc_results{2, 4}, true(size(bsc_results{2, 4})), BSC_gt, '2-DoF-n');
m_2dofn.diff_dB = diff_fit_dB(bsc_results{2, 4}, BSC_gt);
m_2dofn.param   = 'BSC';

MetricsBSC(1:4) = [m_3dof; m_2dofa; m_2dofb; m_2dofn]; 

Tbsc        = struct2table(MetricsBSC);
Tbsc.method = categorical(Tbsc.method);
Tbsc.param  = categorical(Tbsc.param);

%%  PLOT DELTA BSC

% Define properties for customization
xlim_range = pars.bw+ 0.05*[-1 1]; % X-axis limits
ylim_range = [10^-0.65 10^1.1]; % Y-axis limits
line_width = 3.5; % Set line width
font_size  = 32; % Adjust font size

% GoF
bsc_results{1, 1} = sprintf('3-DoF     (GoF_{dB} = %.2f) \n', MetricsBSC(1).diff_dB);
bsc_results{1, 2} = sprintf('2-DoF_{b,n} (GoF_{dB} = %.2f) \n', MetricsBSC(2).diff_dB);
bsc_results{1, 3} = sprintf('2-DoF_{n,a} (GoF_{dB} = %.2f) \n', MetricsBSC(3).diff_dB);
bsc_results{1, 4} = sprintf('2-DoF_{b,a} (GoF_{dB} = %.2f) \n', MetricsBSC(4).diff_dB);

% Convert HEX colors to RGB (MATLAB requires values between 0 and 1)
colot_gt = '#000000';  % Black
color_1  = '#FF0000';  % 3dof
color_2  = '#D95319';  % 2dof a
color_3  = '#0072BD';  % 2dof b
color_4  = '#77AC30';  % 2dof n

% Create figure and plot data
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]); % [x, y, width, height] in pixels

% semilogy(freq, BSC_gt, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
semilogy(freq, bsc_results{2, 1}, '--', 'Color', hex2rgb(color_1), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 1});
hold on;
semilogy(freq, bsc_results{2, 2}, '--', 'Color', hex2rgb(color_2), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 2});
semilogy(freq, bsc_results{2, 3}, '--', 'Color', hex2rgb(color_3), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 3});
semilogy(freq, bsc_results{2, 4}, '--', 'Color', hex2rgb(color_4), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 4});

semilogy(freq, BSC_gt, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
hold off;

% Customize plot
grid on;
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('\DeltaBSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
ylim(ylim_range);
xlim(xlim_range);
% title('\DeltaBSC', 'FontSize', font_size + 2);
title('Power Law Model', 'FontSize', font_size + 2);
legend('Location', 'best', 'FontSize', font_size-7);
set(gca, 'FontSize', font_size);
hold off;

% Function to convert hexadecimal to RGB
function rgb = hex2rgb(hex)
    hex = char(hex); % Ensure it's a string
    if hex(1) == '#'
        hex(1) = []; % Remove the '#' symbol
    end
    rgb = reshape(sscanf(hex, '%2x') / 255, 1, 3);
end

%% METRICS MAPS TABLE FORM 

% Metricas a
m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, logical(ones(size(acs_sam))), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, logical(ones(size(acs_sam))), NaN, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{2, 3}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{2, 4}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-n');

% Extract field names
fields = fieldnames(m_3dof);

% Create a table
Ta = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

% Metrics b dB
m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{3, 3}, logical(ones(size(acs_sam))), NaN, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{3, 4}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '2-DoF-n');

Tb = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

% Metrics n
m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, logical(ones(size(acs_sam))), delta_n_theo, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, logical(ones(size(acs_sam))), delta_n_theo, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{4, 3}, logical(ones(size(acs_sam))), delta_n_theo, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{4, 4}, logical(ones(size(acs_sam))), NaN, '2-DoF-n');

Tn = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), 'RowNames', fields, 'VariableNames', methods);

clear m_3dof m_2dofa m_2dofb m_2dofn

% Add a new column to each table indicating its group
Ta.Group = repmat("a", height(Ta), 1);
Tb.Group = repmat("b", height(Tb), 1);
Tn.Group = repmat("n", height(Tn), 1);

% Reorder columns so "Group" is the first column
Ta = movevars(Ta, 'Group', 'Before', Ta.Properties.VariableNames{1});
Tb = movevars(Tb, 'Group', 'Before', Tb.Properties.VariableNames{1});
Tn = movevars(Tn, 'Group', 'Before', Tn.Properties.VariableNames{1});

%  Modify row names to avoid duplicates
Ta.Properties.RowNames = strcat(Ta.Properties.RowNames, " a");
Tb.Properties.RowNames = strcat(Tb.Properties.RowNames, " b");
Tn.Properties.RowNames = strcat(Tn.Properties.RowNames, " n");

T_combined = [Ta; Tb; Tn];

%% SAVE OUTCOMES

if (saveOutcomes)

% SAVE METRICS (EXCEL)

%%%%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%
if methodsRegu;   nameExcel = 'metricsRegu_gaussSimu.xlsx'; 
else     nameExcel = 'metricsNoRegu_gaussSimu.xlsx'; 
end

% Define output file name
excelFile = fullfile(dirOutcomes, nameExcel);

% Write to Excel
writetable(T_combined, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

fprintf('Table Maps Metrics saved to %s\n', excelFile);

%%%%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% BSC PLOT %%%%%%%%%%%%%%%%%%%%
% Write to Excel
if methodsRegu;   nameExcel = 'metricsBSCregu_gaussSimu.xlsx'; 
else     nameExcel = 'metricsBSCnoregu_gaussSimu.xlsx'; 
end

excelFile = fullfile(dirOutcomes, nameExcel);

writetable(Tbsc, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

fprintf('Table BSC saved to %s\n', excelFile);

% SAVE FIGURES (PNG)

titleFigout = 'simuG_Fig';
% save_all_figures_to_directory(dirOutcomes, titleFigout) % 
save_all_figures_to_directory(dirOutcomes, titleFigout, 'svg') % 
fprintf('Figures saved %s\n', excelFile);

end