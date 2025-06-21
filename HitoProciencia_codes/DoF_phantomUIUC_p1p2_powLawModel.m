% General Description: 
% *** Script for DoF test in Power Law Simulation using Regularized Power Law TV (by H. Chahuara)
% *** All results (metrics Table and figures can be saved in "dirOutcomes"
% Data Description:
% *** RF data from UIUC acq by SonixOne p1-p2 from 10.1109/TUFFC.2023.3245988
% *** theoretical BSC by Faran Theory 
% *** (check function duke_bsc_faranTheory, this requires faran_coeffs.mat)
% *** Data is saved in NAS2 'Q:\emiranda\data4Prociencia\phantoms'
% Requirements
% *** utilsPowerLaw, utilsRPLTV, utilsBSC, utilsUS, utilsMetrics, init.m
% Author:
% *** EAMZ based on LIM codes
% Date: W3M06Y25

%%%%%% MANUAL FIT FROM REFERENCE PHANTOM METHOD
% -----RPM PowLaw (b.(f.^n))-----
% Δb           = 15.294266
% b_s/b_r [dB] = 11.845286
% Δn           = -1.175598
% ---------
% -----RPM Gauss (g.exp(-s.f^2))-----
% d_g          = 6.016284
% g_s/g_r [dB] = 7.793283
% Δs           = 0.032991

%% INITIALIZATION
init
warning('off');

%% PARAMETERS

% Important constants
Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB;
range_bmode     = [-80 0];

% Options
saveOutcomes    = true; % save data

methodsRegu     = true;  % regularization activater RPL
plotBmode       = false; % plot bmode sample and reference ROI
plotBmodeOverl  = true;  % plot bmode overlay with colorimage
plotBSCdB       = true;  % plot \Delta b in dB
plotMaps        = false; % plot maps using imagesc
bsc_gt_by_rpm   = true;  % calculate GT by Reference Phantom Method
manualroi       = false;

% Directory outcomes (can be changed)
if saveOutcomes % (FIGURES & METRICS)
    dirOutcomes = './out/PROCIENCIA_Hito/phantomUIUC/';
    if (~exist(dirOutcomes)); mkdir (dirOutcomes); end
end

% Data Directories (can be changed)
dirData         = 'Q:\emiranda\data4Prociencia\phantoms'; % NAS2
folderDataSam   = 'UIUC';
folderDataRef   = 'UIUC';

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

% This label for visualization and standarizing the name
label_methods = {'3-DoF', '2-DoF_{b,n}', '2-DoF_{n,a}', '2-DoF_{b,a}'}; 

% First row for headers, second for data
bsc_results  = cell(2, length(methods)); 
maps_results = cell(2, length(methods));

bsc_results(1, :)  = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'),  sprintf('2-DoF "n"')};
maps_results(1, :) = {sprintf('3-DoF'), sprintf('2-DoF_{b,n}'), sprintf('2-DoF_{n,a}'),  sprintf('2-DoF_{b,a}')};

%% LOAD DATA
numPhantomSam = '1'; 
numPhantomRef = '2';

j_sam = 1.0;
j_ref = 1.0;

% SAMPLE SPECS
switch str2double(numPhantomSam)
    case 1
        alpha_sam   = 0.4; % [dB/cm/MHz]
    case 2
        alpha_sam   = 0.1; % [dB/cm/MHz]
end
% REFERENCE SPECS
switch str2double(numPhantomRef)
    case 1
        alpha_ref   = 0.4; % [dB/cm/MHz]
    case 2
        alpha_ref   = 0.1; % [dB/cm/MHz]
end

delta_alpha_prior   = alpha_sam - alpha_ref; % [dB/cm/MHz]
% delta_b_prior       = log(db2pow(11.845286)); % from RPM method
delta_b_prior       = log(db2pow(15.1)); % from RPM method the best
delta_n_prior       = -1.175598; % from RPM method

% STABLE
samName = ['rf_phantom', numPhantomSam];
SAM     = load (fullfile(dirData, folderDataSam, samName));
% SAM.rf  = SAM.rf(:,:,end-9:end); % Last frames 
SAM.acs = alpha_sam;

refName = ['rf_phantom', numPhantomRef];
REF     = load (fullfile(dirData, folderDataRef, refName));
% REF.rf  = REF.rf(:,:,end-9:end); % Last frames 
REF.acs = alpha_ref;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf(:,:,1)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf(:,:,1)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% ROI SELECTION

if ~manualroi

    pars.x_roi       = [-15 15]*1E-3; % [m] 
    pars.z_roi       = [5 35]*1E-3; % [m] % Same than sonix
    
else 

    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(SAM.x*1E3, SAM.z*1E3,bmode_sam,range_bmode);
    colormap gray; clim(range_bmode);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('Lateral [mm]'), ylabel('Depth [mm]'); 
    title('Bmode')
    
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

    pars.x_roi     = [rect(1), rect(1)+rect(3)]*1E-3; % [m]
    pars.z_roi     = [rect(2), rect(2)+rect(4)]*1E-3; % [m]
end

%% SPECTRAL METHOD PARAMETERS

pars.P           = 512; % NFFT for 10wl is 256, 20wl 512
pars.bw          = [2 7.5]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 10; % wavelengths
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

if (plotBmode)
deadBand = 0.1e-2;
figure,

subplot(121), 
imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image"), hold on;
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title(['SAM', ' P', numPhantomSam])
% ylim([deadBand*1000 50])
colormap('gray')

subplot(122), 
imagesc(REF.x*1E3, REF.z*1E3, bmode_ref, range_bmode), axis("image");
rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
clim(range_bmode)
cb = colorbar;
cb.Label.String = 'dB'; % Add the label "dB"
xlabel('Lateral [mm]'), ylabel('Depth [mm]');
title(['REF', ' P', numPhantomRef])
colormap('gray')
% ylim([deadBand*1000 50])
end

%% POWER SPECTRA ESTIMATION
% spectralData_sam = calc_powerSpectra(SAM, pars);
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars); % @
S_sam = spectralData_sam.powerSpectra;

num_ref = 1;

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
fprintf('d_g          = %f\n', d_g);
fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
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
par_rpl.df_op      = 0;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

%% FOR BUCLE

for iMet = 1:length(methods)

estim_method = methods{iMet};

%% COMPENSATE GAUSS ATTEMPT 2-DoF-a
if strcmp(estim_method, '2-DoF-a')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
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
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 

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

%% COMPENSATE GAUSS ATTEMPT 2-DoF-n
elseif strcmp(estim_method, '2-DoF-n')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E5]; % [mu_b, mu_n, mu_a]
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
W = kron( speye(p*q), -4*f );

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

% Prior "s"
s = +delta_n_prior*ones(p*q, 1);

% utils 
z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
dz = reshape(Dy*z(:),p,q);
dz(end,:) = dz(end-1,:);  

a_Np2dB = Np2dB*Dy*a./dz(:);

%% COMPENSATE GAUSS ATTEMPT 2-DoF-b
elseif strcmp(estim_method, '2-DoF-b')

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_b, mu_n, mu_a]
    else              mu_rpl_tv    = [0.001 0.001 0.001];
    end

band    = spectralData_sam.band;
depth   = spectralData_sam.depth;
[r,p,q] = size(SR);
  
comp_ref    = comp_ref_b_bsc(delta_b_prior);
comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

SR_comp = SR .* comp_ref .*comp_freq_a;

% indices initialization
f = band(:); % always column vector
% [r,p,q] = size(SR_comp);

% log-spectrum Ratio Y_g = Z.s + W.a 
Y = log(SR_comp);

% matrices for RPL-based algorithms
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 
W = kron( speye(p*q), -4*f );

% initialization for RPL-based methods
u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
    
% RPL estimation
[u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);

dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
Dy = sparse(kron(speye(q),dy));

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

    if (methodsRegu); mu_rpl_tv    = [1E3; 1E3; 10^3.5]; % [mu_b, mu_n, mu_a]
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
Z = kron( speye(p*q), log(f) ); % EMZ PowLaw  Size: [p*q*r, p*q] 
% Z = kron( speye(p*q), -f.^2 ); % EMZ Gauss  Size: [p*q*r, p*q] 
% W = kron( speye(p*q), -4*f.^j_sam );
W = kron( speye(p*q), -4*f );

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

%% FIGURE SIMPLE OVERLAY BMODE, ACS, DELTA BSC and DELTA N

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
colorImg        = bigImg(acs_sam, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [0.1 1.2];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);   
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
    title('$\alpha_s$', 'Interpreter', 'latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta alpha (ACS) %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta b (BSC) %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = bigImg(g_ratio, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

    if plotBSCdB 
       colorImg = bigImg(g_ratio_dB, spectralData_sam.rf_roi);
    end

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = '';
        if plotBSCdB 
            hColor.Label.String ='dB';
        end
    % title('$\frac{b_s}{b_r}$', 'Interpreter','latex')
    title('$\Delta b$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta g ratio in dB %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% Delta n %%%%%%%%%%%%%%%%%%%%%%%%%%

units           = 1E2;
bmodeFull       = bmode_sam;
colorImg        = bigImg(s_ratio, spectralData_sam.rf_roi);
range_bmode     = [-60 0];
range_img       = [];
transparency    = 0.65;
x_img           = spectralData_sam.x_roi*units;
z_img           = spectralData_sam.z_roi*units;
xFull           = SAM.x*units;
zFull           = SAM.z*units;

t = nexttile;
    [~,hB,hColor] = imOverlaySimple(bmodeFull, colorImg, range_bmode, range_img, ...
                        transparency, x_img, z_img, xFull, zFull);
    hold on;
    rectangle('Position', units*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), 
    hold off;
    xlabel('Lateral [cm]'), ylabel('Depth [cm]');
    hColor.Label.String = 'a.u.';
    title('$\Delta n$', 'Interpreter','latex')
    set(gca,'fontsize',fontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%% Delta s %%%%%%%%%%%%%%%%%%%%%%%%%%
end


%% SAVE ALL BSC EST GAUSS

% Delta_b*b_ref*f^(delta_n_prior +n_ref)

freq = spectralData_sam.band; % Given choosen BW

% OPTION B
b_sam_est      = median(g_ratio, 'all');
n_sam_est      = median(s_ratio, 'all');
bsc_est_powlaw = b_sam_est *(freq.^n_sam_est);

bsc_results{2, iMet} = bsc_est_powlaw;

% SAVE ALL MAPS
maps_results{2, iMet} = acs_sam; 
maps_results{3, iMet} = g_ratio_dB; 
maps_results{4, iMet} = s_ratio; 

end

%% BOX PLOT distribution a, delta b and delta n

font_size = 26;

delta_b_theo = exp(delta_b_prior);
delta_n_theo = delta_n_prior;
numMethods = size(maps_results, 2); % Number of methods (iMet values)

% Extract Data
acs_data = cell(1, numMethods);
g_ratio_data = cell(1, numMethods);
s_ratio_data = cell(1, numMethods);

for iMet = 1:numMethods
    acs_data{iMet}     = maps_results{2, iMet}(:);  % Flatten to column
    g_ratio_data{iMet} = maps_results{3, iMet}(:);
    s_ratio_data{iMet} = maps_results{4, iMet}(:);
end

% Convert to matrix for plotting (ensuring correct format)
acs_mat     = padconcatenation(acs_data, NaN, 1); % Pad with NaN for different lengths
g_ratio_mat = padconcatenation(g_ratio_data, NaN, 1);
s_ratio_mat = padconcatenation(s_ratio_data, NaN, 1);

% method_labels = string({maps_results{1, :}}); % Convert first row to string array
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

%%%%%%%%%%%%%%% Box Plot a %%%%%%%%%%%%%%%

% Exclude the second column in plot a
acs_mat_filtered = acs_mat(:, [1, 3, 4]);
method_labels_a = method_labels([1, 3, 4]);


figure; 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(acs_mat_filtered, 'Labels', method_labels_a);
% axis("image")
yline(alpha_sam, 'k--')

ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
xt = get(ax, 'XTick');
for i = 1:length(method_labels_a)
    text(xt(i), ax.YLim(1)-0.07*diff(ax.YLim), ['$' method_labels_a{i} '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end

% if (methodsRegu); ylim([-0.19 1.81])
% else              ylim([-80 80])
% end
title('\alpha');
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
set(gca, 'FontSize', font_size, 'FontWeight','normal');
%%%%%%%%%%%%%%% Box Plot a %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Box Plot b %%%%%%%%%%%%%%%
% Exclude the third column in plot b
g_ratio_mat_filtered = g_ratio_mat(:, [1, 2, 4]);
method_labels_b = method_labels([1, 2, 4]);

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(g_ratio_mat_filtered, 'Labels', method_labels_b);
% axis("image")
yline(10*log10(delta_b_theo), 'k--')
ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
xt = get(ax, 'XTick');
for i = 1:length(method_labels_b)
    text(xt(i), ax.YLim(1)-0.175*diff(ax.YLim), ['$', method_labels_b{i}, '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end
ylim([5 16])

title('\Deltab');
ylabel('\Deltab [dB]');
set(gca, 'FontSize', font_size);
%%%%%%%%%%%%%%% Box Plot b %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Box Plot n %%%%%%%%%%%%%%%

% Exclude the fourth column in plot n
s_ratio_mat_filtered = s_ratio_mat(:, [1, 2, 3]);
method_labels_n = method_labels([1, 2, 3]);

figure;
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
box on;
boxplot(s_ratio_mat_filtered, 'Labels', method_labels_n);
% axis("image")
yline(delta_n_theo, 'k--')
ax = gca;
ax.XTickLabel = {''}; % Remove default labels
ax.XTickLabelMode = 'manual';
xt = get(ax, 'XTick');
for i = 1:length(method_labels_n)
    text(xt(i), ax.YLim(1)-0.07*diff(ax.YLim), ['$' method_labels_n{i} '$'], ...
        'Interpreter','latex', 'HorizontalAlignment','center', 'FontSize', font_size, 'FontWeight','bold')
end
% if (methodsRegu); ylim([-0.25 0.05]); %yticks([0:0.01:0.05])
% else              ylim([-15 15])
% end
title('\Deltan');
ylabel('\Deltan [a.u.]');
set(gca, 'FontSize', font_size);
%%%%%%%%%%%%%%% Box Plot n %%%%%%%%%%%%%%%

function M = padconcatenation(C, padval, dim)
    % C: Cell array to concatenate
    % padval: Value used to pad (e.g., NaN)
    % dim: Dimension along which to concatenate (1 = rows, 2 = columns)
    max_length = max(cellfun(@numel, C));
    M = cellfun(@(x) padarray(x, [max_length - numel(x), 0], padval, 'post'), C, 'UniformOutput', false);
    M = cell2mat(M);
end

%% METRICS MAPS TABLE FORM 

clear m_3dof m_2dofa m_2dofb m_2dofn MetricsParam

%%%%%%%%%%%%%%%%%%% Metricas a %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, true(size(maps_results{2, 1})), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, true(size(maps_results{2, 1})), NaN, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{2, 3}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{2, 4}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-n');
m_3dof.param  = 'a';
m_2dofa.param = 'a';
m_2dofb.param = 'a';
m_2dofn.param = 'a';

MetricsParam(1:4) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas a %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metricas b %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, true(size(maps_results{3, 1})), pow2db(delta_b_theo), '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, true(size(maps_results{3, 2})), pow2db(delta_b_theo), '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{3, 3}, true(size(maps_results{3, 3})), NaN, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{3, 4}, true(size(maps_results{3, 4})), pow2db(delta_b_theo), '2-DoF-n');
m_3dof.param  = 'b';
m_2dofa.param = 'b';
m_2dofb.param = 'b';
m_2dofn.param = 'b';

MetricsParam(5:8) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas b %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metricas n %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, true(size(maps_results{4, 1})), delta_n_theo, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, true(size(maps_results{4, 2})), delta_n_theo, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{4, 3}, true(size(maps_results{4, 3})), delta_n_theo, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{4, 4}, true(size(maps_results{4, 4})), NaN, '2-DoF-n');
m_3dof.param  = 'n';
m_2dofa.param = 'n';
m_2dofb.param = 'n';
m_2dofn.param = 'n';

MetricsParam(9:12) = [m_3dof; m_2dofa; m_2dofb; m_2dofn];
%%%%%%%%%%%%%%%%%%% Metricas n %%%%%%%%%%%%%%%%%%%

T_combined        = struct2table(MetricsParam);
T_combined.method = categorical(T_combined.method);
T_combined.param  = categorical(T_combined.param);


%% METRICS BSC 

% BSC GT
BSC_gt          = bsc_rpm;
BSC_gt_interp   = interp1(BSC.band, BSC_gt, spectralData_ref.band, 'linear', 'extrap');
BSC_gt          = BSC_gt_interp;
freq            = spectralData_ref.band;

diff_fit_dB = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );

clear m_3dof m_2dofa m_2dofb m_2dofn MetricsBSC

m_3dof  = get_metrics_homo_gt(bsc_results{2, 1}, true(size(bsc_results{2, 1})), BSC_gt, '3-DoF');
m_3dof.diff_dB = diff_fit_dB(bsc_results{2, 1}, BSC_gt);
m_3dof.param = 'BSC';

m_2dofa  = get_metrics_homo_gt(bsc_results{2, 2}, true(size(bsc_results{2, 2})), BSC_gt, '2-DoF-a');
m_2dofa.diff_dB = diff_fit_dB(bsc_results{2, 2}, BSC_gt);
m_2dofa.param = 'BSC';

m_2dofb  = get_metrics_homo_gt(bsc_results{2, 3}, true(size(bsc_results{2, 3})), BSC_gt, '2-DoF-b');
m_2dofb.diff_dB = diff_fit_dB(bsc_results{2, 3}, BSC_gt);
m_2dofb.param = 'BSC';

m_2dofn  = get_metrics_homo_gt(bsc_results{2, 4}, true(size(bsc_results{2, 4})), BSC_gt, '2-DoF-n');
m_2dofn.diff_dB = diff_fit_dB(bsc_results{2, 4}, BSC_gt);
m_2dofn.param = 'BSC';

MetricsBSC(1:4) = [m_3dof; m_2dofa; m_2dofb; m_2dofn]; 

Tbsc        = struct2table(MetricsBSC);
Tbsc.method = categorical(Tbsc.method);
Tbsc.param  = categorical(Tbsc.param);

%%  PLOT DELTA BSC

% Define properties for customization
xlim_range = pars.bw + [-0.01 0.01];
ylim_range = [10^-1 10^1]; % Y-axis limits
line_width = 2.85; % Set line width
font_size = 16; % Adjust font size

bsc_results{1, 1} = sprintf('3-DoF     (GoF_{dB} = %.2f) \n', MetricsBSC(1).diff_dB);
bsc_results{1, 2} = sprintf('2-DoF_{b,n} (GoF_{dB} = %.2f) \n', MetricsBSC(2).diff_dB);
bsc_results{1, 3} = sprintf('2-DoF_{n,a} (GoF_{dB} = %.2f) \n', MetricsBSC(3).diff_dB);
bsc_results{1, 4} = sprintf('2-DoF_{b,a} (GoF_{dB} = %.2f) \n', MetricsBSC(4).diff_dB);

% Convert hexadecimal colors to RGB (MATLAB requires values between 0 and 1)
colot_gt = '#000000';  % Black
color_1  = '#FF0000';  % 3dof
color_2  = '#D95319';  % 2dof a
color_3  = '#0072BD';  % 2dof b
color_4  = '#77AC30';  % 2dof n

% Create figure and plot data
figure, 
set(gcf, 'Units', 'pixels', 'Position', [100, 100, 700, 700]); % [x, y, width, height] in pixels
semilogy(freq, BSC_gt, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
hold on;
semilogy(freq, bsc_results{2, 1}, '--', 'Color', hex2rgb(color_1), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 1});
semilogy(freq, bsc_results{2, 2}, '--', 'Color', hex2rgb(color_2), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 2});
semilogy(freq, bsc_results{2, 3}, '--', 'Color', hex2rgb(color_3), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 3});
semilogy(freq, bsc_results{2, 4}, '--', 'Color', hex2rgb(color_4), 'LineWidth', line_width, 'DisplayName', bsc_results{1, 4});
hold off;

% Customize plot
grid on;
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
% ylim(ylim_range);
xlim(xlim_range);
title('Backscatter Coefficient (BSC)', 'FontSize', font_size + 2);
legend('Location', 'best', 'FontSize', font_size + 2);
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

%% SAVE OUTCOMES

if (saveOutcomes)

% SAVE METRICS (EXCEL)

%%%%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%
if methodsRegu;   nameExcel = 'metricsRegu_p1p2UIUC.xlsx'; 
else     nameExcel = 'metricsNoRegu_p1p2UIUC.xlsx'; 
end

% Define output file name
excelFile = fullfile(dirOutcomes, nameExcel);

% Write to Excel
writetable(T_combined, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

fprintf('Table Maps Metrics saved to %s\n', excelFile);

%%%%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% BSC PLOT %%%%%%%%%%%%%%%%%%%%
% Write to Excel
if methodsRegu;   nameExcel = 'metricsBSCregu_p1p2UIUC.xlsx'; 
else     nameExcel = 'metricsBSCnoregu_p1p2UIUC.xlsx'; 
end

excelFile = fullfile(dirOutcomes, nameExcel);

writetable(Tbsc, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

fprintf('Table BSC saved to %s\n', excelFile);

% SAVE FIGURES (PNG)

titleFigout = 'UIUCp1p2_Fig';
save_all_figures_to_directory(dirOutcomes, titleFigout) % 
% save_all_figures_to_directory(dirOutcomes, titleFigout, 'svg') % 
fprintf('Figures saved %s\n', excelFile);

end