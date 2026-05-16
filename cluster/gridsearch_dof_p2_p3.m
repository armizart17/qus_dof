function gridsearch_dof_p2_p3()
% GRIDSEARCH_DOF_POWLAW_SAVE_CASES
% Full grid search over mu_a, mu_b, mu_n for the DoF Power Law simulation.
% Saves each case separately.
%
% Requirements:
%   utilsPowerLaw, utilsRPLTV, utilsBSC, utilsUS, utilsMetrics, init.m
%   and all functions already used by your original script:
%   calc_powerSpectra_vSimple, comp_ref_a, comp_mod_freq_a, initialize_rpl,
%   rpl_tv, get_metrics_homo_gt, etc.
%
% Author:
%   Adapted from EAMZ script for per-case grid search saving (P1 - P3)

clear; clc; close all;
warning('off');

cd('../'); % when you run it usually goes inside ./parentFolder/script.m
disp(pwd)
addpath(genpath(pwd))
addpath(genpath('/opt/MATLAB Add-Ons'))

%% ========================= USER CONFIG =========================
% Choose the mu grids
list_mu_a = 10.^(-3:0.5:6);
list_mu_b = 10.^(-3:0.5:6);
list_mu_n = 10.^(-3:0.5:6);

% list_mu_a = 10^4.1;
% list_mu_b = [1000, 10]; % simple test
% list_mu_n = 10^3;


% Output directory
if ispc
    % Windows (local)
    dirOut = './out/dofJournal/p2p3/';
elseif isunix
    % Linux cluster
    dirOut = '/mnt/nfs2/emiranda/proj/dof26/out/p2p3';  % explicit
else
    error('Unknown operating system');
end


if ~exist(dirOut, 'dir')
    mkdir(dirOut);
end

timestamp = datestr(datetime('now'), 'yyyy-mm-dd_HH-MM');
dirCases = fullfile(dirOut, ['gridSearch_cases_' timestamp]);
if ~exist(dirCases, 'dir')
    mkdir(dirCases);
end

% Skip case if file already exists
skipExisting = true;

% Keep plots disabled for grid search
plotBmode       = false;
plotBmodeOverl  = false;
plotBSCdB       = true;
plotMaps        = false;
saveOutcomes    = false;
methodsRegu     = true;

% Save large map arrays too
saveMapsAndBSC  = true;

% ===============================================================

%% Build all combinations
[MU_A, MU_B, MU_N] = ndgrid(list_mu_a, list_mu_b, list_mu_n);

comb_mu_a = MU_A(:);
comb_mu_b = MU_B(:);
comb_mu_n = MU_N(:);

nComb = numel(comb_mu_a);
fprintf('Total combinations: %d\n', nComb);

fmtExp = @(x) strrep(sprintf('%.1f', log10(x)), '.', 'p');

%% Summary structure
% Summary(nComb,1) = struct( ...
%     'idx', [], ...
%     'mu_a', [], ...
%     'mu_b', [], ...
%     'mu_n', [], ...
%     'log10_mu_a', [], ...
%     'log10_mu_b', [], ...
%     'log10_mu_n', [], ...
%     'fileName', "", ...
%     'filePath', "", ...
%     'status', "", ...
%     'elapsed_sec', NaN, ...
%     'mpe_a_3dof', NaN, ...
%     'mpe_b_3dof', NaN, ...
%     'mpe_n_3dof', NaN, ...
%     'mpe_a_2dofa', NaN, ...
%     'mpe_b_2dofa', NaN, ...
%     'mpe_n_2dofa', NaN, ...
%     'mpe_a_2dofb', NaN, ...
%     'mpe_b_2dofb', NaN, ...
%     'mpe_n_2dofb', NaN, ...
%     'mpe_a_2dofn', NaN, ...
%     'mpe_b_2dofn', NaN, ...
%     'mpe_n_2dofn', NaN, ...
%     'gof_bsc_3dof', NaN, ...
%     'gof_bsc_2dofa', NaN, ...
%     'gof_bsc_2dofb', NaN, ...
%     'gof_bsc_2dofn', NaN);
% 
% summaryMatFile  = fullfile(dirOut, 'summaryTable.mat');
% summaryXlsxFile = fullfile(dirOut, 'summaryTable.xlsx');

%% Main loop
for iComb = 1:nComb

    mu_a = comb_mu_a(iComb);
    mu_b = comb_mu_b(iComb);
    mu_n = comb_mu_n(iComb);

    tag_mu_b = fmtExp(mu_b);
    tag_mu_n = fmtExp(mu_n);
    tag_mu_a = fmtExp(mu_a);

    fileName = sprintf('grid_muB_%s_muN_%s_muA_%s.mat', ...
        tag_mu_b, tag_mu_n, tag_mu_a);
    filePath = fullfile(dirCases, fileName);

    % Summary(iComb).idx        = iComb;
    % Summary(iComb).mu_a       = mu_a;
    % Summary(iComb).mu_b       = mu_b;
    % Summary(iComb).mu_n       = mu_n;
    % Summary(iComb).log10_mu_a = log10(mu_a);
    % Summary(iComb).log10_mu_b = log10(mu_b);
    % Summary(iComb).log10_mu_n = log10(mu_n);
    % Summary(iComb).fileName   = string(fileName);
    % Summary(iComb).filePath   = string(filePath);

    fprintf('\n====================================================\n');
    fprintf('Case %d / %d\n', iComb, nComb);
    fprintf('mu_b = 10^(%.2f), mu_n = 10^(%.2f), mu_a = 10^(%.2f)\n', ...
        log10(mu_b), log10(mu_n), log10(mu_a));
    fprintf('%s\n', fileName);

    if skipExisting && exist(filePath, 'file')
        fprintf('Already exists. Skipping.\n');
        % Summary(iComb).status = "skipped_existing";
        % update_summary_files(Summary, summaryMatFile, summaryXlsxFile); % better not
        continue;
    end

    tStart = tic;

    try
        result = run_one_case(mu_b, mu_n, mu_a, ...
            plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
            saveOutcomes, methodsRegu, saveMapsAndBSC);

        % Summary(iComb).elapsed_sec = toc(tStart);
        % Summary(iComb).status      = "ok";

        % Extract metrics from T_combined
        % Summary(iComb).mpe_a_3dof  = get_metric_from_table(result.T_combined, '3-DoF',   'a', 'mpe_homo');
        % Summary(iComb).mpe_b_3dof  = get_metric_from_table(result.T_combined, '3-DoF',   'b', 'mpe_homo');
        % Summary(iComb).mpe_n_3dof  = get_metric_from_table(result.T_combined, '3-DoF',   'n', 'mpe_homo');
        % 
        % Summary(iComb).mpe_a_2dofa = get_metric_from_table(result.T_combined, '2-DoF-a', 'a', 'mpe_homo');
        % Summary(iComb).mpe_b_2dofa = get_metric_from_table(result.T_combined, '2-DoF-a', 'b', 'mpe_homo');
        % Summary(iComb).mpe_n_2dofa = get_metric_from_table(result.T_combined, '2-DoF-a', 'n', 'mpe_homo');
        % 
        % Summary(iComb).mpe_a_2dofb = get_metric_from_table(result.T_combined, '2-DoF-b', 'a', 'mpe_homo');
        % Summary(iComb).mpe_b_2dofb = get_metric_from_table(result.T_combined, '2-DoF-b', 'b', 'mpe_homo');
        % Summary(iComb).mpe_n_2dofb = get_metric_from_table(result.T_combined, '2-DoF-b', 'n', 'mpe_homo');
        % 
        % Summary(iComb).mpe_a_2dofn = get_metric_from_table(result.T_combined, '2-DoF-n', 'a', 'mpe_homo');
        % Summary(iComb).mpe_b_2dofn = get_metric_from_table(result.T_combined, '2-DoF-n', 'b', 'mpe_homo');
        % Summary(iComb).mpe_n_2dofn = get_metric_from_table(result.T_combined, '2-DoF-n', 'n', 'mpe_homo');
        % 
        % Summary(iComb).gof_bsc_3dof  = get_bsc_metric(result.Tbsc, '3-DoF',   'diff_dB');
        % Summary(iComb).gof_bsc_2dofa = get_bsc_metric(result.Tbsc, '2-DoF-a', 'diff_dB');
        % Summary(iComb).gof_bsc_2dofb = get_bsc_metric(result.Tbsc, '2-DoF-b', 'diff_dB');
        % Summary(iComb).gof_bsc_2dofn = get_bsc_metric(result.Tbsc, '2-DoF-n', 'diff_dB');

        % Save one case per file
        save(filePath, '-struct', 'result', '-v7.3');
        fprintf('Saved case: %s\n', fileName);

    catch ME
        % Summary(iComb).elapsed_sec = toc(tStart);
        % Summary(iComb).status      = "error";
        fprintf(2, 'Error in case %d: %s\n', iComb, ME.message);

        errFile = fullfile(dirCases, sprintf('ERROR_%s.mat', fileName(1:end-4)));
        save(errFile, 'mu_a', 'mu_b', 'mu_n', 'ME');
    end

    % update_summary_files(Summary, summaryMatFile, summaryXlsxFile); % better not
end

fprintf('\nAll done.\n');
fprintf('Cases folder: %s\n', dirCases);
% fprintf('Summary MAT:  %s\n', summaryMatFile);
% fprintf('Summary XLSX: %s\n', summaryXlsxFile);

end

%% ========================================================================
function result = run_one_case(mu_b, mu_n, mu_a, ...
    plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
    saveOutcomes, methodsRegu, saveMapsAndBSC)

close all;
warning('off');

%% MU SETUP [mu_b, mu_n, mu_a]
mu_setup = [mu_b, mu_n, mu_a; ... % 2-DoF-a
            mu_b, mu_n, mu_a; ... % 2-DoF-n
            mu_b, mu_n, mu_a; ... % 2-DoF-b
            mu_b, mu_n, mu_a];    % 3-DoF

%% INITIALIZATION
Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB; %#ok<NASGU>
range_bmode     = [-80 0];

% Directory outcomes (not used in grid search, kept for compatibility)
if saveOutcomes
    dirOutcomes = './out/dofJournal/p1p3/';
    if methodsRegu
        dirOutcomes = fullfile(dirOutcomes, 'regu');
    else
        dirOutcomes = fullfile(dirOutcomes, 'no_regu');
    end
    if ~exist(dirOutcomes, 'dir')
        mkdir(dirOutcomes);
    end
end

% Data directories
if ispc
    % Windows (local)
    dirData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\data4Prociencia\phantoms';
elseif isunix
    % Linux cluster
    dirData = '/mnt/nfs2/emiranda/data4Prociencia/phantoms';  % <-- TO CHANGE ^^
else
    error('Unknown operating system');
end

folderDataSam   = 'UIUC';
folderDataRef   = 'UIUC';

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
label_methods = {'3-DoF', '2-DoF_{b,n}', '2-DoF_{n,a}', '2-DoF_{b,a}'};

bsc_results     = cell(4, length(methods));
maps_results    = cell(4, length(methods));

bsc_results(1, :)   = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'), sprintf('2-DoF "n"')};
maps_results(1, :)  = {sprintf('3-DoF'), sprintf('2-DoF "a"'), sprintf('2-DoF "b"'), sprintf('2-DoF "n"')};

%% LOAD DATA

j_sam = 1.0;
j_ref = 1.0;

numPhantomSam = '2'; 
numPhantomRef = '3';
idName = ['p', numPhantomSam, 'p', numPhantomRef];

% SAMPLE SPECS
switch str2double(numPhantomSam)
    case 1
        alpha_sam   = 0.4; % [dB/cm/MHz]
        % BW 3.5-7.5
        % delta_b_prior       = log(db2pow(20)); % from RPM 18.842471
        % delta_n_prior       = -2.28; % from RPM method -2.285622   

        delta_b_prior       = log(db2pow(18.8)); % from RPM new BW ^^
        delta_n_prior       = -2.28; % from RPM method -2.285622     
    case 2
        alpha_sam   = 0.1; % [dB/cm/MHz]
        % BW 3.5-7.5
        delta_b_prior       = log(db2pow(-4.5)); % from RPM -3.602607
        delta_n_prior       = 0.38; % from RPM method 0.267361   
    case 3
        alpha_sam   = 0.1; % [dB/cm/MHz]
end
% REFERENCE SPECS
switch str2double(numPhantomRef)
    case 1
        alpha_ref   = 0.4; % [dB/cm/MHz]
    case 2
        alpha_ref   = 0.1; % [dB/cm/MHz]
    case 3
        alpha_ref   = 0.1; % [dB/cm/MHz]
end

delta_alpha_prior   = alpha_sam - alpha_ref; % [dB/cm/MHz]

% STABLE
samName = ['rf_phantom', numPhantomSam];
SAM     = load (fullfile(dirData, folderDataSam, samName));
SAM.rf  = SAM.rf(:,:,end-9:end); % Last frames 
SAM.acs = alpha_sam;

refName = ['rf_phantom', numPhantomRef];
REF     = load (fullfile(dirData, folderDataRef, refName));
REF.rf  = REF.rf(:,:,end-9:end); % Last frames 
REF.acs = alpha_ref;

% B-MODE CHECK
bmode_sam = db(hilbert(SAM.rf(:,:,1)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(hilbert(REF.rf(:,:,1)));
bmode_ref = bmode_ref - max(bmode_ref(:));


%% SPECTRAL METHOD PARAMETERS
clear pars

pars.P           = 512; % NFFT for 10wl is 256, 20wl 512
pars.bw          = [3.9 7.6]; % [MHz]
pars.overlap     = 0.8;
pars.blocksize   = 14; % wavelengths
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

pars.x_roi       = [-15 15]*1E-3; % [m] 
pars.z_roi       = [5 35]*1E-3; % [m] % Same than sonix


if plotBmode
    figure;

    subplot(121)
    imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image"), hold on;
    rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
    clim(range_bmode)
    cb = colorbar;
    cb.Label.String = 'dB';
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('SAM')
    colormap('gray')

    subplot(122)
    imagesc(REF.x*1E3, REF.z*1E3, bmode_ref, range_bmode), axis("image");
    rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
    clim(range_bmode)
    cb = colorbar;
    cb.Label.String = 'dB';
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('REF')
    colormap('gray')
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


%% bsc
% pathBSC_RPM = '/mnt/nfs2/emiranda/phantoms/UIUC/RPM';
% pathBSC_RPM = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\data4Prociencia\phantoms\UIUC\RPM';

% Data directories
if ispc
    % Windows (local)
    pathBSC_RPM = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\data4Prociencia\phantoms\UIUC\RPM';
elseif isunix
    % Linux cluster
    pathBSC_RPM = '/mnt/nfs2/emiranda/data4Prociencia/phantoms/UIUC/RPM';  % <-- TO CHANGE ^^
else
    error('Unknown operating system');
end

nameCaseRPM = sprintf('p%sp%s_RPM', numPhantomSam, numPhantomRef);
% BSC     = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD** *better preload)
load(fullfile(pathBSC_RPM, nameCaseRPM));


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




%% GENERAL REGULARIZATION SETTINGS
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0;
par_rpl.ini_tol    = 1e-5;
par_rpl.df_op      = 0;
par_rpl.ini_method = 1;

%% FOR METHOD LOOP
for iMet = 1:length(methods)

    estim_method = methods{iMet};

    %% COMPENSATE 2-DoF-a
    if strcmp(estim_method, '2-DoF-a')

        if methodsRegu
            mu_rpl_tv = mu_setup(1,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR); %#ok<ASGLU>

        comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron( speye(p*q), ones(size(f)) );
        Z = kron( speye(p*q), log(f) );

        u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);

        if par_rpl.df_op == 1
            dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
            dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
            Dy = sparse(kron(speye(q),dy));
        else
            dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
            Dy = sparse(kron(speye(q),dy));
        end

        g = u_opt(1:p*q);
        s = u_opt(p*q+1:2*p*q);
        a_Np2dB = delta_alpha_prior*ones(p*q, 1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

    %% COMPENSATE 2-DoF-n
    elseif strcmp(estim_method, '2-DoF-n')

        if methodsRegu
            mu_rpl_tv = mu_setup(2,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR); %#ok<ASGLU>

        comp_ref    = comp_ref_n_bsc(delta_n_prior, band, p, q);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron( speye(p*q), ones(size(f)) );
        W = kron( speye(p*q), -4*f );

        u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));

        g = u_opt(1:p*q);
        a = u_opt(p*q+1:2*p*q);

        s = delta_n_prior*ones(p*q, 1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% COMPENSATE 2-DoF-b
    elseif strcmp(estim_method, '2-DoF-b')

        if methodsRegu
            mu_rpl_tv = mu_setup(3,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR); %#ok<ASGLU>

        comp_ref    = comp_ref_b_bsc(delta_b_prior);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        Z = kron( speye(p*q), log(f) );
        W = kron( speye(p*q), -4*f );

        u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);

        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));

        s = u_opt(1:p*q);
        a = u_opt(p*q+1:2*p*q);

        g = delta_b_prior*ones(p*q, 1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% COMPENSATE 3-DoF
    elseif strcmp(estim_method, '3-DoF')

        if methodsRegu
            mu_rpl_tv = mu_setup(4,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR); %#ok<ASGLU>

        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR.*comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron( speye(p*q), ones(size(f)) );
        Z = kron( speye(p*q), log(f) );
        W = kron( speye(p*q), -4*f.^j_sam );

        u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

        if par_rpl.df_op == 1
            dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
            dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
            Dy = sparse(kron(speye(q),dy));
        else
            dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
            Dy = sparse(kron(speye(q),dy));
        end

        g = u_opt(1:p*q);
        s = u_opt(p*q+1:2*p*q);
        a = u_opt(2*p*q+1:3*p*q);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);
    end

    %% QUS PARAMETERS
    b_ratio     = reshape(exp(g), p, q);
    b_ratio_dB  = 10*log10(b_ratio);
    alpha_ratio = reshape(a_Np2dB, p, q);
    s_ratio     = reshape(s, p, q);

    acs_sam   = alpha_ratio + alpha_ref;

    calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

    [m_b, s_b, cv_b] = deal(calc2dStats{1}(b_ratio), calc2dStats{2}(b_ratio), calc2dStats{3}(b_ratio));
    if plotBSCdB
        [m_b, s_b, cv_b] = deal(calc2dStats{1}(b_ratio_dB), calc2dStats{2}(b_ratio_dB), calc2dStats{3}(b_ratio_dB));
    end

    [m_n, s_n, cv_n] = deal(calc2dStats{1}(s_ratio),   calc2dStats{2}(s_ratio),   calc2dStats{3}(s_ratio));
    [m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam),   calc2dStats{2}(acs_sam),   calc2dStats{3}(acs_sam));

    fprintf('-----%s---\n', estim_method);
    fprintf('alpha_s    : %.3f +- %.4f, %%CV = %.4f\n', round(m_a, 3), round(s_a, 4), round(cv_a, 4));
    if plotBSCdB
        fprintf('Delta b[dB]: %.3f +- %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
    else
        fprintf('Delta b    : %.3f +- %.4f, %%CV = %.4f\n', round(m_b, 3), round(s_b, 4), round(cv_b, 4));
    end
    fprintf('Delta n    : %.4f +- %.4f, %%CV = %.4f\n', round(m_n, 4), round(s_n, 4), round(cv_n, 4));
    fprintf('--------\n');

    %% OPTIONAL MAP PLOTS
    Xaxis   = spectralData_ref.lateral;
    Zaxis   = spectralData_ref.depth;
    cm      = 1e2;
    fontSize = 16;

    if plotMaps
        figure;
        set(gcf,'units','normalized','outerposition',[0 0.15 1 0.75]); box on;
        sgtitle(label_methods{iMet}, 'FontSize', fontSize+2, 'FontWeight', 'bold');

        subplot(1,3,1)
        imagesc(Xaxis*cm, Zaxis*cm, acs_sam), colorbar
        axis("image");
        xlabel('Lateral [cm]'), ylabel('Depth [cm]'), colormap turbo;
        title({'$\alpha_s$:', ...
            [num2str(round(m_a, 3)), ' $\pm$ ', num2str(round(s_a, 3)), ', CV = ', num2str(round(cv_a, 3))]}, ...
            'Interpreter', 'latex');

        h2 = colorbar;
        ylabel(h2,'dB\cdotcm^{-1}\cdotMHz^{-1}','FontSize', fontSize);
        set(gca,'fontsize',fontSize)

        subplot(1,3,2)
        imagesc(Xaxis*cm, Zaxis*cm, b_ratio)
        h2 = colorbar;
        if plotBSCdB
            imagesc(Xaxis*cm, Zaxis*cm, b_ratio_dB)
            h2 = colorbar;
            ylabel(h2,'dB','FontSize', fontSize);
        end
        axis("image");
        xlabel('Lateral [cm]'), colormap turbo;
        title({'$\Delta b$:', ...
            [num2str(round(m_b, 3)), ' $\pm$ ', num2str(round(s_b, 3)), ', CV = ', num2str(round(cv_b, 3))]}, ...
            'Interpreter', 'latex');
        set(gca,'fontsize',fontSize)

        subplot(1,3,3)
        imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
        axis("image");
        xlabel('Lateral [cm]'), colormap turbo;
        title({'$\Delta n$:', ...
            [num2str(round(m_n, 3)), ' $\pm$ ', num2str(round(s_n, 3)), ', CV = ', num2str(round(cv_n, 3))]}, ...
            'Interpreter', 'latex');
        set(gca,'fontsize',fontSize)
    end

    %% OPTIONAL OVERLAY PLOTS
    if plotBmodeOverl
        fontSize = 16;

        figure;
        set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(label_methods{iMet}, 'FontSize', fontSize+2, 'FontWeight', 'bold');

        % alpha
        units           = 1E2;
        bmodeFull       = bmode_sam;
        colorImg        = acs_sam;
        range_img       = [0.1 1.2];
        transparency    = 0.65;
        x_img           = spectralData_sam.lateral*units;
        z_img           = spectralData_sam.depth*units;
        xFull           = SAM.x*units;
        zFull           = SAM.z*units;
        [X, Z] = meshgrid(xFull, zFull);
        roi = and(X >= x_img(1), X <= x_img(end)) & and(Z >= z_img(1), Z <= z_img(end));

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
        title('$\alpha_s$', 'Interpreter', 'latex')
        set(gca,'fontsize',fontSize)

        % b
        colorImg = b_ratio;
        range_img = [];
        if plotBSCdB
            colorImg = b_ratio_dB;
        end

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        if plotBSCdB
            hColor.Label.String = 'dB';
        else
            hColor.Label.String = '';
        end
        title('$\Delta b$', 'Interpreter','latex')
        set(gca,'fontsize',fontSize)

        % n
        colorImg = s_ratio;

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, range_img, ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        hColor.Label.String = 'a.u.';
        title('$\Delta s$', 'Interpreter','latex')
        set(gca,'fontsize',fontSize)
    end

    %% BSC RECONSTRUCTION
    g_est = median(b_ratio(:));
    s_est = median(s_ratio(:));

    freq = spectralData_sam.band;
    bsc_est_powlaw = g_est.*(freq.^s_est);
    bsc_results{2, iMet} = bsc_est_powlaw;

    maps_results{2, iMet} = acs_sam;
    maps_results{3, iMet} = b_ratio_dB;
    maps_results{4, iMet} = s_ratio;
end

%% METRICS MAPS TABLE
delta_b_theo = exp(delta_b_prior); % usually before
delta_n_theo = delta_n_prior;
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

clear m_3dof m_2dofa m_2dofb m_2dofn

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



%% OPTIONAL SAVE OUTCOMES
if saveOutcomes
    if methodsRegu
        nameExcel = ['metricsRegu_',idName,'UIUC.xlsx']; 
    else
        nameExcel = ['metricsNoRegu_',idName,'UIUC.xlsx']; 
    end
    excelFile = fullfile(dirOutcomes, nameExcel);
    writetable(T_combined, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

    if methodsRegu
        nameExcel = ['metricsBSCRegu_',idName,'UIUC.xlsx']; 
    else
        nameExcel = ['metricsBSCNoRegu_',idName,'UIUC.xlsx']; 
    end
    excelFile = fullfile(dirOutcomes, nameExcel);
    writetable(Tbsc, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);
end

%% Prepare result struct to save one case
result = struct();

result.mu_a = mu_a;
result.mu_b = mu_b;
result.mu_n = mu_n;
result.mu_setup = mu_setup;

result.methods = methods;
result.label_methods = label_methods;

result.alpha_sam = alpha_sam;
result.alpha_ref = alpha_ref;
% result.b_sam = b_sam;
% result.b_ref = b_ref;
% result.n_sam = n_sam;
% result.n_ref = n_ref;

result.delta_alpha_prior = delta_alpha_prior;
result.delta_b_prior     = delta_b_prior;
result.delta_n_prior     = delta_n_prior;

result.pars = pars;
result.par_rpl = par_rpl;

result.T_combined = T_combined;
result.Tbsc       = Tbsc;

if saveMapsAndBSC
    result.maps_results = maps_results;
    result.bsc_results  = bsc_results;
end

% Compact summary block
result.summaryInfo = struct();
result.summaryInfo.mu_a       = mu_a;
result.summaryInfo.mu_b       = mu_b;
result.summaryInfo.mu_n       = mu_n;
result.summaryInfo.log10_mu_a = log10(mu_a);
result.summaryInfo.log10_mu_b = log10(mu_b);
result.summaryInfo.log10_mu_n = log10(mu_n);

end

%% ========================================================================
function val = get_metric_from_table(T, methodName, paramName, metricName)
%GET_METRIC_FROM_TABLE Extract one scalar metric from long-format T_combined
%
% Inputs:
%   T          : table like struct2table(MetricsParam)
%   methodName : e.g. '3-DoF', '2-DoF-a'
%   paramName  : e.g. 'a', 'b', 'n'
%   metricName : e.g. 'MPE', 'RMSE', 'MAE', 'NRMSE'
%
% Output:
%   val        : scalar metric value, NaN if not found

    val = NaN;

    try
        % Check required columns
        if ~ismember('method', T.Properties.VariableNames)
            warning('Column "method" not found.');
            return;
        end

        if ~ismember('param', T.Properties.VariableNames)
            warning('Column "param" not found.');
            return;
        end

        if ~ismember(metricName, T.Properties.VariableNames)
            warning('Metric column "%s" not found.', metricName);
            return;
        end

        % Convert to string to avoid categorical/char issues
        methodCol = string(T.method);
        paramCol  = string(T.param);

        idx = (methodCol == string(methodName)) & (paramCol == string(paramName));

        if ~any(idx)
            warning('No row found for method="%s", param="%s".', methodName, paramName);
            return;
        end

        tmp = T{find(idx,1,'first'), metricName};

        if iscell(tmp)
            tmp = tmp{1};
        end

        if isnumeric(tmp) && isscalar(tmp)
            val = tmp;
        else
            warning('Metric "%s" for method="%s", param="%s" is not scalar numeric.', ...
                metricName, methodName, paramName);
        end

    catch ME
        warning('Error extracting metric: %s', ME.message);
        val = NaN;
    end
end

%% ========================================================================
function val = get_bsc_metric(Tbsc, methodName, metricName)
val = NaN;
try
    idx = strcmp(string(Tbsc.method), methodName);
    if any(idx) && ismember(metricName, Tbsc.Properties.VariableNames)
        tmp = Tbsc{find(idx,1), metricName};
        if iscell(tmp)
            val = tmp{1};
        else
            val = tmp;
        end
        if ~isscalar(val)
            val = NaN;
        end
    end
catch
    val = NaN;
end
end

%% ========================================================================
function update_summary_files(Summary, summaryMatFile, summaryXlsxFile)
SummaryTable = struct2table(Summary);
save(summaryMatFile, 'Summary', 'SummaryTable', '-v7.3');
try
    writetable(SummaryTable, summaryXlsxFile);
catch
    warning('Could not write summary xlsx in this iteration.');
end
end