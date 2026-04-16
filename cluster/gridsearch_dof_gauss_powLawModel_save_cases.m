function gridsearch_dof_gauss_powLawModel_save_cases()
% GRIDSEARCH_DOF_GAUSS_POWLAWMODEL_SAVE_CASES
% Grid search for:
%   - Gaussian simulation
%   - Power law model
% Saves each case separately as .mat
%
% Based on EAMZ code structure.

clear; clc; close all;
warning('off');

cd('../'); % when you run it usually goes inside ./parentFolder/script.m
disp(pwd)
addpath(genpath(pwd))
addpath(genpath('/opt/MATLAB Add-Ons'))

%% ========================= USER CONFIG =========================
list_mu_a = 10.^(-3:0.5:6);
list_mu_b = 10.^(-3:0.5:6);
list_mu_n = 10.^(-3:0.5:6);

% Small debug example:
% list_mu_a = 10.^4;
% list_mu_b = [1000];
% list_mu_n = 10.^3;

% Output directory
if ispc
    % Windows (local)
    dirOut = './out/dofJournal/simuGauss/powLawModel';
elseif isunix
    % Linux cluster
    dirOut = '/mnt/nfs2/emiranda/proj/dof26/out/simuGauss/powLawModel';  % explicit
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

skipExisting = true;

plotBmode       = false;
plotBmodeOverl  = false;
plotBSCdB       = true;
plotMaps        = false;
plotBoxplots    = false;
plotDeltaBSC    = false;

saveOutcomes    = false;
methodsRegu     = true;
bsc_gt_by_rpm   = true;

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

%% Main loop
for iComb = 5898:nComb
% for iComb = 1:nComb

    mu_a = comb_mu_a(iComb);
    mu_b = comb_mu_b(iComb);
    mu_n = comb_mu_n(iComb);

    tag_mu_b = fmtExp(mu_b);
    tag_mu_n = fmtExp(mu_n);
    tag_mu_a = fmtExp(mu_a);

    fileName = sprintf('grid_muB_%s_muN_%s_muA_%s.mat', ...
        tag_mu_b, tag_mu_n, tag_mu_a);
    filePath = fullfile(dirCases, fileName);

    fprintf('\n====================================================\n');
    fprintf('Case %d / %d\n', iComb, nComb);
    fprintf('mu_b = 10^(%.2f), mu_n = 10^(%.2f), mu_a = 10^(%.2f)\n', ...
        log10(mu_b), log10(mu_n), log10(mu_a));
    fprintf('%s\n', fileName);

    if skipExisting && exist(filePath, 'file')
        fprintf('Already exists. Skipping.\n');
        continue;
    end

    try
        result = run_one_case_gauss_powLaw(mu_b, mu_n, mu_a, ...
            plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
            plotBoxplots, plotDeltaBSC, ...
            saveOutcomes, methodsRegu, bsc_gt_by_rpm, saveMapsAndBSC);

        save(filePath, '-struct', 'result', '-v7.3');
        fprintf('Saved case: %s\n', fileName);

    catch ME
        fprintf(2, 'Error in case %d: %s\n', iComb, ME.message);
        errFile = fullfile(dirCases, sprintf('ERROR_%s.mat', fileName(1:end-4)));
        save(errFile, 'mu_a', 'mu_b', 'mu_n', 'ME');
    end
end

fprintf('\nAll done.\n');
fprintf('Cases folder: %s\n', dirCases);

end

%% ========================================================================
function result = run_one_case_gauss_powLaw(mu_b, mu_n, mu_a, ...
    plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
    plotBoxplots, plotDeltaBSC, ...
    saveOutcomes, methodsRegu, bsc_gt_by_rpm, saveMapsAndBSC)

close all;
warning('off');

%% MU SETUP [mu_b, mu_n, mu_a]
mu_setup = [mu_b, mu_n, mu_a; ... % 2-DoF-a
            mu_b, mu_n, mu_a; ... % 2-DoF-b
            mu_b, mu_n, mu_a; ... % 2-DoF-n
            mu_b, mu_n, mu_a];    % 3-DoF

%% INITIALIZATION

Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB; %#ok<NASGU>
range_bmode     = [-80 0];

if saveOutcomes
    dirOutcomes = './out/dofJournal/simuGauss/powLawModel';
    if ~exist(dirOutcomes, 'dir')
        mkdir(dirOutcomes);
    end
end

%% DATA DIRECTORIES
% Data directories
if ispc
    % Windows (local)
    dirData         = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\dataLIM\data4Prociencia\simus';
elseif isunix
    % Linux cluster
    dirData = '/mnt/nfs2/emiranda/proj/dof26/data/';  % <-- TO CHANGE ^^
else
    error('Unknown operating system');
end

folderDataSam   = 'gaussSimu';
folderDataRef   = 'gaussSimu';

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
label_methods = {'3-DoF', '2-DoF_{b,n}', '2-DoF_{n,a}', '2-DoF_{b,a}'};

bsc_results     = cell(4, length(methods));
maps_results    = cell(4, length(methods));

bsc_results(1, :) = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
maps_results(1, :) = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

%% LOAD DATA
rf_sam_name = 'rf_sd2pcSCALE4_bsc2_gauss_mask7_att0p5';
alpha_sam   = 0.5;
j_sam       = 1.1;

SAM             = load(fullfile(dirData, folderDataSam, rf_sam_name));
SAM.acs         = alpha_sam;
SAM.alpha_power = j_sam;

rf_ref_name = 'rf_sd2pcSCALE4_bsc4_att0p4';
alpha_ref   = 0.4;
j_ref       = 1.1;

REF             = load(fullfile(dirData, folderDataRef, rf_ref_name));
REF.acs         = alpha_ref;
REF.alpha_power = j_ref;

%% PRIORS / GT
delta_alpha_prior = alpha_sam - alpha_ref;
delta_b_prior     = log(db2pow(6.416335));
delta_n_prior     = -1.280963;

delta_b_theo = 4.605340;
delta_n_theo = delta_n_prior;

%% BMODE
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL PARAMETERS
clear pars
pars.P           = 4096;
pars.overlap     = 0.8;
pars.z_roi       = [5 42.5]*1E-3;
pars.x_roi       = [-18 18]*1E-3;
pars.window_type = 3;
pars.saran_layer = false;

pars.bw          = [3 9.2];
pars.blocksize   = 14;

if plotBmode
    figure;

    subplot(121)
    imagesc(SAM.x*1E3, SAM.z*1E3, bmode_sam, range_bmode), axis("image"), hold on;
    rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
    clim(range_bmode)
    cb = colorbar; cb.Label.String = 'dB';
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('SAM')
    colormap('gray')

    subplot(122)
    imagesc(REF.x*1E3, REF.z*1E3, bmode_ref, range_bmode), axis("image"), hold on;
    rectangle('Position', 1E3*[pars.x_roi(1) pars.z_roi(1) pars.x_roi(2)-pars.x_roi(1) pars.z_roi(2)-pars.z_roi(1)], ...
        'EdgeColor','w', 'LineWidth', 2, 'LineStyle','--'), hold off;
    clim(range_bmode)
    cb = colorbar; cb.Label.String = 'dB';
    xlabel('Lateral [mm]'), ylabel('Depth [mm]');
    title('REF')
    colormap('gray')
end

%% POWER SPECTRA
spectralData_sam = calc_powerSpectra_vSimple(SAM, pars);
S_sam = spectralData_sam.powerSpectra;

spectralData_ref = calc_powerSpectra_vSimple(REF, pars);
S_ref = spectralData_ref.powerSpectra;

SR_emz = S_sam ./ S_ref;
SR = permute(SR_emz,[3,1,2]); clear SR_emz

%% BSC GT BY RPM
if bsc_gt_by_rpm
    BSC     = calculateBSC_RPM_fast(SAM, REF, pars);
    bsc_rpm = BSC.BSCcurve_Uni(:,2);
    freq_rpm = BSC.band;

    coeffs   = polyfit(log(freq_rpm), log(bsc_rpm), 1);
    d_n      = coeffs(1);
    ln_db    = coeffs(2);
    d_b      = exp(ln_db);

    fprintf('-----RPM PowLaw (b.(f.^n))-----\n')
    fprintf('Δb           = %f\n', d_b);
    fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b));
    fprintf('Δn           = %f\n', d_n);
    fprintf('---------\n')

    bsc_rpm_powlaw = d_b*(freq_rpm.^d_n); %#ok<NASGU>

    coeffs = polyfit(-freq_rpm.^2, log(bsc_rpm), 1);
    d_s    = coeffs(1);
    ln_dg  = coeffs(2);
    d_g    = exp(ln_dg);

    fprintf('-----RPM Gauss (g.exp(-s.f^2))-----\n')
    fprintf('Δg           = %f\n', d_g);
    fprintf('Δg      [dB] = %f\n', 10*log10(d_g));
    fprintf('Δs           = %f\n', d_s);
    fprintf('---------\n')

    bsc_rpm_gauss = d_g*exp(-d_s*freq_rpm.^2); %#ok<NASGU>
else
    error('This script currently expects bsc_gt_by_rpm = true');
end

%% REGULARIZATION SETTINGS
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0;
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1;

%% METHOD LOOP
for iMet = 1:length(methods)

    estim_method = methods{iMet};

    %% 2-DoF-a
    if strcmp(estim_method, '2-DoF-a')

        if methodsRegu
            mu_rpl_tv = [mu_b; mu_n; mu_a];
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron(speye(p*q), ones(size(f)));
        Z = kron(speye(p*q), log(f));

        u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);

        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));

        g = u_opt(1:p*q);
        s = u_opt(p*q+1:2*p*q);
        a_Np2dB = delta_alpha_prior*ones(p*q,1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

    %% 2-DoF-n
    elseif strcmp(estim_method, '2-DoF-n')

        if methodsRegu
            mu_rpl_tv = [mu_b; mu_n; mu_a];
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_ref    = comp_ref_n_bsc(delta_n_prior, band, p, q);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron(speye(p*q), ones(size(f)));
        W = kron(speye(p*q), -4*f.^j_sam);

        u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));

        g = u_opt(1:p*q);
        a = u_opt(p*q+1:2*p*q);
        s = delta_n_prior*ones(p*q,1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% 2-DoF-b
    elseif strcmp(estim_method, '2-DoF-b')

        if methodsRegu
            mu_rpl_tv = [mu_b; mu_n; mu_a];
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_ref    = comp_ref_g_bsc(delta_b_prior);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        Z = kron(speye(p*q), log(f));
        W = kron(speye(p*q), -4*f.^j_sam);

        u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);

        if par_rpl.df_op == 1
            dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
            dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
            Dy = sparse(kron(speye(q),dy));
        else
            dy = diag(ones(p-1,1),1) - diag([ones(p-1,1);0]);
            Dy = sparse(kron(speye(q),dy));
        end

        s = u_opt(1:p*q);
        a = u_opt(p*q+1:2*p*q);
        g = delta_b_prior*ones(p*q,1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% 3-DoF
    elseif strcmp(estim_method, '3-DoF')

        if methodsRegu
            mu_rpl_tv = [mu_b; mu_n; mu_a];
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        X = kron(speye(p*q), ones(size(f)));
        Z = kron(speye(p*q), log(f));
        W = kron(speye(p*q), -4*f.^j_sam);

        u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
        [u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));

        g = u_opt(1:p*q);
        s = u_opt(p*q+1:2*p*q);
        a = u_opt(2*p*q+1:3*p*q);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);
    end

    %% QUS PARAMETERS
    g_ratio     = reshape(exp(g), p, q);
    g_ratio_dB  = 10*log10(g_ratio);
    alpha_ratio = reshape(a_Np2dB, p, q);
    s_ratio     = reshape(s, p, q);

    acs_sam = alpha_ratio + alpha_ref;

    calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100*std(x(:))/mean(x(:))};

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

    %% OPTIONAL MAPS
    Xaxis = spectralData_ref.lateral;
    Zaxis = spectralData_ref.depth;
    cm = 1e2;
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
        imagesc(Xaxis*cm, Zaxis*cm, g_ratio)
        h2 = colorbar;
        if plotBSCdB
            imagesc(Xaxis*cm, Zaxis*cm, g_ratio_dB)
            h2 = colorbar;
            ylabel(h2,'dB','FontSize', fontSize);
        end
        axis("image");
        xlabel('Lateral [cm]'), colormap turbo;
        title({'$\Delta b$:', ...
            [num2str(round(m_g, 3)), ' $\pm$ ', num2str(round(s_g, 3)), ', CV = ', num2str(round(cv_g, 3))]}, ...
            'Interpreter', 'latex');
        set(gca,'fontsize',fontSize)

        subplot(1,3,3)
        imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
        axis("image");
        xlabel('Lateral [cm]'), colormap turbo;
        title({'$\Delta n$:', ...
            [num2str(round(m_s, 3)), ' $\pm$ ', num2str(round(s_s, 3)), ', CV = ', num2str(round(cv_s, 3))]}, ...
            'Interpreter', 'latex');
        set(gca,'fontsize',fontSize)
    end

    %% OPTIONAL OVERLAY
    if plotBmodeOverl
        fontSize = 16;

        figure;
        set(gcf,'units','normalized','outerposition',[0 0.1 1 0.8]); box on;
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        sgtitle(label_methods{iMet}, 'FontSize', fontSize+2, 'FontWeight', 'bold');

        units           = 1E2;
        bmodeFull       = bmode_sam;
        x_img           = spectralData_sam.lateral*units;
        z_img           = spectralData_sam.depth*units;
        xFull           = SAM.x*units;
        zFull           = SAM.z*units;
        [X, Z] = meshgrid(xFull, zFull);
        roi = and(X >= x_img(1), X <= x_img(end)) & and(Z >= z_img(1), Z <= z_img(end));
        transparency = 0.65;

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, acs_sam, range_bmode, [0.1 1.2], ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        hColor.Label.String = 'dB\cdotcm^{-1}\cdotMHz^{-1}';
        title('$\alpha_s$', 'Interpreter', 'latex')
        set(gca,'fontsize',fontSize)

        colorImg = g_ratio;
        if plotBSCdB
            colorImg = g_ratio_dB;
        end
        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, colorImg, range_bmode, [], ...
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

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, s_ratio, range_bmode, [], ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        hColor.Label.String = 'a.u.';
        title('$\Delta n$', 'Interpreter','latex')
        set(gca,'fontsize',fontSize)
    end

    %% BSC RECONSTRUCTION USING POWER LAW
    freq = spectralData_sam.band;
    g_est = median(g_ratio(:));
    s_est = median(s_ratio(:));

    bsc_est_powlaw = g_est .* (freq.^s_est);

    bsc_results{2, iMet} = bsc_est_powlaw;
    maps_results{2, iMet} = acs_sam;
    maps_results{3, iMet} = g_ratio_dB;
    maps_results{4, iMet} = s_ratio;
end

%% OPTIONAL BOXPLOTS
if plotBoxplots
    font_size  = 34;
    numMethods = size(maps_results, 2);

    acs_data     = cell(1, numMethods);
    b_ratio_data = cell(1, numMethods);
    n_ratio_data = cell(1, numMethods);

    for iMet = 1:numMethods
        acs_data{iMet}     = maps_results{2, iMet}(:);
        b_ratio_data{iMet} = maps_results{3, iMet}(:);
        n_ratio_data{iMet} = maps_results{4, iMet}(:);
    end

    acs_mat     = padconcatenation(acs_data, NaN, 1);
    b_ratio_mat = padconcatenation(b_ratio_data, NaN, 1);
    n_ratio_mat = padconcatenation(n_ratio_data, NaN, 1);

    method_labels = { ...
        '\mathrm{3\textrm{-}DoF}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
    };

    % alpha
    acs_mat_filtered = acs_mat(:, [1, 3, 4]);
    method_labels_a  = method_labels([1, 3, 4]);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);
    box on;
    boxplot(acs_mat_filtered, 'Labels', method_labels_a);
    yline(alpha_sam, 'k--')
    ax = gca;
    ax.XTickLabel = {''};
    ax.XTickLabelMode = 'manual';
    if methodsRegu, ylim([0.05 1.01]); else, ylim([-20 20]); end
    xt = get(ax, 'XTick');
    for i = 1:length(method_labels_a)
        text(xt(i), ax.YLim(1)*0.005, ['$' method_labels_a{i} '$'], ...
            'Interpreter','latex', 'HorizontalAlignment','center', ...
            'FontSize', font_size, 'FontWeight','bold')
    end
    title('\bf\alpha')
    ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
    set(gca, 'FontSize', font_size);

    % b
    b_ratio_mat_filtered = b_ratio_mat(:, [1, 2, 4]);
    method_labels_b      = method_labels([1, 2, 4]);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);
    box on;
    boxplot(b_ratio_mat_filtered, 'Labels', method_labels_b);
    yline(10*log10(delta_b_theo), 'k--')
    ax = gca;
    ax.XTickLabel = {''};
    ax.XTickLabelMode = 'manual';
    xt = get(ax, 'XTick');
    if methodsRegu, ylim([-4 7.5]); else, ylim([-60 20]); end
    for i = 1:length(method_labels_b)
        text(xt(i), 1.15*ax.YLim(1), ['$' method_labels_b{i} '$'], ...
            'Interpreter','latex', 'HorizontalAlignment','center', ...
            'FontSize', font_size, 'FontWeight','bold')
    end
    title('\bf\Deltab');
    ylabel('\Deltab [dB]');
    set(gca, 'FontSize', font_size);

    % n
    n_ratio_mat_filtered = n_ratio_mat(:, [1, 2, 3]);
    method_labels_n      = method_labels([1, 2, 3]);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);
    box on;
    boxplot(n_ratio_mat_filtered, 'Labels', method_labels_n);
    yline(delta_n_theo, 'k--')
    ax = gca;
    ax.XTickLabel = {''};
    ax.XTickLabelMode = 'manual';
    xt = get(ax, 'XTick');
    if ~methodsRegu, ylim([-15 15]); end
    for i = 1:length(method_labels_n)
        text(xt(i), 1.15*ax.YLim(1), ['$' method_labels_n{i} '$'], ...
            'Interpreter','latex', 'HorizontalAlignment','center', ...
            'FontSize', font_size, 'FontWeight','bold')
    end
    title('\bf\Deltan');
    ylabel('\Deltan [a.u.]');
    set(gca, 'FontSize', font_size);
end

%% METRICS BSC
BSC_gt = bsc_rpm;
freq = spectralData_sam.band;
diff_fit_dB = @(bsc_pred, bsc_gt) mean(abs(10*log10(bsc_pred) - 10*log10(bsc_gt)));

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
Tbsc            = struct2table(MetricsBSC);
Tbsc.method     = categorical(Tbsc.method);
Tbsc.param      = categorical(Tbsc.param);

%% OPTIONAL DELTA BSC PLOT
if plotDeltaBSC
    xlim_range = pars.bw + 0.05*[-1 1];
    ylim_range = [10^-0.65 10^1.1];
    line_width = 3.5;
    font_size  = 32;

    bsc_results{1,1} = sprintf('3-DoF     (GoF_{dB} = %.2f)', MetricsBSC(1).diff_dB);
    bsc_results{1,2} = sprintf('2-DoF_{b,n} (GoF_{dB} = %.2f)', MetricsBSC(2).diff_dB);
    bsc_results{1,3} = sprintf('2-DoF_{n,a} (GoF_{dB} = %.2f)', MetricsBSC(3).diff_dB);
    bsc_results{1,4} = sprintf('2-DoF_{b,a} (GoF_{dB} = %.2f)', MetricsBSC(4).diff_dB);

    colot_gt = '#000000';
    color_1  = '#FF0000';
    color_2  = '#D95319';
    color_3  = '#0072BD';
    color_4  = '#77AC30';

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);

    semilogy(freq, bsc_results{2,1}, '--', 'Color', hex2rgb(color_1), 'LineWidth', line_width, 'DisplayName', bsc_results{1,1});
    hold on;
    semilogy(freq, bsc_results{2,2}, '--', 'Color', hex2rgb(color_2), 'LineWidth', line_width, 'DisplayName', bsc_results{1,2});
    semilogy(freq, bsc_results{2,3}, '--', 'Color', hex2rgb(color_3), 'LineWidth', line_width, 'DisplayName', bsc_results{1,3});
    semilogy(freq, bsc_results{2,4}, '--', 'Color', hex2rgb(color_4), 'LineWidth', line_width, 'DisplayName', bsc_results{1,4});
    semilogy(freq, BSC_gt, '-', 'Color', hex2rgb(colot_gt), 'LineWidth', line_width+0.5, 'DisplayName', 'GT');
    hold off;

    grid on;
    xlabel('Frequency [MHz]', 'FontSize', font_size);
    ylabel('\DeltaBSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
    ylim(ylim_range);
    xlim(xlim_range);
    title('Power Law Model', 'FontSize', font_size + 2);
    legend('Location', 'best', 'FontSize', font_size-7);
    set(gca, 'FontSize', font_size);
end

%% T_combined in same format as your current script
% Metricas a
m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, logical(ones(size(acs_sam))), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, logical(ones(size(acs_sam))), NaN, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{2, 3}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{2, 4}, logical(ones(size(acs_sam))), alpha_sam, '2-DoF-n');

fields = fieldnames(m_3dof);

Ta = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), ...
    'RowNames', fields, 'VariableNames', methods);

% Metrics b dB
m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{3, 3}, logical(ones(size(acs_sam))), NaN, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{3, 4}, logical(ones(size(acs_sam))), pow2db(delta_b_theo), '2-DoF-n');

Tb = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), ...
    'RowNames', fields, 'VariableNames', methods);

% Metrics n
m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, logical(ones(size(acs_sam))), delta_n_theo, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, logical(ones(size(acs_sam))), delta_n_theo, '2-DoF-a');
m_2dofb = get_metrics_homo_gt(maps_results{4, 3}, logical(ones(size(acs_sam))), delta_n_theo, '2-DoF-b');
m_2dofn = get_metrics_homo_gt(maps_results{4, 4}, logical(ones(size(acs_sam))), NaN, '2-DoF-n');

Tn = table(struct2cell(m_3dof), struct2cell(m_2dofa), ...
    struct2cell(m_2dofb), struct2cell(m_2dofn), ...
    'RowNames', fields, 'VariableNames', methods);

clear m_3dof m_2dofa m_2dofb m_2dofn

Ta.Group = repmat("a", height(Ta), 1);
Tb.Group = repmat("b", height(Tb), 1);
Tn.Group = repmat("n", height(Tn), 1);

Ta = movevars(Ta, 'Group', 'Before', Ta.Properties.VariableNames{1});
Tb = movevars(Tb, 'Group', 'Before', Tb.Properties.VariableNames{1});
Tn = movevars(Tn, 'Group', 'Before', Tn.Properties.VariableNames{1});

Ta.Properties.RowNames = strcat(Ta.Properties.RowNames, " a");
Tb.Properties.RowNames = strcat(Tb.Properties.RowNames, " b");
Tn.Properties.RowNames = strcat(Tn.Properties.RowNames, " n");

T_combined = [Ta; Tb; Tn];

%% OPTIONAL SAVE OUTCOMES
if saveOutcomes
    if methodsRegu
        nameExcel = 'metricsRegu_gaussSimu.xlsx';
    else
        nameExcel = 'metricsNoRegu_gaussSimu.xlsx';
    end
    excelFile = fullfile(dirOutcomes, nameExcel);
    writetable(T_combined, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

    if methodsRegu
        nameExcel = 'metricsBSCregu_gaussSimu.xlsx';
    else
        nameExcel = 'metricsBSCnoregu_gaussSimu.xlsx';
    end
    excelFile = fullfile(dirOutcomes, nameExcel);
    writetable(Tbsc, excelFile, 'Sheet', 'Metrics', 'WriteRowNames', true);

    titleFigout = 'simuG_Fig';
    save_all_figures_to_directory(dirOutcomes, titleFigout, 'svg');
end

%% PREPARE RESULT STRUCT
result = struct();

result.mu_a = mu_a;
result.mu_b = mu_b;
result.mu_n = mu_n;
result.mu_setup = mu_setup;

result.methods = methods;
result.label_methods = label_methods;

result.alpha_sam = alpha_sam;
result.alpha_ref = alpha_ref;
result.j_sam = j_sam;
result.j_ref = j_ref;

result.delta_alpha_prior = delta_alpha_prior;
result.delta_b_prior     = delta_b_prior;
result.delta_n_prior     = delta_n_prior;
result.delta_b_theo      = delta_b_theo;
result.delta_n_theo      = delta_n_theo;

result.pars = pars;
result.par_rpl = par_rpl;

result.T_combined = T_combined;
result.Tbsc       = Tbsc;

if saveMapsAndBSC
    result.maps_results = maps_results;
    result.bsc_results  = bsc_results;
end

result.BSC_gt = BSC_gt;
result.freq   = freq;

result.summaryInfo = struct();
result.summaryInfo.mu_a       = mu_a;
result.summaryInfo.mu_b       = mu_b;
result.summaryInfo.mu_n       = mu_n;
result.summaryInfo.log10_mu_a = log10(mu_a);
result.summaryInfo.log10_mu_b = log10(mu_b);
result.summaryInfo.log10_mu_n = log10(mu_n);

end

%% ========================================================================
function M = padconcatenation(C, padval, dim)
max_length = max(cellfun(@numel, C));
M = cellfun(@(x) padarray(x, [max_length - numel(x), 0], padval, 'post'), ...
    C, 'UniformOutput', false);
M = cell2mat(M);
end

%% ========================================================================
function rgb = hex2rgb(hex)
hex = char(hex);
if hex(1) == '#'
    hex(1) = [];
end
rgb = reshape(sscanf(hex, '%2x') / 255, 1, 3);
end