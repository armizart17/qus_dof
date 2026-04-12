function gridsearch_dof_gauss_save_cases()
% GRIDSEARCH_DOF_GAUSS_SAVE_CASES
% Full grid search over mu_a, mu_g, mu_s for the Gaussian-model DoF simulation.
% Saves each case separately.
%

clear; clc; close all;
warning('off');

%% ========================= USER CONFIG =========================
% Choose the mu grids
% list_mu_a = 10.^(-3:0.5:6);
% list_mu_g = 10.^(-3:0.5:6);
% list_mu_s = 10.^(-3:0.5:6);

% Example smaller test (best):
list_mu_a = 10^4.1;
list_mu_g = [1000] ;
list_mu_s = 10^3;

% Output directory
dirOut = './out/dofJournal/simuGauss/gaussModel/';
if ~exist(dirOut, 'dir')
    mkdir(dirOut);
end

dirCases = fullfile(dirOut, 'gridSearch_cases');
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
plotBoxplots    = false;
plotDeltaBSC    = false;

saveOutcomes    = false;
methodsRegu     = true;
bsc_gt_by_rpm   = true;

% Save map arrays too
saveMapsAndBSC  = true;

% ===============================================================

%% Build all combinations
[MU_A, MU_G, MU_S] = ndgrid(list_mu_a, list_mu_g, list_mu_s);

comb_mu_a = MU_A(:);
comb_mu_g = MU_G(:);
comb_mu_s = MU_S(:);

nComb = numel(comb_mu_a);
fprintf('Total combinations: %d\n', nComb);

fmtExp = @(x) strrep(sprintf('%.1f', log10(x)), '.', 'p');

%% Main loop
for iComb = 1:nComb

    mu_a = comb_mu_a(iComb);
    mu_g = comb_mu_g(iComb);
    mu_s = comb_mu_s(iComb);

    tag_mu_g = fmtExp(mu_g);
    tag_mu_s = fmtExp(mu_s);
    tag_mu_a = fmtExp(mu_a);

    fileName = sprintf('grid_muG_%s_muS_%s_muA_%s.mat', ...
        tag_mu_g, tag_mu_s, tag_mu_a);
    filePath = fullfile(dirCases, fileName);

    fprintf('\n====================================================\n');
    fprintf('Case %d / %d\n', iComb, nComb);
    fprintf('mu_g = 10^(%.2f), mu_s = 10^(%.2f), mu_a = 10^(%.2f)\n', ...
        log10(mu_g), log10(mu_s), log10(mu_a));
    fprintf('%s\n', fileName);

    if skipExisting && exist(filePath, 'file')
        fprintf('Already exists. Skipping.\n');
        continue;
    end

    try
        result = run_one_case_gauss(mu_g, mu_s, mu_a, ...
            plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
            plotBoxplots, plotDeltaBSC, ...
            saveOutcomes, methodsRegu, bsc_gt_by_rpm, saveMapsAndBSC);

        save(filePath, '-struct', 'result', '-v7.3');
        fprintf('Saved case: %s\n', fileName);

    catch ME
        fprintf(2, 'Error in case %d: %s\n', iComb, ME.message);
        errFile = fullfile(dirCases, sprintf('ERROR_%s.mat', fileName(1:end-4)));
        save(errFile, 'mu_a', 'mu_g', 'mu_s', 'ME');
    end
end

fprintf('\nAll done.\n');
fprintf('Cases folder: %s\n', dirCases);

end

%% ========================================================================
function result = run_one_case_gauss(mu_g, mu_s, mu_a, ...
    plotBmode, plotBmodeOverl, plotBSCdB, plotMaps, ...
    plotBoxplots, plotDeltaBSC, ...
    saveOutcomes, methodsRegu, bsc_gt_by_rpm, saveMapsAndBSC)

close all;
warning('off');

%% MU SETUP [mu_g, mu_s, mu_a]
mu_setup = [mu_g, mu_s, mu_a; ... % 2-DoF-a
            mu_g, mu_s, mu_a; ... % 2-DoF-g
            mu_g, mu_s, mu_a; ... % 2-DoF-s
            mu_g, mu_s, mu_a];    % 3-DoF

%% INITIALIZATION
Np2dB           = 20*log10(exp(1));
dB2Np           = 1/Np2dB; %#ok<NASGU>
range_bmode     = [-80 0];

if saveOutcomes
    dirOutcomes = './out/dofJournal/simuGauss/gaussModel';
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
    dirData = '/home/armiz/dataLIM';  % <-- TO CHANGE ^^
else
    error('Unknown operating system');
end


folderDataSam   = 'gaussSimu';
folderDataRef   = 'gaussSimu';

methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};
label_methods = {'3-DoF', '2-DoF_{g,s}', '2-DoF_{s,a}', '2-DoF_{g,a}'};

bsc_results     = cell(4, length(methods));
maps_results    = cell(4, length(methods));

bsc_results(1, :) = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};
maps_results(1, :) = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};

%% LOAD DATA
% SAMPLE
rf_sam_name = 'rf_sd2pcSCALE4_bsc2_gauss_mask7_att0p5';
alpha_sam   = 0.5;
j_sam       = 1.1;

SAM             = load(fullfile(dirData, folderDataSam, rf_sam_name));
SAM.acs         = alpha_sam;
SAM.alpha_power = j_sam;

% REFERENCE
rf_ref_name = 'rf_sd2pcSCALE4_bsc4_att0p4';
alpha_ref   = 0.4;
j_ref       = 1.1;

REF             = load(fullfile(dirData, folderDataRef, rf_ref_name));
REF.acs         = alpha_ref;
REF.alpha_power = j_ref;

%% PRIORS / GT
delta_alpha_prior = alpha_sam - alpha_ref;

% BW = 3.7 - 8.2 MHz values from your script
delta_g_prior = log(db2pow(-0.142481));
delta_s_prior = 0.018717;

delta_g_theo = exp(delta_g_prior);
delta_s_theo = delta_s_prior;

%% BMODE
bmode_sam = db(abs(hilbert(SAM.rf)));
bmode_sam = bmode_sam - max(bmode_sam(:));

bmode_ref = db(abs(hilbert(REF.rf)));
bmode_ref = bmode_ref - max(bmode_ref(:));

%% SPECTRAL PARAMETERS
clear pars
pars.P           = 4096;
pars.overlap     = 0.8;
pars.window_type = 3;
pars.saran_layer = false;

pars.z_roi       = [5 42.5]*1E-3;
pars.x_roi       = [-18 18]*1E-3;

% Final active values from your script
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
    bsc_rpm = BSC.BSCcurve_Uni(:,2); % median
    freq_rpm = BSC.band;

    % Power-law fit just for info
    coeffs   = polyfit(log(freq_rpm), log(bsc_rpm), 1);
    d_n      = coeffs(1);
    ln_db    = coeffs(2);
    d_b      = exp(ln_db);

    fprintf('-----RPM PowLaw (b.(f.^n))-----\n')
    fprintf('Δb           = %f\n', d_b);
    fprintf('b_s/b_r [dB] = %f\n', 10*log10(d_b));
    fprintf('Δn           = %f\n', d_n);
    fprintf('---------\n')

    % Gaussian fit
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

%% GENERAL REGULARIZATION SETTINGS
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
            mu_rpl_tv = mu_setup(1,:);
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
        Z = kron(speye(p*q), -f.^2);

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

    %% 2-DoF-s
    elseif strcmp(estim_method, '2-DoF-s')

        if methodsRegu
            mu_rpl_tv = mu_setup(4,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_ref    = comp_ref_s_bsc(delta_s_prior, band, p, q);
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

        s = delta_s_prior*ones(p*q,1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% 2-DoF-g
    elseif strcmp(estim_method, '2-DoF-g')

        if methodsRegu
            mu_rpl_tv = mu_setup(3,:);
        else
            mu_rpl_tv = [0.001 0.001 0.001];
        end

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [~,p,q] = size(SR);

        comp_ref    = comp_ref_g_bsc(delta_g_prior);
        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);
        SR_comp = SR .* comp_ref .* comp_freq_a;

        f = band(:);
        Y = log(SR_comp);

        Z = kron(speye(p*q), -f.^2);
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

        g = delta_g_prior*ones(p*q,1);

        z = 1E2*repmat(depth,1,q);
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);

        a_Np2dB = Np2dB*Dy*a./dz(:);

    %% 3-DoF
    elseif strcmp(estim_method, '3-DoF')

        if methodsRegu
            mu_rpl_tv = mu_setup(2,:);
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
        Z = kron(speye(p*q), -f.^2);
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
        fprintf('Δg [dB]    : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    else
        fprintf('Δg         : %.3f ± %.4f, %%CV = %.4f\n', round(m_g, 3), round(s_g, 4), round(cv_g, 4));
    end
    fprintf('Δs         : %.4f ± %.4f, %%CV = %.4f\n', round(m_s, 4), round(s_s, 4), round(cv_s, 4));
    fprintf('--------\n');

    %% OPTIONAL MAP PLOTS
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
        title({'$\Delta g$:', ...
            [num2str(round(m_g, 3)), ' $\pm$ ', num2str(round(s_g, 3)), ', CV = ', num2str(round(cv_g, 3))]}, ...
            'Interpreter', 'latex');
        set(gca,'fontsize',fontSize)

        subplot(1,3,3)
        imagesc(Xaxis*cm, Zaxis*cm, s_ratio), colorbar
        axis("image");
        xlabel('Lateral [cm]'), colormap turbo;
        title({'$\Delta s$:', ...
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
        title('$\Delta g$', 'Interpreter','latex')
        set(gca,'fontsize',fontSize)

        nexttile;
        [~,~,hColor] = imOverlayInterp(bmodeFull, s_ratio, range_bmode, [], ...
            transparency, x_img, z_img, roi, xFull, zFull);
        hold on; contour(xFull, zFull, roi, 1,'w--', 'LineWidth', 2); hold off;
        xlabel('Lateral [cm]'), ylabel('Depth [cm]');
        hColor.Label.String = 'a.u.';
        title('$\Delta s$', 'Interpreter','latex')
        set(gca,'fontsize',fontSize)
    end

    %% BSC RECONSTRUCTION
    freq = spectralData_sam.band;
    g_est = median(g_ratio(:));
    s_est = median(s_ratio(:));

    bsc_est_gauss = g_est .* exp(-s_est .* freq.^2);

    bsc_results{2, iMet} = bsc_est_gauss;
    maps_results{2, iMet} = acs_sam;
    maps_results{3, iMet} = g_ratio_dB;
    maps_results{4, iMet} = s_ratio;
end

%% OPTIONAL BOXPLOTS
if plotBoxplots
    font_size  = 34;
    numMethods = size(maps_results, 2);

    acs_data     = cell(1, numMethods);
    g_ratio_data = cell(1, numMethods);
    s_ratio_data = cell(1, numMethods);

    for iMet = 1:numMethods
        acs_data{iMet}     = maps_results{2, iMet}(:);
        g_ratio_data{iMet} = maps_results{3, iMet}(:);
        s_ratio_data{iMet} = maps_results{4, iMet}(:);
    end

    acs_mat     = padconcatenation(acs_data, NaN, 1);
    g_ratio_mat = padconcatenation(g_ratio_data, NaN, 1);
    s_ratio_mat = padconcatenation(s_ratio_data, NaN, 1);

    method_labels = { ...
        '\mathrm{3\textrm{-}DoF}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{g,s}}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{s,a}}', ...
        '\mathrm{2\textrm{-}DoF}_{\mathrm{g,a}}' ...
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
    xt = get(ax, 'XTick');
    if methodsRegu
        ylim([0.05 1.01])
        for i = 1:length(method_labels_a)
            text(xt(i), 0.005*ax.YLim(1), ['$' method_labels_a{i} '$'], ...
                'Interpreter','latex', 'HorizontalAlignment','center', ...
                'FontSize', font_size, 'FontWeight','bold')
        end
    end
    title('\bf\alpha')
    ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
    set(gca, 'FontSize', font_size);

    % g
    g_ratio_mat_filtered = g_ratio_mat(:, [1, 2, 4]);
    method_labels_g      = method_labels([1, 2, 4]);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);
    box on;
    boxplot(g_ratio_mat_filtered, 'Labels', method_labels_g);
    yline(10*log10(delta_g_theo), 'k--')
    ax = gca;
    ax.XTickLabel = {''};
    ax.XTickLabelMode = 'manual';
    xt = get(ax, 'XTick');
    if methodsRegu
        ylim([-1.15 0.15])
        for i = 1:length(method_labels_g)
            text(xt(i), 1.06*ax.YLim(1), ['$' method_labels_g{i} '$'], ...
                'Interpreter','latex', 'HorizontalAlignment','center', ...
                'FontSize', font_size, 'FontWeight','bold')
        end
    end
    title('\bf\Deltag');
    ylabel('\Deltag [dB]');
    set(gca, 'FontSize', font_size);

    % s
    s_ratio_mat_filtered = s_ratio_mat(:, [1, 2, 3]);
    method_labels_s      = method_labels([1, 2, 3]);

    figure;
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 800, 800]);
    box on;
    boxplot(s_ratio_mat_filtered, 'Labels', method_labels_s);
    yline(delta_s_theo, 'k--')
    ax = gca;
    ax.XTickLabel = {''};
    ax.XTickLabelMode = 'manual';
    xt = get(ax, 'XTick');
    if methodsRegu
        ylim([-0.008 0.048])
        for i = 1:length(method_labels_s)
            text(xt(i), 1.35*ax.YLim(1), ['$' method_labels_s{i} '$'], ...
                'Interpreter','latex', 'HorizontalAlignment','center', ...
                'FontSize', font_size, 'FontWeight','bold')
        end
    end
    title('\bf\Deltas');
    ylabel('\Deltas [a.u.]');
    set(gca, 'FontSize', font_size);
end

%% METRICS BSC
BSC_gt = bsc_rpm;
freq = spectralData_sam.band;
diff_fit_dB = @(bsc_pred, bsc_gt) mean(abs(10*log10(bsc_pred) - 10*log10(bsc_gt)));

clear m_3dof m_2dofa m_2dofg m_2dofs MetricsBSC
m_3dof          = get_metrics_homo_gt(bsc_results{2, 1}, true(size(bsc_results{2, 1})), BSC_gt, '3-DoF');
m_3dof.diff_dB  = diff_fit_dB(bsc_results{2, 1}, BSC_gt);
m_3dof.param    = 'BSC';

m_2dofa         = get_metrics_homo_gt(bsc_results{2, 2}, true(size(bsc_results{2, 2})), BSC_gt, '2-DoF-a');
m_2dofa.diff_dB = diff_fit_dB(bsc_results{2, 2}, BSC_gt);
m_2dofa.param   = 'BSC';

m_2dofg         = get_metrics_homo_gt(bsc_results{2, 3}, true(size(bsc_results{2, 3})), BSC_gt, '2-DoF-g');
m_2dofg.diff_dB = diff_fit_dB(bsc_results{2, 3}, BSC_gt);
m_2dofg.param   = 'BSC';

m_2dofs         = get_metrics_homo_gt(bsc_results{2, 4}, true(size(bsc_results{2, 4})), BSC_gt, '2-DoF-s');
m_2dofs.diff_dB = diff_fit_dB(bsc_results{2, 4}, BSC_gt);
m_2dofs.param   = 'BSC';

MetricsBSC(1:4) = [m_3dof; m_2dofa; m_2dofg; m_2dofs];
Tbsc            = struct2table(MetricsBSC);
Tbsc.method     = categorical(Tbsc.method);
Tbsc.param      = categorical(Tbsc.param);

%% OPTIONAL PLOT DELTA BSC
if plotDeltaBSC
    xlim_range = pars.bw + 0.05*[-1 1];
    ylim_range = [0.05 1];
    line_width = 3.5;
    font_size  = 32;

    bsc_results{1,1} = sprintf('3-DoF     (GoF_{dB} = %.2f)', MetricsBSC(1).diff_dB);
    bsc_results{1,2} = sprintf('2-DoF_{g,s} (GoF_{dB} = %.2f)', MetricsBSC(2).diff_dB);
    bsc_results{1,3} = sprintf('2-DoF_{s,a} (GoF_{dB} = %.2f)', MetricsBSC(3).diff_dB);
    bsc_results{1,4} = sprintf('2-DoF_{g,a} (GoF_{dB} = %.2f)', MetricsBSC(4).diff_dB);

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
    title('Gaussian Model', 'FontSize', font_size + 2);
    legend('Location', 'best', 'FontSize', font_size-7);
    set(gca, 'FontSize', font_size);
end

%% METRICS MAPS TABLE - LONG FORMAT
clear m_3dof m_2dofa m_2dofg m_2dofs MetricsParam

%%%%%%%%%%%%%%%%%%% Metrics a %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{2, 1}, true(size(maps_results{2, 1})), alpha_sam, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{2, 2}, true(size(maps_results{2, 1})), NaN,       '2-DoF-a');
m_2dofg = get_metrics_homo_gt(maps_results{2, 3}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-g');
m_2dofs = get_metrics_homo_gt(maps_results{2, 4}, true(size(maps_results{2, 1})), alpha_sam, '2-DoF-s');
m_3dof.param  = 'a';
m_2dofa.param = 'a';
m_2dofg.param = 'a';
m_2dofs.param = 'a';

MetricsParam(1:4) = [m_3dof; m_2dofa; m_2dofg; m_2dofs];
%%%%%%%%%%%%%%%%%%% Metrics a %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metrics g %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{3, 1}, true(size(maps_results{3, 1})), pow2db(delta_g_theo), '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{3, 2}, true(size(maps_results{3, 2})), pow2db(delta_g_theo), '2-DoF-a');
m_2dofg = get_metrics_homo_gt(maps_results{3, 3}, true(size(maps_results{3, 3})), NaN,                  '2-DoF-g');
m_2dofs = get_metrics_homo_gt(maps_results{3, 4}, true(size(maps_results{3, 4})), pow2db(delta_g_theo), '2-DoF-s');
m_3dof.param  = 'g';
m_2dofa.param = 'g';
m_2dofg.param = 'g';
m_2dofs.param = 'g';

MetricsParam(5:8) = [m_3dof; m_2dofa; m_2dofg; m_2dofs];
%%%%%%%%%%%%%%%%%%% Metrics g %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Metrics s %%%%%%%%%%%%%%%%%%%
m_3dof  = get_metrics_homo_gt(maps_results{4, 1}, true(size(maps_results{4, 1})), delta_s_theo, '3-DoF');
m_2dofa = get_metrics_homo_gt(maps_results{4, 2}, true(size(maps_results{4, 2})), delta_s_theo, '2-DoF-a');
m_2dofg = get_metrics_homo_gt(maps_results{4, 3}, true(size(maps_results{4, 3})), delta_s_theo, '2-DoF-g');
m_2dofs = get_metrics_homo_gt(maps_results{4, 4}, true(size(maps_results{4, 4})), NaN,          '2-DoF-s');
m_3dof.param  = 's';
m_2dofa.param = 's';
m_2dofg.param = 's';
m_2dofs.param = 's';

MetricsParam(9:12) = [m_3dof; m_2dofa; m_2dofg; m_2dofs];
%%%%%%%%%%%%%%%%%%% Metrics s %%%%%%%%%%%%%%%%%%%

clear m_3dof m_2dofa m_2dofg m_2dofs

T_combined        = struct2table(MetricsParam);
T_combined.method = categorical(T_combined.method);
T_combined.param  = categorical(T_combined.param);

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
result.mu_g = mu_g;
result.mu_s = mu_s;
result.mu_setup = mu_setup;

result.methods = methods;
result.label_methods = label_methods;

result.alpha_sam = alpha_sam;
result.alpha_ref = alpha_ref;
result.j_sam = j_sam;
result.j_ref = j_ref;

result.delta_alpha_prior = delta_alpha_prior;
result.delta_g_prior     = delta_g_prior;
result.delta_s_prior     = delta_s_prior;
result.delta_g_theo      = delta_g_theo;
result.delta_s_theo      = delta_s_theo;

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
result.summaryInfo.mu_g       = mu_g;
result.summaryInfo.mu_s       = mu_s;
result.summaryInfo.log10_mu_a = log10(mu_a);
result.summaryInfo.log10_mu_g = log10(mu_g);
result.summaryInfo.log10_mu_s = log10(mu_s);

end

%% ========================================================================
function val = get_metric_from_table(T, methodName, paramName, metricName)
val = NaN;

try
    if ~ismember('method', T.Properties.VariableNames), return; end
    if ~ismember('param',  T.Properties.VariableNames), return; end
    if ~ismember(metricName, T.Properties.VariableNames), return; end

    methodCol = string(T.method);
    paramCol  = string(T.param);

    idx = (methodCol == string(methodName)) & (paramCol == string(paramName));
    if ~any(idx), return; end

    tmp = T{find(idx,1,'first'), metricName};

    if iscell(tmp)
        tmp = tmp{1};
    end

    if isnumeric(tmp) && isscalar(tmp)
        val = tmp;
    end
catch
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