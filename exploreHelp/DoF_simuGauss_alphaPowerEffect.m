%% SPECTRAL PARAMETERS
pars.P           = 4096; % NFFT only for calculate BSC_RPM_ok 10wl
pars.bw          = [3.7 8.2]; % [MHz] % last
pars.overlap     = 0.8;
pars.blocksize   = 12; % wavelengths

% new
pars.z_roi       = [5 42.5]*1E-3; % all
pars.x_roi       = [-18 18]*1E-3;
pars.window_type = 3; %  (1) Hanning, (2) Tuckey, (3) Hamming, (4) Tchebychev
pars.saran_layer = false;

%% GENERAL REGULARIZTATION SETTINGS
% Implementation parameters
par_rpl.tol        = 1e-16;
par_rpl.kmax       = 100;
par_rpl.eps_f      = 1e-16;
par_rpl.m_est      = 0; %Robust
par_rpl.ini_tol    = 1e-16;
par_rpl.df_op      = 1;
par_rpl.ini_method = 1; % METHOD LEAST SQUARES INITIALIZATION 

% Parameters for RPL-TV
% mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_g, mu_s, mu_a]

%% UTILS
Np2dB       = 20*log10(exp(1));
dB2Np       = 1/Np2dB;
calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};

pathData = 'D:\emirandaz\qus\data\NonuniformBSC2D_2025\differentAlphaPower';

j_sam_values = [1.1, 1.2, 1.3, 1.4, 1.5];
j_ref_values = [1.1, 1.2, 1.3, 1.4, 1.5];

alpha_sam = 0.5;
alpha_ref = 0.4;

% Preallocate results
num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

%% CALCULATE GROUNDTRUTH RPM METHOD


methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};

clear maps_results_all bsc_results_all

maps_results_all = cell(num_sam, num_ref);
bsc_results_all  = cell(num_sam, num_ref);

tic
% Loop over all j_sam and j_ref combinations
for iSam = 1:num_sam
    for iRef = 1:num_ref
        j_sam = j_sam_values(iSam);
        j_ref = j_ref_values(iRef);

        % Update folder and file names
        folderDataSam = sprintf('');
        % folderDataSam = strrep(folderDataSam, '.', 'p');
        rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
        rf_sam_name = strrep(rf_sam_name, '.', 'p');
        rf_sam_name = strcat(rf_sam_name, '.mat');
        SAM = load(fullfile(pathData, folderDataSam, rf_sam_name));

        folderDataRef = sprintf('');
        % folderDataRef = strrep(folderDataRef, '.', 'p');
        rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
        rf_ref_name = strrep(rf_ref_name, '.', 'p');
        rf_ref_name = strcat(rf_ref_name, '.mat');
        REF = load(fullfile(pathData, folderDataRef, rf_ref_name));

        % Set values
        SAM.alpha_power = j_sam;
        SAM.acs = alpha_sam; % [dB/cm/MHz] 
        REF.alpha_power = j_ref; 
        REF.acs = alpha_ref; % [dB/cm/MHz]

        BSC     = calculateBSC_RPM_fast(SAM, REF, pars); % fast TBD**
        bsc_rpm = BSC.BSCcurve_Uni(:,2); % median

        % Perform linear regression  bsc = d_s . f^2 + ln(d_g) 
        freq     = BSC.band;
        coeffs   = polyfit(-freq.^2, log(bsc_rpm), 1); % Fit y = mx + c
        d_s      = coeffs(1); % Slope = d_s -0.0319 (mean), -0.0317 (median)
        d_g      = exp(coeffs(2)); % Intercept = ln(d_g) 

        % Option 2: BSC = g.exp(-s.f^2) or b.f^n
        % Matrix [log(BSC)] = [ log(f) | 1] * [ n ; log(b) ] rr = MM * qq % POW LAW
        % Matrix [log(BSC)] = [ -f.^2  | 1] * [ s ; log(b) ] rr = MM * qq % GAUSS
        % MM           = [ log(freq), ones( size(freq) ) ]; % POW LAW
        % MM           = [ freq.^2,   ones( size(freq) ) ]; % GAUSS 
        % rr           = log(bsc);
        % qq           = cgs(MM' * MM, MM' * rr, 1e-16, 100); % qq = MM \ rr;
        % d_n          = qq(1);
        % d_b           = exp(qq(2));
       
        % Display results
        
        fprintf('======RPM Gauss (g.exp(-s.f^2))=====\n')
        fprintf('j_sam=%.1f, j_ref=%.1f\n', j_sam, j_ref);
        fprintf('d_g          = %f\n', d_g);
        fprintf('g_s/g_r [dB] = %f\n', 10*log10(d_g));
        fprintf('Δs           = %f\n', d_s);
        fprintf('===================================\n')
        
        bsc_rpm_gauss = d_g*exp(-d_s* freq.^2);

        bsc_results_all{iSam, iRef} = {bsc_rpm, bsc_rpm_gauss, 10*log10(d_g), d_s};
                      
    end
end
tt = toc;
fprintf('Elapsed time %.4f\n', tt)

%%
dirResults = 'D:\emirandaz\qus\qus_dof\qus_dof\out\alpha_power';
if (~exist(dirResults)) mkdir(dirResults); end
fileName = 'resultsRPM_1p1_1p5';

save(fullfile(dirResults, fileName+".mat"), ...
    'bsc_results_all', ...
    'par_rpl', 'j_sam_values', 'j_ref_values', 'alpha_sam', 'alpha_ref');
% save('bsc_maps_results.mat', 'bsc_results_all', 'maps_results_all');

%% DOF
calc2dStats = {@(x) mean(x(:)), @(x) std(x(:)), @(x) 100 * std(x(:)) / mean(x(:))};
methods = {'3-DoF', '2-DoF-a', '2-DoF-g', '2-DoF-s'};

maps_results_dof = cell(num_sam, num_ref, length(methods));
bsc_results_dof  = cell(num_sam, num_ref, length(methods));

% PRIORS
% Compute delta priors
delta_alpha_prior   = alpha_sam - alpha_ref; 
delta_g_prior       = log(1);
delta_s_prior       = 0.0245;

deltaPriorFix = false;

tic
% Loop over all j_sam and j_ref combinations
for iSam = 1:num_sam
    for iRef = 1:num_ref
        j_sam = j_sam_values(iSam);
        j_ref = j_ref_values(iRef);

        % Update folder and file names
        folderDataSam = sprintf('');
        % folderDataSam = strrep(folderDataSam, '.', 'p');
        rf_sam_name = strcat('rfsam_', sprintf('%.3f', j_sam));
        rf_sam_name = strrep(rf_sam_name, '.', 'p');
        rf_sam_name = strcat(rf_sam_name, '.mat');
        SAM         = load(fullfile(pathData, folderDataSam, rf_sam_name));

        folderDataRef = sprintf('');
        % folderDataRef = strrep(folderDataRef, '.', 'p');
        rf_ref_name = strcat('rfref_', sprintf('%.3f', j_ref));
        rf_ref_name = strrep(rf_ref_name, '.', 'p');
        rf_ref_name = strcat(rf_ref_name, '.mat');
        REF         = load(fullfile(pathData, folderDataRef, rf_ref_name));

        % Set values
        SAM.alpha_power = j_sam;
        SAM.acs         = alpha_sam; % [dB/cm/MHz] 
        REF.alpha_power = j_ref; 
        REF.acs         = alpha_ref; % [dB/cm/MHz]

        if ~deltaPriorFix
            delta_g_prior       = log(db2pow(cell2mat(bsc_results_all{iSam, iRef}(3))));
            delta_s_prior       = cell2mat(bsc_results_all{iSam, iRef}(4));
        end

        % Power spectra estimation
        spectralData_sam = calc_powerSpectra_vSimple(SAM, pars);
        S_sam = spectralData_sam.powerSpectra;
        spectralData_ref = calc_powerSpectra_vSimple(REF, pars);
        S_ref = spectralData_ref.powerSpectra;

        % Compute Spectral Ratio
        SR_emz = S_sam ./ S_ref;
        SR = permute(SR_emz, [3,1,2]); clear SR_emz

        band    = spectralData_sam.band;
        depth   = spectralData_sam.depth;
        [r,p,q] = size(SR);
        f = band(:); 

        comp_freq_a = comp_mod_freq_a(alpha_ref,j_sam,j_ref,band,depth,q);

        % Matrix
        X = kron( speye(p*q), ones(size(f)) );
        Z = kron( speye(p*q), -f.^2 );
        W = kron( speye(p*q), -4*f.^j_sam );
        
        dy = 0.5*(diag(ones(p-1,1),1) - diag(ones(p-1,1),-1));
        dy(1,1) = -1; dy(1,2) = 1; dy(end,end) = 1; dy(end,end-1) = -1;
        Dy = sparse(kron(speye(q),dy));
        z = 1E2*repmat(depth,1,q); % 1E2*spectralData_sam.depth * ones(1, q); % 2d array
        dz = reshape(Dy*z(:),p,q);
        dz(end,:) = dz(end-1,:);  

        % Loop over methods
        for iMet = 1:length(methods)
            estim_method = methods{iMet};

            if strcmp(estim_method, '2-DoF-a')
                
                mu_rpl_tv    = [1E3; 1E3; 10^4.1]; % [mu_b, mu_n, mu_a]

                comp_ref    = comp_ref_a(-delta_alpha_prior,j_ref,band,depth,q);                
                SR_comp     = SR .* comp_ref .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_a_prior(Y, X, Z, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_a_prior(Y, X, Z, mu_rpl_tv, u_0, par_rpl);
                
                g = u_opt(1:p*q);
                s = u_opt(p*q+1:2*p*q);
                a_Np2dB = delta_alpha_prior*ones(p*q, 1);

            elseif strcmp( estim_method, '2-DoF-s')

                mu_rpl_tv    = [1E3; 1E3; 10^4.2];

                comp_ref    = comp_ref_s_bsc(delta_s_prior, band, p, q);
                SR_comp     = SR .* comp_ref .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_n_prior(Y, X, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_n_prior(Y, X, W, mu_rpl_tv, u_0, par_rpl);

                g = u_opt(1:p*q);
                a = u_opt(p*q+1:2*p*q);
                s = delta_s_prior*ones(p*q, 1);
                a_Np2dB = Np2dB*Dy*a./dz(:);

            elseif strcmp( estim_method, '2-DoF-g')
                mu_rpl_tv    = [1E3; 1E3; 10^4.2];

                comp_ref    = comp_ref_g_bsc(delta_g_prior);
                SR_comp     = SR .* comp_ref .*comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl_b_prior(Y, Z, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv_b_prior(Y, Z, W, mu_rpl_tv, u_0, par_rpl);
                
                s = u_opt(1:p*q);
                a = u_opt(p*q+1:2*p*q);
                g = delta_g_prior*ones(p*q, 1);
                a_Np2dB = Np2dB*Dy*a./dz(:);

            elseif strcmp(estim_method, '3-DoF')
                mu_rpl_tv    = [1E3; 1E3; 1E4]; % [mu_g, mu_s, mu_a]

                SR_comp = SR .* comp_freq_a;
                Y = log(SR_comp);

                u_0 = initialize_rpl(Y, X, Z, W, mu_rpl_tv, par_rpl);
                [u_opt,~] = rpl_tv(Y, X, Z, W, mu_rpl_tv, u_0, par_rpl);

                g = u_opt(1:p*q);
                s = u_opt(p*q+1:2*p*q);
                a = u_opt(2*p*q+1:3*p*q);
                a_Np2dB = Np2dB*Dy*a./dz(:);
            end

            % Compute final parameters
            g_ratio     = reshape(exp(g), p, q);
            g_ratio_dB  = 10*log10(g_ratio);
            alpha_ratio = reshape(a_Np2dB, p, q);
            s_ratio     = reshape(s, p, q); 
            acs_sam     = alpha_ratio + alpha_ref;

            % Compute statistics
            [m_a, s_a, cv_a] = deal(calc2dStats{1}(acs_sam), calc2dStats{2}(acs_sam), calc2dStats{3}(acs_sam));
            [m_g, s_g, cv_g] = deal(calc2dStats{1}(g_ratio_dB), calc2dStats{2}(g_ratio_dB), calc2dStats{3}(g_ratio_dB));
            [m_s, s_s, cv_s] = deal(calc2dStats{1}(s_ratio), calc2dStats{2}(s_ratio), calc2dStats{3}(s_ratio));

            % Save maps results
            maps_results_dof{iSam, iRef, iMet} = {acs_sam, g_ratio_dB, s_ratio};

            % Compute and save BSC results
            g_est = median(g_ratio(:));
            s_est = median(s_ratio(:));
            bsc_est_gauss = g_est .* exp(-s_est .* band.^2);

            coeffs_pl = polyfit(log(spectralData_sam.band), log(bsc_est_gauss), 1);
            d_n_pl = coeffs_pl(1);
            d_b_pl = exp(coeffs_pl(2));
            bsc_fit_powlaw = d_b_pl * band.^d_n_pl;

            bsc_results_dof{iSam, iRef, iMet} = {bsc_est_gauss, bsc_fit_powlaw};

            fprintf('=============== Method: %s ===============\n', estim_method);
            fprintf('j_sam=%.1f, j_ref=%.1f, α_s    : %.3f ± %.4f\n', j_sam, j_ref, m_a, s_a);
            fprintf('j_sam=%.1f, j_ref=%.1f, Δb [dB]: %.3f ± %.4f\n', j_sam, j_ref, m_g, s_g);
            fprintf('j_sam=%.1f, j_ref=%.1f, Δn     : %.3f ± %.4f\n', j_sam, j_ref, m_s, s_s);
   
        end
    end
end
tt = toc;
fprintf('Elapsed time %.4f\n', tt)
%% Save results

dirResults = 'D:\emirandaz\qus\qus_dof\qus_dof\out\alpha_power';
if (~exist(dirResults)) mkdir(dirResults); end
fileName = 'resultsDoF_mov_1p1_1p5';

save(fullfile(dirResults, fileName+".mat"), ...
    'bsc_results_dof', 'maps_results_dof', ...
    'par_rpl', 'j_sam_values', 'j_ref_values', 'alpha_sam', 'alpha_ref');
% save('bsc_maps_results.mat', 'bsc_results_all', 'maps_results_all');

%%




%% AC METRICS

methods = {'3-DoF', '2-DoF-b', '2-DoF-n'}; % Excluding 2-DoF-a

% Load Data
% pathData = 'C:\Users\armiz\OneDrive\Documentos\MATLAB\qus-lim\out\alpha_power\';

dirResults      = 'D:\emirandaz\qus\qus_dof\qus_dof\out\alpha_power';
pathSaveData    = dirResults;

fileNameRPLTV   = 'resultsDoF_mov_1p1_1p5';
rpltv           = load(fullfile(pathSaveData, fileNameRPLTV));

fileNamefreq    = 'freq';
load(fullfile(pathSaveData, fileNamefreq));

% Extract necessary data
maps_results_all    = rpltv.maps_results_dof;
j_sam_values        = rpltv.j_sam_values;
j_ref_values        = rpltv.j_ref_values;
alpha_sam           = rpltv.alpha_sam;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

% Preallocate metrics storage
metrics_all = cell(num_sam, num_ref, length(methods));

myColormap = 'sky'; % hot, cool, summer, autumn, bone, copper, pink, sky
% myColormap = 'hot';
% myColormap = 'autumn';
font_size = 30;

% Compute metrics for selected methods
tic
for iMethod = 1:length(methods)
    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Extract the a parameter (attenuation coefficient map)
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}))
            if ~isempty(maps_results_all{iSam, iRef, method_index}) 
                a_map = maps_results_all{iSam, iRef, method_index}{1}; 

                % Define a mask (since all values are used, logical ones)
                mask_homo = logical(ones(size(a_map)));

                % Compute metrics
                metrics = get_metrics_homo_gt(a_map, mask_homo, alpha_sam, methods{iMethod});

                % Store metrics
                metrics_all{iSam, iRef, iMethod} = metrics;
            end
        end
    end
end
fprintf('AC Elapsed time: %.2f seconds\n', toc)

%% PLOT AC NRMSE HEATMAP FOR SELECTED METHODS (3-DoF, 2-DoF-b, 2-DoF-n)

% myColormap = 'sky'; % hot, cool, summer, autumn, bone, copper, pink, sky
myColormap = 'hot';
% myColormap = 'autumn';
font_size = 30;

j_sam_values = 1.1:0.1:1.5;
j_ref_values = 1.1:0.1:1.3;

methodsTitle = {sprintf('3-DoF'), sprintf('2-DoF_{a,n}'),  sprintf('2-DoF_{b,a}')};
for iMethod = 1:length(methods)
    metric_name = 'nrmse_homo';  
    metric_matrix = nan(length(j_ref_values), length(j_sam_values));

    % Extract metric values into a matrix
    for iSam = 1:length(j_sam_values)
        for iRef = 1:length(j_ref_values)
            if ~isempty(metrics_all{iSam, iRef, iMethod})
                metric_matrix(iRef, iSam ) = metrics_all{iSam, iRef, iMethod}.(metric_name);
            end
        end
    end

    % Plot the heatmap
    % subplot(1,3,iMethod
    figure,
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 900, 800]); % [x, y, width, height] in pixels
    imagesc(j_sam_values, j_ref_values, flipud(metric_matrix) .* 100); % Convert to percentage
    axis("image")
    colormap(myColormap);
    % clim([0 35])
    hb2 = colorbar; ylabel(hb2, 'NRMSE (%)', 'FontSize', font_size+2);
    
    title(sprintf('AC: %s', methodsTitle{iMethod}), 'FontWeight', 'bold');
    set(gca, 'FontSize', font_size);
    xlabel('$m_s$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    ylabel('$m_r$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    
    set(gca, 'YDir', 'normal'); % Ensures correct orientation
end

%%
% Generate axis labels
xlabels = arrayfun(@(x) sprintf('%.2f', x), j_sam_values, 'UniformOutput', false);
ylabels = arrayfun(@(x) sprintf('%.2f', x), flip(j_ref_values), 'UniformOutput', false); % flipped!

% Create heatmap with flipped matrix
figure,
h = heatmap(xlabels, ylabels, flipud(metric_matrix)*100, ...
            'Colormap', hot, ...
            'ColorbarVisible', 'on');

% Add titles and labels
h.Title = sprintf('AC: %s', methodsTitle{iMethod});
h.XLabel = 'j_{sam}';
h.YLabel = 'j_{ref}';
h.CellLabelFormat = '%.2f';




%% BSC METRICS
diff_fit_dB     = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}; % Excluding 2-DoF-a

fileNameRPM     = 'resultsRPM_1p1_1p5';
rpm             = load(fullfile(pathSaveData, fileNameRPM));

fileNameRPLTV   = 'resultsDoF_mov_1p1_1p5';
rpltv           = load(fullfile(pathSaveData, fileNameRPLTV));

fileNamefreq    = 'freq';
load(fullfile(pathSaveData, fileNamefreq));

% Extract necessary data
bsc_results_all     = rpltv.bsc_results_dof;
j_sam_values        = rpltv.j_sam_values;
j_ref_values        = rpltv.j_ref_values;
alpha_sam           = rpltv.alpha_sam;

bsc_results_all_rpm = rpm.bsc_results_all;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

% Preallocate metrics storage
metrics_all_bsc = cell(num_sam, num_ref, length(methods));

clear metrics

% Compute metrics for selected methods
tic
for iMethod = 1:length(methods)
    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Extract the a parameter (attenuation coefficient map)
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}));
            if ~isempty(bsc_results_all{iSam, iRef, method_index}) 
                bsc_rpl = bsc_results_all{iSam, iRef, method_index}{1}; 

                bsc_rpm = bsc_results_all_rpm{iSam, iRef}{1}; 

                % Define a mask (since all values are used, logical ones)
                mask_homo = logical(ones(size(bsc_rpl)));

                % Compute metrics
                metrics          = get_metrics_homo_gt(bsc_rpl, mask_homo, bsc_rpm, methods{iMethod});
                metrics.diff_dB  = diff_fit_dB(bsc_rpl, bsc_rpm);

                % Store metrics
                metrics_all_bsc{iSam, iRef, iMethod} = metrics;
            end
        end
    end
end
fprintf('BSC Elapsed time: %.2f seconds\n', toc)
%% PLOT SEMILOGY BSC
% Number of methods
nMethods = length(methods);

% Extract necessary data
bsc_results_all     = rpltv.bsc_results_dof;
j_sam_values        = rpltv.j_sam_values;
j_ref_values        = rpltv.j_ref_values;
alpha_sam           = rpltv.alpha_sam;

bsc_results_all_rpm = rpm.bsc_results_all;

freq = band;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

pars.bw    = [3.7 8.2]; % [MHz] % last
xlim_range = pars.bw + 0.05*[-1 1]; % X-axis limits
ylim_range = [0.05 1]; % Y-axis limits

methodsTitle = {sprintf('3-DoF'), sprintf('2-DoF_{a,n}'),  sprintf('2-DoF_{b,a}')};

% Loop over methods
for iMethod = 1:nMethods
    figure('Name', ['BSC Comparison - Method: ', methods{iMethod}], 'Color', 'w');

    % Define subplot grid
    nRows = num_sam;
    nCols = num_ref;

    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Subplot index
            idx = (iSam - 1) * num_ref + iRef;

            clear bsc_rpl bsc_rpm

            % Extract BSC estimates
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}));
            bsc_rpl = bsc_results_all{iSam, iRef, method_index};
            bsc_rpm = bsc_results_all_rpm{iSam, iRef};

            % Continue only if both exist
            %if isempty(bsc_rpl) || isempty(bsc_rpm)
             %   continue;
            %end

            % Extract BSC curves (assume first cell contains freq vs bsc)
            bsc_rpl = bsc_rpl{1};  % 1D array same size as freq
            bsc_rpm = bsc_rpm{1};

            % Plot in subplot
            subplot(nRows, nCols, idx);
            semilogy(band, bsc_rpl, 'b-', 'LineWidth', 1.5); hold on;
            semilogy(band, bsc_rpm, 'k-', 'LineWidth', 1.5);
            grid on;
            title(sprintf('j_{sam}=%.2f, j_{ref}=%.2f', j_sam_values(iSam), j_ref_values(iRef)), 'FontSize', 8);
            if iSam == nRows
                xlabel('Frequency (MHz)', 'FontSize', 8);
            end
            if iRef == 1
                ylabel('BSC (1/(sr·cm))', 'FontSize', 8);
            end
             ylim_range = [0.05 1]; % Y-axis limits
             ylim(ylim_range)
             xlim(xlim_range)
            set(gca, 'FontSize', 7);
        end
    end

    % Add legend outside subplots
    legend({'RPLTV', 'RPM'}, 'Position', [0.9, 0.93, 0.1, 0.05]);
    sgtitle(['BSC Comparison - ', methods{iMethod}], 'FontSize', 12);
end


%% PLOT NRMSE BSC HEATMAP FOR SELECTED METHODS (3-DoF, 2-DoF-b, 2-DoF-n)

j_sam_values = 1.1:0.1:1.5;
j_ref_values = j_sam_values;

% methodsTitle = {sprintf('3-DoF'), sprintf('2-DoF_{a,n}'),  sprintf('2-DoF_{b,a}')};

for iMethod = 1:length(methods)
    %metric_name = 'nrmse_homo'; 
    metric_name = 'diff_dB'; 
    metric_matrix = nan(length(j_sam_values), length(j_ref_values));

    % Extract metric values into a matrix
    for iSam = 1:length(j_sam_values)
        for iRef = 1:length(j_ref_values)
            if ~isempty(metrics_all_bsc{iSam, iRef, iMethod})
                metric_matrix(iSam, iRef) = metrics_all_bsc{iSam, iRef, iMethod}.(metric_name);
            end
        end
    end

    % Plot the heatmap
    figure,
    set(gcf, 'Units', 'pixels', 'Position', [100, 100, 920, 800]); % [x, y, width, height] in pixels
    imagesc(j_sam_values, j_ref_values, metric_matrix); 
    axis("image")
    colormap(myColormap);
    %clim([0 5])
    hb2 = colorbar; ylabel(hb2, 'NRMSE (%)', 'FontSize', font_size+2);
    hb2 = colorbar; ylabel(hb2, 'GoF (dB)', 'FontSize', font_size+2);
    
    title(sprintf('BSC: %s', methods{iMethod}), 'FontWeight', 'bold');
    set(gca, 'FontSize', font_size);
    xlabel('$m_s$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    ylabel('$m_r$', 'Interpreter', 'latex', 'FontSize', font_size+15);
    set(gca, 'YDir', 'normal'); % Ensures correct orientation
    
end

%% SAVE FIGURES

keyboard
dirFigout = '.\TUFFC25\simuGauss\ac0p5\gaussianmodel\alpha_powerv2';
if (~exist(dirFigout)); mkdir (dirFigout); end
titleFigout = 'Fig';
save_all_figures_to_directory(dirFigout, titleFigout, 'svg')

%%
% Number of methods
nMethods = length(methods);

diff_fit_dB     = @(bsc_pred, bsc_gt) mean ( abs ( 10*log10(bsc_pred) - 10*log10(bsc_gt) ) );

methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}; % Excluding 2-DoF-a

fileNameRPM     = 'resultsRPM_1p1_1p5';
rpm             = load(fullfile(pathSaveData, fileNameRPM));

fileNameRPLTV   = 'resultsDoF_mov_1p1_1p5';
rpltv           = load(fullfile(pathSaveData, fileNameRPLTV));

fileNamefreq    = 'freq';
load(fullfile(pathSaveData, fileNamefreq));

% Extract necessary data
bsc_results_all     = rpltv.bsc_results_dof;
j_sam_values        = rpltv.j_sam_values;
j_ref_values        = rpltv.j_ref_values;
alpha_sam           = rpltv.alpha_sam;

bsc_results_all_rpm = rpm.bsc_results_all;

freq = band;

num_sam = length(j_sam_values);
num_ref = length(j_ref_values);

% Loop over methods
for iMethod = 1:nMethods
    figure('Name', ['BSC Comparison - Method: ', methods{iMethod}], 'Color', 'w');

    % Define subplot grid
    nRows = num_sam;
    nCols = num_ref;

    for iSam = 1:num_sam
        for iRef = 1:num_ref
            % Subplot index
            idx = (iSam - 1) * num_ref + iRef;

            % Extract BSC estimates
            method_index = find(strcmp({'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'}, methods{iMethod}));
            bsc_rpl = bsc_results_all{iSam, iRef, method_index};
            bsc_rpm = bsc_results_all_rpm{iSam, iRef};

            % Continue only if both exist
            if isempty(bsc_rpl) || isempty(bsc_rpm)
                continue;
            end

            % Extract BSC curves (assume first cell contains freq vs bsc)
            bsc_rpl = bsc_rpl{1};  % 1D array same size as freq
            bsc_rpm = bsc_rpm{1};

            % Plot in subplot
            subplot(nRows, nCols, idx);
            plot(freq, bsc_rpl, 'b-', 'LineWidth', 1.5); hold on;
            plot(freq, bsc_rpm, 'r--', 'LineWidth', 1.5);
            grid on;
            title(sprintf('j_{sam}=%.2f, j_{ref}=%.2f', j_sam_values(iSam), j_ref_values(iRef)), 'FontSize', 8);
            if iSam == nRows
                xlabel('Frequency (MHz)', 'FontSize', 8);
            end
            if iRef == 1
                ylabel('BSC (1/(sr·cm))', 'FontSize', 8);
            end
            set(gca, 'FontSize', 7);
        end
    end

    % Add legend outside subplots
    legend({'RPLTV', 'RPM'}, 'Position', [0.9, 0.93, 0.1, 0.05]);
    sgtitle(['BSC Comparison - ', methods{iMethod}], 'FontSize', 12);
end
