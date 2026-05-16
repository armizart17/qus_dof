% Script to evaluate the gridSearch generated in LIM cluster and check the
% minimum bias (MPE)

%%
reinit

%% DIRECTORIES

dirCases    = '..\..\data\out\simuGauss\gaussModel\gridSearch_cases_2026-04-13_07-44'; % Gauss- Gauss
files       = dir(fullfile(dirCases, 'grid_*.mat'));

% methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

%%


%% PLOT 3D ONLY

rangeMetric = [0 2];

% ------ 3-DoF ------
methodName = '3-DoF';
T = struct2table(summary_3dof);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)


rangeMetric = [0 2];

% ------ 2-DoF-a ------
methodName = '2-DoF-a';
T = struct2table(summary_2dofa);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% ------ 2-DoF-n ------
methodName = '2-DoF-n';
T = struct2table(summary_2dofn);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% %%%%%%%%%%%%% MAE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'mae';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethod(files, '3-DoF',   metricName);
summary_2dofa   = genSummaryMethod(files, '2-DoF-a', metricName);
summary_2dofb   = genSummaryMethod(files, '2-DoF-b', metricName);
summary_2dofn   = genSummaryMethod(files, '2-DoF-n', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini)); clear tini

% ------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');

% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\gaussModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all


%% %%%%%%%%%%%%% NRMSE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'nrmse';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethod(files, '3-DoF',   metricName);
summary_2dofa   = genSummaryMethod(files, '2-DoF-a', metricName);
summary_2dofb   = genSummaryMethod(files, '2-DoF-b', metricName);
summary_2dofn   = genSummaryMethod(files, '2-DoF-n', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini));

% ------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');

% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\gaussModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all

%% %%%%%%%%%%%%% MPE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'mpe';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethodGG(files, '3-DoF',   metricName);
summary_2dofa   = genSummaryMethodGG(files, '2-DoF-a', metricName);

%%
summary_2dofb   = genSummaryMethodGG(files, '2-DoF-g', metricName);
summary_2dofn   = genSummaryMethodGG(files, '2-DoF-s', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini));

%% ------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');

%% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\gaussModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all


%%
dirCases    =  '..\..\data\out\simuGauss\powLawModel\gridSearch_cases_2026-04-14_08-12'; % Gauss- PowLaw

files       = dir(fullfile(dirCases, 'grid_*.mat'));

% %%%%%%%%%%%%% NRMSE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'nrmse';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethodGP(files, '3-DoF',   metricName);
%%
summary_2dofa   = genSummaryMethodGP(files, '2-DoF-a', metricName);
summary_2dofb   = genSummaryMethodGP(files, '2-DoF-b', metricName);
summary_2dofn   = genSummaryMethodGP(files, '2-DoF-n', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini));

%% %------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');


%%

%% PLOT 3D ONLY

rangeMetric = [0 2];

% ------ 3-DoF ------
methodName = '3-DoF';
T = struct2table(summary_3dof);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); 
clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% ------ 2-DoF-b ------
methodName = '2-DoF-b';
T = struct2table(summary_2dofb);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); 
% clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% ------ 2-DoF-n ------
methodName = '2-DoF-n';
T = struct2table(summary_2dofn);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); 
% clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% ------ 2-DoF-a ------
methodName = '2-DoF-a';
T = struct2table(summary_2dofa);

figure;
scatter3(log10(T.mu_b), log10(T.mu_n), log10(T.mu_a), ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); 
% clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

%% PLOT SLICES ONLY
fixedVar = 'mu_a';
% fixedVal = 10^-3; rangeMetric = [0 0.15];
% fixedVal = 10^3; rangeMetric = [0 0.15];
fixedVal = 10^6; rangeMetric = [0 0.15];

plotFixedSlice(struct2table(summary_3dof), fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    '3-DoF', metricName, ...
    rangeMetric, 'turbo');


plotFixedSlice(struct2table(summary_2dofa), fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    '2-DoF-a', metricName, ...
    rangeMetric, 'turbo');

plotFixedSlice(struct2table(summary_2dofb), fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    '2-DoF-b', metricName, ...
    rangeMetric, 'turbo');

plotFixedSlice(struct2table(summary_2dofn), fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    '2-DoF-n', metricName, ...
    rangeMetric, 'turbo');


%%
% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\powLawModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all

% %%%%%%%%%%%%% MAE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'mae';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethod(files, '3-DoF',   metricName);
summary_2dofa   = genSummaryMethod(files, '2-DoF-a', metricName);
summary_2dofb   = genSummaryMethod(files, '2-DoF-b', metricName);
summary_2dofn   = genSummaryMethod(files, '2-DoF-n', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini)); clear tini

% ------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');

% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\gaussModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all

% %%%%%%%%%%%%% MPE %%%%%%%%%%%%%

% ------ Evaluated metric ------
metricName  = 'mpe';

% ------ Read .mat files and generate summary ------
tini = tic;
summary_3dof    = genSummaryMethod(files, '3-DoF',   metricName);
summary_2dofa   = genSummaryMethod(files, '2-DoF-a', metricName);
summary_2dofb   = genSummaryMethod(files, '2-DoF-b', metricName);
summary_2dofn   = genSummaryMethod(files, '2-DoF-n', metricName);
fprintf('Time all methods: %.2f secs \n', toc(tini));

% ------ Generate plots ------
genPlotsMethodMetric(summary_3dof,  '3-DoF');
genPlotsMethodMetric(summary_2dofa, '2-DoF-a');
genPlotsMethodMetric(summary_2dofb, '2-DoF-b');
genPlotsMethodMetric(summary_2dofn, '2-DoF-n');

% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuGauss\powLawModel';
folderOut  = 'metrics_v1';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all



%%


function Summary = genSummaryMethodGG(files, methodName, metricName)

    nFiles      = numel(files);
    
    Summary(nFiles,1) = struct( ...
        'fileName', "", ...
        'mu_a', NaN, ...
        'mu_b', NaN, ...
        'mu_n', NaN, ...
        'metricName', NaN, ...
        'metric_a', NaN, ...
        'metric_b', NaN, ...
        'metric_n', NaN, ...
        'metric_avg', NaN, ...
        'mean_a', NaN, ...
        'std_a', NaN, ...
        'mean_b', NaN, ...
        'std_b', NaN, ...
        'mean_n', NaN, ...
        'std_n', NaN );

    hWait = waitbar(0, sprintf('Summary of %s - %s. Processing 0/%d', methodName, metricName, nFiles));
       
    for i = 1:nFiles

        % ---- Update waitbar ----  
        waitbar(i / nFiles, hWait, sprintf('Summary of %s - %s. Processing %d/%d', methodName, metricName, i, nFiles));
    
        filePath = fullfile(files(i).folder, files(i).name);
        S = load(filePath, 'mu_a', 'mu_b', 'mu_n', 'T_combined', 'Tbsc', ...
                                   'mu_g', 'mu_s');
    
        T = S.T_combined;
    
        % ---- Summary table ----
        Summary(i).fileName     = string(files(i).name);
    
        Summary(i).mu_a         = S.mu_a;
        % Summary(i).mu_b         = S.mu_b;
        % Summary(i).mu_n         = S.mu_n;
        if isfield(S,'mu_b')
            Summary(i).mu_b = S.mu_b;
            varB = 'b';
        elseif isfield(S,'mu_g')
            Summary(i).mu_b = S.mu_g;
            varB = 'g';
        else
            Summary(i).mu_b = [];
        end
        
        if isfield(S,'mu_n')
            Summary(i).mu_n = S.mu_n;
            varN = 'n';
        elseif isfield(S,'mu_s')
            Summary(i).mu_n = S.mu_s;
            varN = 's';
        else
            Summary(i).mu_n = [];
        end
    
        % ---- extract metric ----
        Summary(i).metricName   = metricName;
        Summary(i).metric_a     = extract_metric(T, methodName, 'a', metricName);
        Summary(i).metric_b     = extract_metric(T, methodName, varB, metricName);
        Summary(i).metric_n     = extract_metric(T, methodName, varN, metricName); 
        Summary(i).metric_avg   = mean(abs([ ...
                        Summary(i).metric_a ...
                        Summary(i).metric_b ...
                        Summary(i).metric_n]), ...
                        'omitnan');
        Summary(i).mean_a       = extract_metric(T, methodName, 'a', 'mean');
        Summary(i).std_a        = extract_metric(T, methodName, 'a', 'std');
        Summary(i).mean_b       = extract_metric(T, methodName, varB, 'mean');
        Summary(i).std_b        = extract_metric(T, methodName, varB, 'std');
        Summary(i).mean_n       = extract_metric(T, methodName, varN, 'mean');
        Summary(i).std_n        = extract_metric(T, methodName, varN, 'std');
    end

    close(hWait)
end

function Summary = genSummaryMethodGP(files, methodName, metricName)

    nFiles      = numel(files);
    
    Summary(nFiles,1) = struct( ...
        'fileName', "", ...
        'mu_a', NaN, ...
        'mu_b', NaN, ...
        'mu_n', NaN, ...
        'metricName', NaN, ...
        'metric_a', NaN, ...
        'metric_b', NaN, ...
        'metric_n', NaN, ...
        'metric_avg', NaN, ...
        'mean_a', NaN, ...
        'std_a', NaN, ...
        'mean_b', NaN, ...
        'std_b', NaN, ...
        'mean_n', NaN, ...
        'std_n', NaN );

    hWait = waitbar(0, sprintf('Summary of %s - %s. Processing 0/%d', methodName, metricName, nFiles));
       
    for i = 1:nFiles

        % ---- Update waitbar ----  
        waitbar(i / nFiles, hWait, sprintf('Summary of %s - %s. Processing %d/%d', methodName, metricName, i, nFiles));
    
        filePath = fullfile(files(i).folder, files(i).name);
        S = load(filePath, 'mu_a', 'mu_b', 'mu_n', 'T_combined', 'Tbsc', ...
                                   'mu_g', 'mu_s');
    
        T = S.T_combined;
    
        % ---- Summary table ----
        Summary(i).fileName     = string(files(i).name);
    
        Summary(i).mu_a         = S.mu_a;
        % Summary(i).mu_b         = S.mu_b;
        % Summary(i).mu_n         = S.mu_n;
        if isfield(S,'mu_b')
            Summary(i).mu_b = S.mu_b;
            varB = 'b';
        elseif isfield(S,'mu_g')
            Summary(i).mu_b = S.mu_g;
            varB = 'g';
        else
            Summary(i).mu_b = [];
        end
        
        if isfield(S,'mu_n')
            Summary(i).mu_n = S.mu_n;
            varN = 'n';
        elseif isfield(S,'mu_s')
            Summary(i).mu_n = S.mu_s;
            varN = 's';
        else
            Summary(i).mu_n = [];
        end
    
        % ---- extract metric ----
        Summary(i).metricName   = metricName;
        Summary(i).metric_a     = extract_metric_wide(T, methodName, 'a', metricName);
        Summary(i).metric_b     = extract_metric_wide(T, methodName, varB, metricName);
        Summary(i).metric_n     = extract_metric_wide(T, methodName, varN, metricName); 
        Summary(i).metric_avg   = mean(abs([ ...
                        Summary(i).metric_a ...
                        Summary(i).metric_b ...
                        Summary(i).metric_n]), ...
                        'omitnan');
        Summary(i).mean_a       = extract_metric_wide(T, methodName, 'a', 'mean');
        Summary(i).std_a        = extract_metric_wide(T, methodName, 'a', 'std');
        Summary(i).mean_b       = extract_metric_wide(T, methodName, varB, 'mean');
        Summary(i).std_b        = extract_metric_wide(T, methodName, varB, 'std');
        Summary(i).mean_n       = extract_metric_wide(T, methodName, varN, 'mean');
        Summary(i).std_n        = extract_metric_wide(T, methodName, varN, 'std');
    end

    close(hWait)
end


function val = extract_metric_wide(T, methodName, paramName, metricName)
% Extract scalar metric from wide table with row names.
%
% Example:
%   mpe_a = extract_metric_wide(T, '3-DoF', 'a', 'mpe');
%   mae_b = extract_metric_wide(T, '2-DoF-b', 'b', 'mae');

val = NaN;

if isempty(T) || ~istable(T)
    return;
end

metricRowName = sprintf('%s_homo %s', metricName, paramName);

% Check row exists
idxRow = strcmp(string(T.Properties.RowNames), string(metricRowName));

if ~any(idxRow)
    warning('Row "%s" not found.', metricRowName);
    return;
end

% Check method column exists
if ~ismember(methodName, T.Properties.VariableNames)
    warning('Method column "%s" not found.', methodName);
    return;
end

tmp = T{idxRow, methodName};

% Your table values are inside cells: {[0.2150]}
if iscell(tmp)
    tmp = tmp{1};
end

if isnumeric(tmp) && isscalar(tmp)
    val = tmp;
end

end
%%
function val = extract_metric(T, methodName, paramName, metricName)
% Extract any scalar metric from T_combined.
%
% Example:
%   mpe_a = extract_metric(T, '3-DoF', 'a', 'mpe');
%   mae_b = extract_metric(T, '3-DoF', 'b', 'mae');

val = NaN;

metricName   = strcat(metricName, '_homo'); % just because the metrics are mentioned like that 

if isempty(T) || ~istable(T);  return; end

requiredCols = {'method','param',metricName};


if ~all(ismember(requiredCols, T.Properties.VariableNames)); return; end

idx = string(T.method) == string(methodName) & ...
      string(T.param)  == string(paramName);

if ~any(idx);    return; end

tmp = T{find(idx,1,'first'), metricName};

if iscell(tmp);  tmp = tmp{1}; end

if isnumeric(tmp) && isscalar(tmp);    val = tmp; end

end


%% NEW FUNCTIONS
function genPlotsMethodMetric(summary, methodName)

T           = struct2table(summary);
metricName  = T.metricName{1};

T.log_mu_a  = log10(T.mu_a);
T.log_mu_b  = log10(T.mu_b);
T.log_mu_n  = log10(T.mu_n);

% ------- (1) 3D Scatter -------
figure;
scatter3(T.log_mu_b, T.log_mu_n, T.log_mu_a, ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); clim([0 2]);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% Best overall for slices
[~,idxBest] = min(T.metric_avg);
best        = T(idxBest,:);
rangeMetric = [0 2];

% ------- (2) Best 2D slice: mu_b vs mu_a, fixing mu_n -------
fixedVar = 'mu_n';
fixedVal = best.mu_n; % could be fixed
plotFixedSlice(T, fixedVar, fixedVal, ...
    'mu_b', 'mu_a', ... % x, y
    methodName, metricName, ...
    rangeMetric, 'turbo');

% ------- (2) Best 2D slice: mu_n vs mu_a, fixing mu_b -------
fixedVar = 'mu_b';
fixedVal = best.mu_b; % could be fixed
plotFixedSlice(T, fixedVar, fixedVal, ...
    'mu_n', 'mu_a', ... % x, y
    methodName, metricName, ...
    rangeMetric, 'turbo');

% ------- (2) Best 2D slice: mu_b vs mu_n, fixing mu_a -------
fixedVar = 'mu_a';
fixedVal = best.mu_a; % could be fixed
plotFixedSlice(T, fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ... % x, y
    methodName, metricName, ...
    rangeMetric, 'turbo');

% ------- Best candidate table -------
T_sorted = sortrows(T, 'metric_avg', 'ascend');

disp(T_sorted(1:10, ...
    {'fileName','mu_a','mu_b','mu_n', ...
     'metric_a','metric_b','metric_n','metric_avg'}))

end
%%
function plotFixedSlice(T, fixedVar, fixedVal, xVar, yVar, methodName, metricName, climVals, cmapName)

% Floating-point safe comparison
tol = 1e-12;
idx = abs(T.(fixedVar) - fixedVal) < tol;

x = T.(xVar)(idx);
y = T.(yVar)(idx);
v = T.metric_avg(idx);

xUnique = log10(unique(x));
yUnique = log10(unique(y));

V = reshape(v, [numel(yUnique), numel(xUnique)]);

figure;
imagesc(xUnique, yUnique, V);
set(gca, 'YDir', 'normal');
axis image;
colormap(cmapName);
clim(climVals);
cb = colorbar;
cb.Label.String = sprintf('%s', upper(metricName));
xlabel(formatLabel(xVar));
ylabel(formatLabel(yVar));
title(sprintf('%s, %s, fixed %s = %.1f', ...
    methodName, upper(metricName), formatLabel(fixedVar), log10(fixedVal)));
set(gca, 'FontSize', 12)

end

%
function labelStr = formatLabel(varName)

switch varName
    case 'mu_a'
        labelStr = 'log_{10}(\mu_a)';
    case 'mu_b'
        labelStr = 'log_{10}(\mu_b)';
    case 'mu_n'
        labelStr = 'log_{10}(\mu_n)';
    otherwise
        labelStr = varName;
end
end

