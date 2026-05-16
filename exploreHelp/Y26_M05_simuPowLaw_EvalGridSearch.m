% Script to evaluate the gridSearch generated in LIM cluster and check the
% minimum bias (MPE)

%%
reinit

%% DIRECTORIES

dirCases    = '..\..\data\out\simuPowLaw\gridSearch_cases_2026-04-13_04-59';
files       = dir(fullfile(dirCases, 'grid_*.mat'));

% methods = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

%% %%%%%%%%%%%%% MAE %%%%%%%%%%%%%

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
dirFigOut  = '..\..\data\out\simuPowLaw\';
folderOut  = 'metrics_v3';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, 'MAE_');

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
dirFigOut  = '..\..\data\out\simuPowLaw\';
folderOut  = 'metrics_v3';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, 'NRMSE_');

close all

%% %%%%%%%%%%%%% MPE %%%%%%%%%%%%%

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

%% ------ Save Figures ------
dirFigOut  = '..\..\data\out\simuPowLaw\';
folderOut  = 'metrics_v3';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, 'MPE_');

% close all


%% PLOT 3D ONLY

rangeMetric = [0 0.2];

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


%% SAVE FIGURES ONLY

dirFigOut  = '..\..\data\out\simuPowLaw\';
folderOut  = 'metricsMPE';
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, 'MPE_');


%%
fixedVar = 'mu_a';
fixedVal = 10^2;
plotFixedSlice(struct2table(summary_3dof), fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    '3-DoF', metricName, ...
    rangeMetric, 'turbo');

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
cb = colorbar;  colormap('turbo'); clim(rangeMetric);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

fixedVar = 'mu_a';
fixedVal = 10^2;
plotFixedSlice(T, fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    methodName, metricName, ...
    rangeMetric, 'turbo');

%% ------ 2-DoF-n ------
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
%%
fixedVar = 'mu_a';
fixedVal = power(10, 6);
plotFixedSlice(T, fixedVar, fixedVal, ...
    'mu_b', 'mu_n', ...
    methodName', metricName, ...
    rangeMetric, 'turbo');

%% LOOP 
list_mu = 0:6;

% ------ 3-DoF ------
methodName = '3-DoF';
T = struct2table(summary_3dof);

for ii = 1:length(list_mu)
    fixedVar = 'mu_a';
    fixedVal = power(10, list_mu(ii));
    plotFixedSlice(T, fixedVar, fixedVal, ...
        'mu_b', 'mu_n', ...
        methodName', metricName, ...
        rangeMetric, 'turbo');
end

% ------ 2-DoF-b ------
methodName = '2-DoF-b';
T = struct2table(summary_2dofb);

for ii = 1:length(list_mu)
    fixedVar = 'mu_a';
    fixedVal = power(10, list_mu(ii));
    plotFixedSlice(T, fixedVar, fixedVal, ...
        'mu_b', 'mu_n', ...
        methodName', metricName, ...
        rangeMetric, 'turbo');
end

% ------ 2-DoF-n ------
methodName = '2-DoF-n';
T = struct2table(summary_2dofn);

for ii = 1:length(list_mu)
    fixedVar = 'mu_a';
    fixedVal = power(10, list_mu(ii));
    plotFixedSlice(T, fixedVar, fixedVal, ...
        'mu_b', 'mu_n', ...
        methodName', metricName, ...
        rangeMetric, 'turbo');
end
folderOut = strcat(metricName, '_slices');
fullDirOut = fullfile(dirFigOut, folderOut);
if ~exist(fullDirOut,"dir"); mkdir(fullDirOut); disp('Folder created'); end
save_all_figures_to_directory(fullDirOut, strcat(upper(metricName), '_'));

close all


%% TEST
genPlotsMethodMetric(summary_3dof,  '3-DoF');




%% FAST TEST

% MPE
% summary_3dof    = genSummaryMethod(files, '3-DoF',   metricName);
% summary_2dofn   = genSummaryMethod(files, '2-DoF-n', metricName);

% MAE 
mae_3dof    = genSummaryMethod(files, '3-DoF',   'mae_homo');
mae_2dofn   = genSummaryMethod(files, '2-DoF-n', 'mae_homo');

mae_3dof  = struct2table(mae_3dof);
mae_2dofn = struct2table(mae_2dofn);

T3 = struct2table(summary_3dof);
Tn = struct2table(summary_2dofa);

T3_sorted = sortrows(T3, 'metric_avg', 'ascend');
Tn_sorted = sortrows(Tn, 'metric_avg', 'ascend');

%% ---- Selec fix two mu
targetVal = 10^0.5;
idx       = (mae_2dofn.mu_b == targetVal) & (mae_2dofn.mu_n == targetVal);

mae_2dofn_sel = mae_2dofn(idx,:);

figure, 
subplot(311);
plot(log10(mae_2dofn_sel.mu_a), mae_2dofn_sel.metric_a), grid on
title('a')

subplot(312);
plot(log10(mae_2dofn_sel.mu_a), mae_2dofn_sel.metric_b), grid on
title('b')

subplot(313);
plot(log10(mae_2dofn_sel.mu_a), mae_2dofn_sel.metric_n), grid on
title('n')

%%
targetVal = 10^0.5;
idx       = (mae_3dof.mu_b == targetVal) & (mae_3dof.mu_n == targetVal);
mae_3dof_sel = mae_3dof(idx,:);

figure, 
subplot(311);
% plot((mae_3dof_sel.mu_a), mae_3dof_sel.metric_a), grid on
semilogx((mae_3dof_sel.mu_a), mae_3dof_sel.metric_a), grid on
title('a')

subplot(312);
plot((mae_3dof_sel.mu_a), mae_3dof_sel.metric_b), grid on
title('b')

subplot(313);
plot((mae_3dof_sel.mu_a), mae_3dof_sel.metric_n), grid on
title('n')


%%
mae_3dof(mae_3dof)


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
cb = colorbar;  colormap('turbo'); clim([0 1]);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% Best overall for slices
[~,idxBest] = min(T.metric_avg);
best        = T(idxBest,:);
rangeMetric = [0 0.3];

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
function genPlotsMethodMetric2(summary, methodName)

T           = struct2table(summary);
metricName  = T.metricName{1};

T.log_mu_a  = log10(T.mu_a);
T.log_mu_b  = log10(T.mu_b);
T.log_mu_n  = log10(T.mu_n);

[~,idxBest] = min(T.metric_avg);
best        = T(idxBest,:);

% ------- (1) 3D Scatter -------
figure;
scatter3(T.log_mu_b, T.log_mu_n, T.log_mu_a, ...
    70, T.metric_avg, 'filled');
grid on; colormap('turbo')
xlabel('log_{10}(\mu_b)');
ylabel('log_{10}(\mu_n)');
zlabel('log_{10}(\mu_a)');
cb = colorbar;  colormap('turbo'); clim([0 1]);
cb.Label.String = sprintf('%s', upper(metricName));
title(sprintf('%s, %s', methodName, upper(metricName)));
set(gca, 'FontSize', 12)

% ------- (2) Best 2D slice: mu_b vs mu_a, fixing mu_n -------
fixed_n = best.log_mu_n;
idx     = T.log_mu_n == fixed_n;
x       = T.log_mu_b(idx);
y       = T.log_mu_a(idx);
v       = T.metric_avg(idx);
[X,Y]   = meshgrid(unique(x), unique(y));
V       = griddata(x,y,v,X,Y);

figure;
imagesc(unique(x), unique(y), V);
set(gca,'YDir','normal'); axis image; 
cb = colorbar; colormap('turbo')
clim([0 0.5]);
cb.Label.String = sprintf('%s', upper(metricName));
xlabel('log_{10}(\mu_b)'); ylabel('log_{10}(\mu_a)'); 
title(sprintf('%s, %s, fixed log_{10}(\\mu_n)=%.1f', methodName, upper(metricName), fixed_n));

% ------- (2) Best 2D slice: mu_n vs mu_a, fixing mu_b -------
fixed_b = best.log_mu_b;
idx     = T.log_mu_b == fixed_b;
x       = T.log_mu_n(idx);
y       = T.log_mu_a(idx);
v       = T.metric_avg(idx);
[X,Y]   = meshgrid(unique(x), unique(y));
V       = griddata(x,y,v,X,Y);

figure;
imagesc(unique(x), unique(y), V);
set(gca,'YDir','normal'); axis image;
cb = colorbar;  colormap('turbo')
clim([0 0.5]);
cb.Label.String = sprintf('%s', upper(metricName));
xlabel('log_{10}(\mu_n)'); ylabel('log_{10}(\mu_a)');
title(sprintf('%s, %s, fixed log_{10}(\\mu_b)=%.1f', methodName, upper(metricName), fixed_b));


% ------- (2) Best 2D slice: mu_b vs mu_n, fixing mu_a -------
fixed_a = best.log_mu_a;
idx     = T.log_mu_a == fixed_a;
x       = T.log_mu_b(idx);
y       = T.log_mu_n(idx);
v       = T.metric_avg(idx);
[X,Y]   = meshgrid(unique(x), unique(y));
V       = griddata(x,y,v,X,Y);

figure;
imagesc(unique(x), unique(y), V);
set(gca,'YDir','normal'); axis image;
cb = colorbar; colormap('turbo')
clim([0 0.5]);
cb.Label.String = sprintf('%s', upper(metricName));
xlabel('log_{10}(\mu_b)'); ylabel('log_{10}(\mu_n)');
title(sprintf('%s, %s, fixed log_{10}(\\mu_a)=%.1f', methodName, upper(metricName), fixed_a));

% ------- (3) Individual metric trends -------
% figure;
% hold on;
% plot(T.log_mu_a, abs(T.metric_a), 'ro', 'DisplayName', sprintf('|%s_a|', upper(metricName)));
% plot(T.log_mu_b, abs(T.metric_b), 'go', 'DisplayName', sprintf('|%s_b|', upper(metricName)));
% plot(T.log_mu_n, abs(T.metric_n), 'bo', 'DisplayName', sprintf('|%s_n|', upper(metricName)));
% grid on;
% xlabel('log_{10}(\mu)');
% ylabel(sprintf('%s', upper(metricName)));
% legend('Location','best');
% title(sprintf('Individual %s trends', upper(metricName)));
% 
% ------- (4) Best candidate table -------
T_sorted = sortrows(T, 'metric_avg', 'ascend');

disp(T_sorted(1:10, ...
    {'fileName','mu_a','mu_b','mu_n', ...
     'metric_a','metric_b','metric_n','metric_avg'}))

% Fix a mu and then heat maps of each other mu
end
%%
function Summary = genSummaryMethod(files, methodName, metricName)

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
        S = load(filePath, 'mu_a', 'mu_b', 'mu_n', 'T_combined', 'Tbsc');
    
        T = S.T_combined;
    
        % ---- Summary table ----
        Summary(i).fileName     = string(files(i).name);
    
        Summary(i).mu_a         = S.mu_a;
        Summary(i).mu_b         = S.mu_b;
        Summary(i).mu_n         = S.mu_n;
    
        % ---- extract metric ----
        Summary(i).metricName   = metricName;
        Summary(i).metric_a     = extract_metric(T, methodName, 'a', metricName);
        Summary(i).metric_b     = extract_metric(T, methodName, 'b', metricName);
        Summary(i).metric_n     = extract_metric(T, methodName, 'n', metricName); 
        Summary(i).metric_avg   = mean(abs([ ...
                        Summary(i).metric_a ...
                        Summary(i).metric_b ...
                        Summary(i).metric_n]), ...
                        'omitnan');
        Summary(i).mean_a       = extract_metric(T, methodName, 'a', 'mean');
        Summary(i).std_a        = extract_metric(T, methodName, 'a', 'std');
        Summary(i).mean_b       = extract_metric(T, methodName, 'b', 'mean');
        Summary(i).std_b        = extract_metric(T, methodName, 'b', 'std');
        Summary(i).mean_n       = extract_metric(T, methodName, 'n', 'mean');
        Summary(i).std_n        = extract_metric(T, methodName, 'n', 'std');
    end

    close(hWait)
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



