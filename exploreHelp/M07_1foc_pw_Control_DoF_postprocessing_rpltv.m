% ====================================================================== %
% Post processing of Liver Healthy volunteers 
% BF Data (Acq in 2025 in Avendano CONTROL JULY 10), reference Wisconsin
% Processed by DoF method
% Post-processing all cases dispersion boxcharts, violins, etc.
% ====================================================================== %
%%
% init
warning('off'); %% Important to turn warns in loop off
%% DIRECTORIES

baseDir     = 'D:\emirandaz\qus\data\liver\bf_SavedDataCurvilinearPW\';

resultsDir  = fullfile(baseDir,'resultsControl_RPL_lu'); 

metric_name = 'std';  % Change this to 'mean', 'std', 'median' or 'iqr' as needed

% Result files (Processed by DoF method)
resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles      = length(resultFiles);
%% METHOD LABELS

methods      = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

nMethods = length(methods);

varNames = {'num',  'method','mean',  'std',   'median','iqr',   'patient','sample'};
varTypes = {'int16','string','double','double','double','double','int16',  'int16'};
T       = table('Size',[nFiles*(nMethods-1), length(varNames)], ...
            'VariableTypes',varTypes,'VariableNames',varNames);

Ta = T;
Tb = T;
Tn = T;
clear T


%% ================== alpha metrics ==================

indices_alpha = [1, 3, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    patient     = str2double(fileName(1:3));
    sample      = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_alpha(jj);
        img_map = maps_results_dof{iMethod}.alpha;

        Ta(row,:) = {row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    patient, sample};
        row = row + 1;
    end
end
% writetable(Ta,fullfile(resultsDir,'meanACS.xlsx'))

% ================== Delta b metrics ==================

indices_b = [1, 2, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    patient     = str2double(fileName(1:3));
    sample      = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_b(jj);
        img_map = maps_results_dof{iMethod}.b_dB;

        Tb(row,:) = {row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    patient, sample};
        row = row + 1;
    end
end

% ================== Delta n metrics ==================

indices_n = [1, 2, 3];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    patient     = str2double(fileName(1:3));
    sample      = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_n(jj);
        img_map = maps_results_dof{iMethod}.n;

        Tn(row,:) = {row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    patient, sample};
        row = row + 1;
    end
end


%% BOX PLOTS
% =============== Box alpha ===============

methods_unique = unique(Ta.method, 'stable');  % Keep original order
num_methods = length(methods_unique);

% Prepare data
data_all  = [];
group_all = [];

for i = 1:num_methods
    this_method = methods_unique(i);
    this_data   = Ta{Ta.method == this_method, metric_name};  % dynamic field access

    data_all  = [data_all; this_data];
    group_all = [group_all; repmat(this_method, length(this_data), 1)];
end

% Create box plot
figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 200, 500, 400]); % [x, y, width, height]
boxplot(data_all, group_all);
xlabel('Reconstruction Method');
ylabel('\alpha [dB/cm/MHz]' );
title(sprintf('%s Dispersion of \\alpha per Method', upper(metric_name)));
grid minor;

% Plot using boxchart
% figure;
% boxchart(categorical(group_all, methods_unique), data_all);
% xlabel('Reconstruction Method');
% ylabel('\alpha [dB/cm/MHz]' );
% title(sprintf('%s Dispersion of \\alpha per Method', upper(metric_name)));
% grid on;
%
% =============== Box Delta b ===============

methods_unique = unique(Tb.method, 'stable');  % Keep original order
num_methods = length(methods_unique);

% Prepare data
data_all  = [];
group_all = [];

for i = 1:num_methods
    this_method = methods_unique(i);
    this_data   = Tb{Tb.method == this_method, metric_name};  % dynamic field access

    data_all  = [data_all; this_data];
    group_all = [group_all; repmat(this_method, length(this_data), 1)];
end

% % Create box plot
figure;
set(gcf, 'Units', 'pixels', 'Position', [650, 200, 500, 400]); % [x, y, width, height]
boxplot(data_all, group_all);
xlabel('Reconstruction Method');
ylabel('\Deltab [dB]' );
title(sprintf('%s Dispersion of \\Deltab per Method', upper(metric_name)));
grid minor;

% Plot using boxchart
% figure;
% boxchart(categorical(group_all, methods_unique), data_all);
% xlabel('Reconstruction Method');
% ylabel('\Deltrab [dB]' );
% title(sprintf('%s Dispersion of \\Deltab per Method', upper(metric_name)));
% grid on;

% =============== Box Delta n ===============


methods_unique = unique(Tn.method, 'stable');  % Keep original order
num_methods = length(methods_unique);

% Prepare data and group labels
data_all    = [];
group_all   = [];

for i = 1:num_methods
    this_method = methods_unique(i);
    this_data   = Tn{Tn.method == this_method, metric_name};  % dynamic field access

    data_all    = [data_all; this_data];
    group_all   = [group_all; repmat(this_method, length(this_data), 1)];
end

% Create box plot
figure;
set(gcf, 'Units', 'pixels', 'Position', [1250, 200, 500, 400]); % [x, y, width, height]
boxplot(data_all, group_all);
xlabel('Reconstruction Method');
ylabel('\Deltan [a.u.]' );
title(sprintf('%s Dispersion of \\Deltan per Method', upper(metric_name)));
grid minor;

% Plot using boxchart
% figure;
% boxchart(categorical(group_all, methods_unique), data_all);
% xlabel('Reconstruction Method');
% ylabel(sprintf('%s of \\Deltan map', upper(metric_name)));
% title(sprintf('%s Dispersion of \\Deltan per Method', upper(metric_name)));
% grid on;

