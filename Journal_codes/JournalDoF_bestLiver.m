%%
% init
warning('off'); %% Important to turn warns in loop off

font_size = 30; 
%% BEST CASES

% healthyNames = {'smintercostal3_F', ...
%     'arintercostal2_F', ...
%     % 'arsubcostal_F', ...  % commented lines are fine
%     % 'arsubcostal2_F', ...
%     'ASIntercostal_F_2' ...
%     };

healthyNames = { ...
    'smintercostal3_F', ...
    'arintercostal_F_2', ...
    '73089254_IOLAI_F', ...
    % '73089254_LHI_F'
    % 'ASIntercostal_F_2' ...
        % 'arintercostal2_F', ...

};

steatoticNames = { ...
            
            % '70897309_IOLAI_F2', ...
            % '80296823_IOLAI_F1', ...
            % '80296823_IOLAI_F2', ...
            % '000345400_IOLAI_F1', ... % double check
            '000345400_IHR_F1', ... % maybe
            '70141854_PL_F2', ...
            '41681509_IHR_F1', ...
            % '46747475_IOLAI_F1', ...
            % '41693885_IHR_F1', ...
            };  
% 
% folderSteatotic = 'RPL_mu1_lu_patient_out';
% folderControl   = 'RPL_mu1v2_lu_control_out';

% folderSteatotic = 'RPL_mu10_lu_patient_out';
% folderControl   = 'RPL_mu1v3_lu_control_out';

folderControl   = 'RPL_mu50_lu_control_out'
folderSteatotic = 'RPL_mu50_lu_patient_out';
%folderSteatotic = 'RPL_mu10_lu_patient_out';

% folderControl   = 'RPL_mu10_lu_control_out';
% folderSteatotic = 'RPL_mu10_lu_patient_out';


% folderControl   = 'RPL_lu_control_out';
% folderSteatotic = 'RPL_lu_Patientfilt1_out';


%% DIRECTORIES HEALTHY CONTROL

baseDir     = 'D:\emirandaz\qus\data\liver\bf_SavedDataCurvilinearPW\';
resultsDir  = fullfile(baseDir, folderControl); 


metric_name = 'std';  % Change this to 'mean', 'std', 'median' or 'iqr' as needed

% Result files (Processed by DoF method)
resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles      = length(resultFiles);

% METHOD LABELS

methods      = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

nMethods = length(method_labels);

varNames = {'idName','num',  'method','mean',  'std',   'median','iqr','map'};
varTypes = {'string','int16','string','double','double','double','double', 'cell'};
T       = table('Size',[nFiles*(nMethods-1), length(varNames)], ...
            'VariableTypes',varTypes,'VariableNames',varNames);

Ta = T;
Tb = T;
Tn = T;
clear T

% ================== alpha metrics ==================

indices_alpha = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, healthyNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_alpha(jj);
        img_map = maps_results_dof{iMethod}.alpha;

        Ta(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map}, ...
                    };
        row = row + 1;
    end
end
%
% writetable(Ta,fullfile(resultsDir,'meanACS.xlsx'))

% ================== Delta b metrics ==================

indices_b = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, healthyNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_b(jj);
        img_map = maps_results_dof{iMethod}.b_dB;

        Tb(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map},...
                    };
        row = row + 1;
    end
end

% ================== Delta n metrics ==================

indices_n = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, healthyNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_n(jj);
        img_map = maps_results_dof{iMethod}.n;

        Tn(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map},...
                    };
        row = row + 1;
    end
end

%% Healthy (ISOLATE METHODS)

Ha_3 = Ta(Ta.method == '3-DoF', :);
Ha_n = Ta(Ta.method == '2-DoF-n', :);

Hb_3 = Tb(Tb.method == '3-DoF', :);
Hb_n = Tb(Tb.method == '2-DoF-n', :);

Hn_3 = Tn(Tn.method == '3-DoF', :);
Hn_n = Tn(Tn.method == '2-DoF-n', :);

%% PATIENTS STEATOTIC

baseDir     = 'D:\emirandaz\qus\data\liver\patients_IdepFoci';
resultsDir  = fullfile(baseDir, folderSteatotic); 


% Result files (Processed by DoF method)
resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles      = length(resultFiles);

methods      = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};

method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

nMethods = length(method_labels);

varNames = {'idName','num',  'method','mean',  'std',   'median','iqr','map'};
varTypes = {'string','int16','string','double','double','double','double', 'cell'};
T       = table('Size',[nFiles*(nMethods-1), length(varNames)], ...
            'VariableTypes',varTypes,'VariableNames',varNames);

Ta = T;
Tb = T;
Tn = T;
clear T

% ================== alpha metrics ==================

indices_alpha = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, steatoticNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_alpha(jj);
        img_map = maps_results_dof{iMethod}.alpha;

        Ta(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map}, ...
                    };
        row = row + 1;
    end
end
%
% writetable(Ta,fullfile(resultsDir,'meanACS.xlsx'))

% ================== Delta b metrics ==================

indices_b = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, steatoticNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_b(jj);
        img_map = maps_results_dof{iMethod}.b_dB;

        Tb(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map},...
                    };
        row = row + 1;
    end
end

% ================== Delta n metrics ==================

indices_n = [1, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    % Only continue if this file is in healthyNames
    if ~ismember(idName, steatoticNames)
        continue
    end
    
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods % always 2
        iMethod = indices_n(jj);
        img_map = maps_results_dof{iMethod}.n;

        Tn(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    {img_map},...
                    };
        row = row + 1;
    end
end

%% STEATOTIC (ISOLATE METHODS)

Sa_3 = Ta(Ta.method == '3-DoF', :);
Sa_n = Ta(Ta.method == '2-DoF-n', :);

Sb_3 = Tb(Tb.method == '3-DoF', :);
Sb_n = Tb(Tb.method == '2-DoF-n', :);

Sn_3 = Tn(Tn.method == '3-DoF', :);
Sn_n = Tn(Tn.method == '2-DoF-n', :);

%% 
% keyboard


%% Select rows matching wantedNames

Ha_n = Ha_n(ismember(Ha_n.idName, healthyNames), :);
Ha_3 = Ha_3(ismember(Ha_3.idName, healthyNames), :);

Hb_n = Hb_n(ismember(Hb_n.idName, healthyNames), :);
Hb_3 = Hb_3(ismember(Hb_3.idName, healthyNames), :);

Hn_n = Hn_n(ismember(Hn_n.idName, healthyNames), :);
Hn_3 = Hn_3(ismember(Hn_3.idName, healthyNames), :);

Sa_n = Sa_n(ismember(Sa_n.idName, steatoticNames), :);
Sa_3 = Sa_3(ismember(Sa_3.idName, steatoticNames), :);

Sb_n = Sb_n(ismember(Sb_n.idName, steatoticNames), :);
Sb_3 = Sb_3(ismember(Sb_3.idName, steatoticNames), :);

Sn_n = Sn_n(ismember(Sn_n.idName, steatoticNames), :);
Sn_3 = Sn_3(ismember(Sn_3.idName, steatoticNames), :);


Ha_dof = [Ha_3; Ha_n];
Hb_dof = [Hb_3; Hb_n];
Hn_dof = [Hn_3; Hn_n];

Sa_dof = [Sa_3; Sa_n];
Sb_dof = [Sb_3; Sb_n];
Sn_dof = [Sn_3; Sn_n];

%% ALPHA from maps

% =====================
% 1. Prepare data from tables (extract all pixel values)
% =====================
% Healthy cases: 3-DoF and 2-DoF_{b,a}
data_3dof_healthy = cellfun(@(m) m(:), Ha_3.map, 'UniformOutput', false);
data_2dof_healthy = cellfun(@(m) m(:), Ha_n.map, 'UniformOutput', false);

% Steatotic cases: 3-DoF and 2-DoF_{b,a}
data_3dof_steato = cellfun(@(m) m(:), Sa_3.map, 'UniformOutput', false);
data_2dof_steato = cellfun(@(m) m(:), Sa_n.map, 'UniformOutput', false);

% =====================
% 2. Combine all data
% =====================
% =====================
% 2. Interleave the pairs
% =====================
all_data = {};
for i = 1:3
    % Healthy
    all_data{end+1} = data_3dof_healthy{i};
    all_data{end+1} = data_2dof_healthy{i};
end
for i = 1:3
    % Steatotic
    all_data{end+1} = data_3dof_steato{i};
    all_data{end+1} = data_2dof_steato{i};
end

data_vector = vertcat(all_data{:});
group_number = repelem(1:12, cellfun(@numel, all_data));

% =====================
% 3. Labels and positions
% =====================
xtick_labels = {'H1','H1','H2','H2','H3','H3', ...
                'S1','S1','S2','S2','S3','S3'};
xtick_pos = 1:12;

% Colors
bright_blue   = [0 0.4470 0.7410];   % 3-DoF healthy
bright_orange = [0.8500 0.3250 0.0980]; % 2-DoF healthy
dark_blue     = bright_blue * 1;     % darker for steatotic
dark_orange   = bright_orange * 1;

% =====================
% 4. Create boxplot
% =====================
figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 800]);

boxplot(data_vector, group_number, 'Widths', 0.6, 'Colors', 'k');
set(gca, 'xtick', xtick_pos, 'xticklabel', xtick_labels);
xlim([0.5, 12.5]);
ylim([-0.5 2.5])
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
title('Estimation of \alpha for Human Liver Cases');
set(gca, 'FontSize', font_size);
grid on;

% =====================
% 5. Color the boxes
% =====================
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    idx = length(h) - j + 1; % reverse order for boxplot handles
    if idx <= 6
        % Healthy
        c = bright_blue;
        if mod(idx, 2) == 0
            c = bright_orange;
        end
    else
        % Steatotic
        c = dark_blue;
        if mod(idx, 2) == 0
            c = dark_orange;
        end
    end
    patch(get(h(j), 'XData'), get(h(j), 'YData'), c, 'FaceAlpha', 0.5);
end

% =====================
% 6. Divider between healthy and steatotic
% =====================
hold on;
xline(6.5, '--k', 'LineWidth', 1.2);

% =====================
% 7. Legend
% =====================
h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_blue, ...
    'MarkerEdgeColor', bright_blue, 'DisplayName', '3-DoF');
h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_orange, ...
    'MarkerEdgeColor', bright_orange, 'DisplayName', '2-DoF_{b,a}');
legend([h1, h2], 'Location', 'southwest');

%% \Delta b

% =====================
% 1. Prepare data from tables (extract all pixel values)
% =====================
% Healthy cases: 3-DoF and 2-DoF_{b,a}
data_3dof_healthy = cellfun(@(m) m(:), Hb_3.map, 'UniformOutput', false);
data_2dof_healthy = cellfun(@(m) m(:), Hb_n.map, 'UniformOutput', false);

% Steatotic cases: 3-DoF and 2-DoF_{b,a}
data_3dof_steato = cellfun(@(m) m(:), Sb_3.map, 'UniformOutput', false);
data_2dof_steato = cellfun(@(m) m(:), Sb_n.map, 'UniformOutput', false);

% =====================
% 2. Combine all data
% =====================
% =====================
% 2. Interleave the pairs
% =====================
all_data = {};
for i = 1:3
    % Healthy
    all_data{end+1} = data_3dof_healthy{i};
    all_data{end+1} = data_2dof_healthy{i};
end
for i = 1:3
    % Steatotic
    all_data{end+1} = data_3dof_steato{i};
    all_data{end+1} = data_2dof_steato{i};
end

data_vector = vertcat(all_data{:});
group_number = repelem(1:12, cellfun(@numel, all_data));

% =====================
% 3. Labels and positions
% =====================
xtick_labels = {'H1','H1','H2','H2','H3','H3', ...
                'S1','S1','S2','S2','S3','S3'};
xtick_pos = 1:12;

% Colors
bright_blue   = [0 0.4470 0.7410];   % 3-DoF healthy
bright_orange = [0.8500 0.3250 0.0980]; % 2-DoF healthy
dark_blue     = bright_blue * 1;     % darker for steatotic
dark_orange   = bright_orange * 1;

% =====================
% 4. Create boxplot
% =====================

figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 800]);

boxplot(data_vector, group_number, 'Widths', 0.6, 'Colors', 'k');
set(gca, 'xtick', xtick_pos, 'xticklabel', xtick_labels);
xlim([0.5, 12.5]);
% ylim([-0.5 2.5])
ylabel('\Deltab [dB]');
title('Estimation of \Deltab for Human Liver Cases');
set(gca, 'FontSize', font_size);
grid on;

% =====================
% 5. Color the boxes
% =====================
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    idx = length(h) - j + 1; % reverse order for boxplot handles
    if idx <= 6
        % Healthy
        c = bright_blue;
        if mod(idx, 2) == 0
            c = bright_orange;
        end
    else
        % Steatotic
        c = dark_blue;
        if mod(idx, 2) == 0
            c = dark_orange;
        end
    end
    patch(get(h(j), 'XData'), get(h(j), 'YData'), c, 'FaceAlpha', 0.5);
end

% =====================
% 6. Divider between healthy and steatotic
% =====================
hold on;
xline(6.5, '--k', 'LineWidth', 1.2);

% =====================
% 7. Legend
% =====================
h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_blue, ...
    'MarkerEdgeColor', bright_blue, 'DisplayName', '3-DoF');
h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_orange, ...
    'MarkerEdgeColor', bright_orange, 'DisplayName', '2-DoF_{b,a}');
legend([h1, h2], 'Location', 'southwest');


%% \Delta n

% =====================
% 1. Prepare data from tables (extract all pixel values)
% =====================
% Healthy cases: 3-DoF and 2-DoF_{b,a}
data_3dof_healthy = cellfun(@(m) m(:), Hn_3.map, 'UniformOutput', false);
data_2dof_healthy = cellfun(@(m) m(:), Hn_n.map, 'UniformOutput', false);

% Steatotic cases: 3-DoF and 2-DoF_{b,a}
data_3dof_steato = cellfun(@(m) m(:), Sn_3.map, 'UniformOutput', false);
data_2dof_steato = cellfun(@(m) m(:), Sn_n.map, 'UniformOutput', false);

% =====================
% 2. Combine all data
% =====================
% =====================
% 2. Interleave the pairs
% =====================
all_data = {};
for i = 1:3
    % Healthy
    all_data{end+1} = data_3dof_healthy{i};
    all_data{end+1} = data_2dof_healthy{i};
end
for i = 1:3
    % Steatotic
    all_data{end+1} = data_3dof_steato{i};
    all_data{end+1} = data_2dof_steato{i};
end

data_vector = vertcat(all_data{:});
group_number = repelem(1:12, cellfun(@numel, all_data));

% =====================
% 3. Labels and positions
% =====================
xtick_labels = {'H1','H1','H2','H2','H3','H3', ...
                'S1','S1','S2','S2','S3','S3'};
xtick_pos = 1:12;

% Colors
bright_blue   = [0 0.4470 0.7410];   % 3-DoF healthy
bright_orange = [0.8500 0.3250 0.0980]; % 2-DoF healthy
dark_blue     = bright_blue * 1;     % darker for steatotic
dark_orange   = bright_orange * 1;

% =====================
% 4. Create boxplot
% =====================

figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 1600, 800]);

boxplot(data_vector, group_number, 'Widths', 0.6, 'Colors', 'k');
set(gca, 'xtick', xtick_pos, 'xticklabel', xtick_labels);
xlim([0.5, 12.5]);
% ylim([-0.5 2.5])
ylabel('\Deltan [a.u.]');
title('Estimation of \Deltan for Human Liver Cases');
set(gca, 'FontSize', font_size);
grid on;

% =====================
% 5. Color the boxes
% =====================
h = findobj(gca, 'Tag', 'Box');
for j = 1:length(h)
    idx = length(h) - j + 1; % reverse order for boxplot handles
    if idx <= 6
        % Healthy
        c = bright_blue;
        if mod(idx, 2) == 0
            c = bright_orange;
        end
    else
        % Steatotic
        c = dark_blue;
        if mod(idx, 2) == 0
            c = dark_orange;
        end
    end
    patch(get(h(j), 'XData'), get(h(j), 'YData'), c, 'FaceAlpha', 0.5);
end

% =====================
% 6. Divider between healthy and steatotic
% =====================
hold on;
xline(6.5, '--k', 'LineWidth', 1.2);

% =====================
% 7. Legend
% =====================
h1 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_blue, ...
    'MarkerEdgeColor', bright_blue, 'DisplayName', '3-DoF');
h2 = plot(NaN, NaN, 's', 'MarkerFaceColor', bright_orange, ...
    'MarkerEdgeColor', bright_orange, 'DisplayName', '2-DoF_{b,a}');
legend([h1, h2], 'Location', 'southwest');


%% ALPHA mean ± std from maps
% =====================
% 1. Prepare stats function
% =====================
getStats = @(mapsCell) cellfun(@(m) ...
    [median(m(:),'omitnan'), std(m(:),'omitnan')], ...
    mapsCell, 'UniformOutput', false);

% Healthy
tmp = getStats(Ha_3.map);
stats_3dof_healthy = vertcat(tmp{:});

tmp = getStats(Ha_n.map);
stats_2dof_healthy = vertcat(tmp{:});

% Steatotic
tmp = getStats(Sa_3.map);
stats_3dof_steato = vertcat(tmp{:});

tmp = getStats(Sa_n.map);
stats_2dof_steato = vertcat(tmp{:});

% =====================
% 2. Combine for plotting
% =====================
means_all = [stats_3dof_healthy(:,1)'; stats_2dof_healthy(:,1)'; ...
             stats_3dof_steato(:,1)'; stats_2dof_steato(:,1)'];

stds_all  = [stats_3dof_healthy(:,2)'; stats_2dof_healthy(:,2)'; ...
             stats_3dof_steato(:,2)'; stats_2dof_steato(:,2)'];

% Flatten to vectors for plotting (interleave)
mean_vals = means_all(:);
std_vals  = stds_all(:);

% X positions for groups
x_pos = 1:length(mean_vals);

% =====================
% 3. Colors
% =====================
bright_blue   = [0 0.4470 0.7410];    % 3-DoF healthy
bright_orange = [0.8500 0.3250 0.0980]; % 2-DoF healthy
dark_blue     = bright_blue;      % steatotic
dark_orange   = bright_orange;

colors = [bright_blue; bright_orange; ...
          bright_blue; bright_orange; ...
          bright_blue; bright_orange; ...
          dark_blue; dark_orange; ...
          dark_blue; dark_orange; ...
          dark_blue; dark_orange];

% =====================
% 4. Plot
% =====================
figure;
hold on;
for i = 1:length(mean_vals)
    errorbar(x_pos(i), mean_vals(i), std_vals(i), 's', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'Color', colors(i,:), 'LineWidth', 1.5, 'CapSize', 10);
end

% Formatting
xticks(1:12);
xticklabels({'H1','H1','H2','H2','H3','H3','S1','S1','S2','S2','S3','S3'});
ylabel('\alpha [dB\cdotcm^{-1}\cdotMHz^{-1}]');
title('Mean ± STD of \alpha for Human Liver Cases');
xline(6.5, '--k', 'LineWidth', 1.2); % divider
grid on;

% Legend
h1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', bright_blue, ...
    'MarkerEdgeColor', bright_blue, 'DisplayName', '3-DoF');
h2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', bright_orange, ...
    'MarkerEdgeColor', bright_orange, 'DisplayName', '2-DoF_{b,a}');
legend([h1,h2], 'Location', 'southwest');

set(gca, 'FontSize', 14);
%%

b_wisonsin  = -55; % dB
n_wisconsin = 3.8;

%% --- Params & colors ---
freq = linspace(1.5, 3.5, 30).';       % MHz, Nx1
b_wisconsin  = -55;                  % dB
n_wisconsin  = 3.8;                  % unitless
useMedian = true;                    % true = median curve, false = mean curve
spread    = 'std';                   % 'iqr' or 'std'

bright_blue   = [0 0.4470 0.7410];
bright_orange = [0.8500 0.3250 0.0980];
dark_blue     = [0 0.20 0.40];
dark_orange   = [0.50 0.20 0.00];
line_width = 2;
font_size  = 20;

%% --- Helper: build BSC matrix (Nfreq x Nsubjects) from Δb/Δn tables for a method ---
build_bsc_mat = @(Tb, Tn, methodName) ...
    local_build_bsc(Tb, Tn, methodName, freq, b_wisconsin, n_wisconsin);

%% --- Healthy groups (expects Hb_dof, Hn_dof) ---
BSC_H_3 = build_bsc_mat(Hb_dof, Hn_dof, "3-DoF");      % N x 3
BSC_H_2 = build_bsc_mat(Hb_dof, Hn_dof, "2-DoF-n");    % N x 3   (your tables say 2-DoF-n)

%% --- Steatotic groups (expects Sb_dof, Sn_dof) ---
BSC_S_3 = build_bsc_mat(Sb_dof, Sn_dof, "3-DoF");      % N x 3
BSC_S_2 = build_bsc_mat(Sb_dof, Sn_dof, "2-DoF-n");    % N x 3

%% --- Center & spread across subjects for each group/method ---
[H3_c, H3_lo, H3_hi] = center_spread(BSC_H_3, useMedian, spread);
[H2_c, H2_lo, H2_hi] = center_spread(BSC_H_2, useMedian, spread);
[S3_c, S3_lo, S3_hi] = center_spread(BSC_S_3, useMedian, spread);
[S2_c, S2_lo, S2_hi] = center_spread(BSC_S_2, useMedian, spread);

H3_c = mean(BSC_H_3, 2);
H2_c = mean(BSC_H_2, 2);
S3_c = mean(BSC_S_3, 2);
S2_c = mean(BSC_S_2, 2);


%% --- Plot ---
font_size = 12;
figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 900, 700]); hold on; grid on;

% shaded bands
semilogy(freq, H3_c, '-',  'Color', bright_blue,   'LineWidth', line_width);
semilogy(freq, H2_c, '-',  'Color', bright_orange, 'LineWidth', line_width);
plot_band(freq, H3_lo, H3_hi, bright_blue);
plot_band(freq, H2_lo, H2_hi, bright_orange);
set(gca, 'YScale', 'log');
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title(sprintf('%s Healhty BSC (band = %s)', ternary(useMedian,'Median','Mean'), upper(spread)), 'FontSize', font_size+2);
set(gca, 'FontSize', font_size);
legend('3-DoF', '2-DoF_{b,a}','Location','best'); box on; hold off;


figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 900, 700]); hold on; grid on;
semilogy(freq, S3_c, '--', 'Color', dark_blue,     'LineWidth', line_width);
semilogy(freq, S2_c, '--', 'Color', dark_orange,   'LineWidth', line_width);
plot_band(freq, S3_lo, S3_hi, dark_blue);
plot_band(freq, S2_lo, S2_hi, dark_orange);
set(gca, 'YScale', 'log');
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
title(sprintf('%s Steatotic BSC (band = %s)', ternary(useMedian,'Median','Mean'), upper(spread)), 'FontSize', font_size+2);
set(gca, 'FontSize', font_size);
legend('3-DoF', '2-DoF_{b,a}','Location','best'); box on; hold off;

% median/mean curves

%%
figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 900, 700]); hold on; grid on;

% shaded bands
semilogy(freq, H3_c, '-',  'Color', bright_blue,   'LineWidth', line_width);
semilogy(freq, S3_c, '-',  'Color', bright_orange, 'LineWidth', line_width);
% plot_band(freq, H3_lo, H3_hi, bright_blue);
% plot_band(freq, H2_lo, H2_hi, bright_orange);
set(gca, 'YScale', 'log');
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
ylim([10^-4.5, 10^-1])
xlim([2 3.5])
title(sprintf('3-DoF %s BSC (band = %s)', ternary(useMedian,'Median','Mean'), upper(spread)), 'FontSize', font_size+2);
set(gca, 'FontSize', font_size);
legend('Healthy', 'Steatotic','Location','best'); box on; hold off;


figure;
set(gcf, 'Units', 'pixels', 'Position', [120, 120, 900, 700]); hold on; grid on;
semilogy(freq, H2_c, '-', 'Color', bright_blue,     'LineWidth', line_width);
semilogy(freq, S2_c, '-', 'Color', bright_orange,   'LineWidth', line_width);
% plot_band(freq, S3_lo, S3_hi, dark_blue);
% plot_band(freq, S2_lo, S2_hi, dark_orange);
set(gca, 'YScale', 'log');
xlabel('Frequency [MHz]', 'FontSize', font_size);
ylabel('BSC [cm^{-1}\cdot sr^{-1}]', 'FontSize', font_size);
ylim([10^-4.5, 10^-1])
xlim([2 3.5])
title(sprintf('2-DoF %s Steatotic BSC (band = %s)', ternary(useMedian,'Median','Mean'), upper(spread)), 'FontSize', font_size+2);
set(gca, 'FontSize', font_size);
legend('Healthy', 'Steatotic','Location','best'); box on; hold off;


%% ===== Helpers =====
function BSC_mat = local_build_bsc(Tb, Tn, methodName, freq, bWis_dB, nWis)
    % Filter rows by method
    Tb_m = Tb(Tb.method == methodName, :);
    Tn_m = Tn(Tn.method == methodName, :);
    % Align by idName (safety)
    [commonIDs, ia, ib] = intersect(Tb_m.idName, Tn_m.idName, 'stable');
    Tb_m = Tb_m(ia, :);
    Tn_m = Tn_m(ib, :);
    % Build per-subject BSC curves
    Nf = numel(freq);
    Ns = height(Tb_m);
    BSC_mat = nan(Nf, Ns);
    for k = 1:Ns
        % total b (dB) and n (unitless)
        b_tot_dB = Tb_m.median(k) + bWis_dB;   % Δb + prior
        n_tot    = Tn_m.median(k) + nWis;      % Δn + prior
        b_lin    = 10.^(b_tot_dB/10);          % convert dB -> linear
        BSC_mat(:,k) = b_lin .* (freq(:).^n_tot);
    end
end

function [c, lo, hi] = center_spread(M, useMedian, spread)
    % M is Nfreq x Nsubjects
    if useMedian
        c = median(M, 2, 'omitnan');
    else
        c = mean(M, 2, 'omitnan');
    end
    switch lower(spread)
        case 'iqr'
            q25 = quantile(M, 0.25, 2);
            q75 = quantile(M, 0.75, 2);
            lo = q25; hi = q75;
        case 'std'
            s  = std(M, 0, 2, 'omitnan');
            lo = c - s; hi = c + s;
        otherwise
            error('spread must be ''iqr'' or ''std''.');
    end
end

function plot_band(x, ylow, yhigh, facecol)
    fill([x(:); flipud(x(:))], [ylow(:); flipud(yhigh(:))], facecol, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

function out = ternary(cond, a, b)
if cond, out=a; else, out=b; end
end
