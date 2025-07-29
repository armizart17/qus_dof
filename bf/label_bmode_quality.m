function T = label_bmode_quality(rootDir, labelsFile, rfExt, imgExt)
% function T = label_bmode_quality(rootDir, labelsFile, rfExt, imgExt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LABEL_BMODE_QUALITY  Annotate B-mode quality (LOW, MEDIUM, HIGH).
% Reads images ".png" or ".mat" files and shows B-mode with interface.
% Usage:
%   label_bmode_quality(rootDir, labelsFile, rfExt, imgExt)
%   Example:
%     label_bmode_quality('C:\data\', 'labels.xlsx', '.mat', '.png')
%
% Features:
% - Interactive labeling (HIGH/MEDIUM/LOW/EXIT)
% - Excel or CSV file support
% - Examiner name and timestamp saved
% - Supports .mat (with SAM.bMode) and .png files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % === Inputs ===
    if nargin < 1 || isempty(rootDir)
        rootDir = uigetdir(pwd, 'Select folder with RF or image files');
        if isequal(rootDir, 0), return; end
    end
    if nargin < 2 || isempty(labelsFile)
        [file, path] = uiputfile({'*.xlsx;*.csv','Excel or CSV (*.xlsx, *.csv)'}, ...
                                 'Choose file to save labels');
        if isequal(file,0), return; end
        labelsFile = fullfile(path, file);
    end
    if nargin < 3 || isempty(rfExt),  rfExt = ".mat"; end
    if nargin < 4 || isempty(imgExt), imgExt = ".png"; end

    [~,~,ext] = fileparts(labelsFile);
    if ~ismember(lower(ext), {'.xlsx', '.csv'})
        error('Label file must end in .xlsx or .csv');
    end

    % Ask examiner initials
    examiner = input('Enter your initials (e.g., EMZ): ', 's');

    % Load or initialize label table
    if isfile(labelsFile)
        T = readtable(labelsFile, 'TextType','string');
    else
        nFiles = numel(dir(fullfile(rootDir, ['*' rfExt])));
        T = table(strings(nFiles,1), zeros(nFiles,1), zeros(nFiles,1), zeros(nFiles,1), ...
                  repmat(string(examiner), nFiles,1), NaT(nFiles,1), ...
                  'VariableNames', {'nameFile','LOW','MEDIUM','HIGH','examiner','timestamp'});
        T(1,:) = [];
    end

    % Find candidate files to label
    allFiles = dir(fullfile(rootDir, ['*' rfExt]));
    if isempty(allFiles)
        error('No %s files found in %s', rfExt, rootDir);
    end

    for k = 1:numel(allFiles)
        fileName = allFiles(k).name;
        [~, baseName] = fileparts(fileName);

        if any(T.nameFile == baseName)
            fprintf('Skipping %s (already labeled)\n', baseName);
            continue;
        end

        %% === Display B-mode ===
        fig = figure('Name', baseName, ...
                     'Units','normalized','Position',[0.1 0.1 0.7 0.6]);

        caption = strrep(baseName, '_', ' ');
        fontSize = 12;

        if endsWith(fileName, '.mat')
            SAM = load(fullfile(rootDir, fileName));
            if ~isfield(SAM, 'bMode'), SAM.bmode = my_RF2Bmode(SAM.rf(:,:,1)); end

            if (SAM.xr(end) < 1 && SAM.xr(1)> -1)
                SAM.xr = SAM.xr*1e3;
            end

            t = tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

            nexttile
            imagesc(SAM.xr, SAM.zr*1e3, SAM.bMode);
            colorbar; clim([-55 0]); colormap gray;
            ylabel('Depth [mm]', 'FontSize', 14);
            xlabel('\theta [°]', 'FontSize', 14);
            axis image
            title(caption, 'FontSize', fontSize);

            nexttile
            pcolor(SAM.xp*1e3, SAM.zp*1e3, SAM.bMode); shading interp;
            colorbar; clim([-55 0]); colormap gray;
            ylabel('Depth [mm]', 'FontSize', 14);
            xlabel('Lateral [mm]', 'FontSize', 14);
            axis equal ij tight;
            title(caption, 'FontSize', fontSize);
            set(gca, 'Color', 'k');

        elseif endsWith(fileName, imgExt)
            I = imread(fullfile(rootDir, fileName));
            imshow(I, []);
            title(caption, 'Interpreter','none', 'FontSize', fontSize);
        else
            warning('Unsupported file: %s', fileName);
            close(fig);
            continue;
        end

        drawnow;

        %% === Ask for quality ===
        confirmed = false;
        while ~confirmed
            choice = questdlg(sprintf('Quality of "%s"?', baseName), ...
                              'Select Quality', ...
                              'HIGH','MEDIUM','LOW','EXIT');

            if isempty(choice) || strcmp(choice, 'EXIT')
                fprintf('\nLabeling exited by user. Saving progress...\n');
                writetable(T, labelsFile);
                close(fig);
                return;
            end

            confirm = questdlg(sprintf('Are you sure you want to label "%s" as "%s"?', baseName, choice), ...
                               'Confirm', 'Yes', 'No', 'Yes');
            if strcmp(confirm, 'Yes')
                confirmed = true;
            end
        end

        close(fig);

        %% === Store label ===
        newRow = { baseName, ...
                   double(strcmp(choice,'LOW')), ...
                   double(strcmp(choice,'MEDIUM')), ...
                   double(strcmp(choice,'HIGH')), ...
                   string(examiner), datetime('now') };
        T = [T; newRow]; %#ok<AGROW>

        if mod(height(T),10) == 0
            writetable(T, labelsFile);
            fprintf('Progress saved to %s (%d entries).\n', labelsFile, height(T));
        end
    end

    %% Final save
    writetable(T, labelsFile);
    fprintf('\nDone! Labels saved to:\n  %s\n', labelsFile);
end
