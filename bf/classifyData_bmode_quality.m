%% INSPECT DATA
% For .mat files (e.g., SAM.bMode inside)
baseDir     = 'D:\emirandaz\qus\data\liver\patients_IdepFoci\';
sourceDir   = fullfile(baseDir, 'Beamformed_ClinicalData_Patients_IndepFoci');
tableFile   = 'labelsEMZ.xlsx';

T = label_bmode_quality(fullfile(baseDir,sourceDir), tableFile, '.mat', '.png');

%% COPY FILES IN ANOTHER FOLDER

% Parameters
baseDir     = 'D:\emirandaz\qus\data\liver\patients_IdepFoci\';
sourceDir   = fullfile(baseDir, 'Beamformed_ClinicalData_Patients_IndepFoci');
targetDir   = fullfile(baseDir, 'Beamformed_ClinicalData_Patients_IndepFoci_filt1');
labelsFile  = 'labelsEMZ.xlsx';  % File with label table

% Load label table and filter
T = readtable(labelsFile);
T_selected = T(T.MEDIUM == 1 | T.HIGH == 1, :);

% File types you want to copy
fileTypes = {'.mat', '.png'};  % Add more if needed

% Create destination folder if it doesn't exist
if ~exist(targetDir, 'dir')
    mkdir(targetDir);
end

% Loop through all selected files and copy both .mat and .png
for i = 1:height(T_selected)
    name = T_selected.nameFile{i};

    for j = 1:length(fileTypes)
        ext = fileTypes{j};
        src = fullfile(sourceDir, [name ext]);
        dst = fullfile(targetDir, [name ext]);

        if isfile(src)
            copyfile(src, dst);
        else
            warning('Missing: %s', src);
        end
    end
end
fprintf('Copied %d sets of files (.mat + .png) to:\n %s\n', height(T_selected), targetDir);