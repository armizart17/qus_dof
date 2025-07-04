% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
% setup,
% warning('off'); %% Important to turn warns in loop off

% baseDir = 'LiverAcquisition_24_09_07';
baseDir = pwd;
resultsDir = fullfile(baseDir,'results');

%% Sebastian Merino

resultFiles = dir(fullfile(resultsDir,'*001-dc*_results.mat'));
nFiles = length(resultFiles);

varNames = {'acs','method','patient','sample','dc'};
varTypes = {'double','string','int16','int16','logical'};
T = table('Size',[nFiles*3,5], 'VariableTypes',varTypes,'VariableNames',varNames);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample,true};
    T(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample,true};
    T(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample,true};
end


resultFiles = dir(fullfile(resultsDir,'001-0*_results.mat'));
nFiles = length(resultFiles);

varNames = {'acs','method','patient','sample','dc'};
varTypes = {'double','string','int16','int16','logical'};
T2 = table('Size',[nFiles*3,5], 'VariableTypes',varTypes,'VariableNames',varNames);

resultFiles = dir(fullfile(resultsDir,'*_results.mat'));
nFiles = length(resultFiles);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T2(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample,false};
    T2(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample,false};
    T2(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample,false};
end
T = [T;T2];

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(T.method, {'RSLD','SWTV','SWIFT'}), ...
    T.acs, 'MarkerStyle','o', 'GroupByColor',T.dc);
ylim([-0.4,1.6])
grid on
title('ACS, SM, after eating')
legend


%% JT
resultFiles = dir(fullfile(resultsDir,'*017-DC*_results.mat'));
nFiles = length(resultFiles);

varNames = {'acs','method','patient','sample','dc'};
varTypes = {'double','string','int16','int16','logical'};
T = table('Size',[nFiles*3,5], 'VariableTypes',varTypes,'VariableNames',varNames);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample,true};
    T(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample,true};
    T(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample,true};
end


resultFiles = dir(fullfile(resultsDir,'017-0*_results.mat'));
nFiles = length(resultFiles);

varNames = {'acs','method','patient','sample','dc'};
varTypes = {'double','string','int16','int16','logical'};
T2 = table('Size',[nFiles*3,5], 'VariableTypes',varTypes,'VariableNames',varNames);

resultFiles = dir(fullfile(resultsDir,'*_results.mat'));
nFiles = length(resultFiles);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T2(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample,false};
    T2(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample,false};
    T2(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample,false};
end
T = [T;T2];

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(T.method, {'RSLD','SWTV','SWIFT'}), ...
    T.acs, 'MarkerStyle','o', 'GroupByColor',T.dc);
ylim([-0.4,1.6])
grid on
title('ACS, JT, after eating')
legend
