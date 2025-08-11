% ====================================================================== %
% Post processing of Liver Healthy volunteers 
% Data (Acq in 2024 in LIM, 20 "healthy" volunteers), reference 544
% Processed by DoF method
% Post-processing all cases dispersion boxcharts, violins, etc.
% ====================================================================== %
%%
init
warning('off'); %% Important to turn warns in loop off
%% DIRECTORIES

baseDir     = 'D:\emirandaz\qus\data\liver\healthy';
% resultsDir  = fullfile(baseDir,'results_LSv2'); % results by LS method no regu
resultsDir  = fullfile(baseDir,'results_RPL_TV'); % results by LS method no regu

% Result files (Processed by DoF method)
resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles      = length(resultFiles);
%% METHOD LABELS
%% METHOD LABELS

methods      = {'3-DoF', '2-DoF-a', '2-DoF-b', '2-DoF-n'};
method_labels = { ...
    '\mathrm{3\textrm{-}DoF}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,n}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{n,a}}', ...
    '\mathrm{2\textrm{-}DoF}_{\mathrm{b,a}}' ...
};

nMethods = length(methods);

varNames = {'idName','num',  'method','mean',  'std',   'median','iqr'};
varTypes = {'string','int16','string','double','double','double','double'};
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
    idName      = fileName(1:end-4);
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_alpha(jj);
        img_map = maps_results_dof{iMethod}.alpha;

        Ta(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    };
        row = row + 1;
    end
end
%
% writetable(Ta,fullfile(resultsDir,'meanACS.xlsx'))

% ================== Delta b metrics ==================

indices_b = [1, 2, 4];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_b(jj);
        img_map = maps_results_dof{iMethod}.b_dB;

        Tb(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    };
        row = row + 1;
    end
end

% ================== Delta n metrics ==================

indices_n = [1, 2, 3];  % 3-DoF, 2-Dof-a, 2-DoF-b, 2-DoF-n

row = 1;
for ii = 0:nFiles-1
    fileName    = resultFiles(ii+1).name;
    idName      = fileName(1:end-4);
    load(fullfile(resultsDir,fileName));
    load(fullfile(resultsDir,fileName));
    
    for jj = 1:nMethods-1 % always 3 methods 
        iMethod = indices_n(jj);
        img_map = maps_results_dof{iMethod}.n;

        Tn(row,:) = {idName, row, methods{iMethod}, ...
                    mean(img_map(:), 'omitnan'), ...
                    std(img_map(:), 'omitnan'), ...
                    median(img_map(:), 'omitnan'), ...
                    iqr(img_map(:)), ...
                    };
        row = row + 1;
    end
end


%% BOX PLOTS
% =============== Box alpha ===============

metric_name = 'std';  % Change this to 'mean', 'std', or 'iqr' as needed

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


%%
figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(T.method, {'RSLD','SWTV','SWIFT'}), ...
    T.coeff, 'MarkerStyle','o');
ylim([0.2,1])
grid on
title('ACS')
% legend



%%
Tswift = T(T.method=='SWIFT',:);

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(Tswift.patient), ...
    Tswift.coeff, 'MarkerStyle','o');
ylim([0.2,1])
grid on
title('ACS')
% yline(0.64, 'k--')
xlabel('Patient')
ylabel('ACS [dB/cm/MHz]')

disp(mean(Tswift.coeff))
disp(std(Tswift.coeff))

%%
Tswtv = T(T.method=='SWTV',:);

% figure('Units','centimeters', 'Position', [5 5 20 10]),
% boxchart(categorical(Tswtv.patient), ...
%     Tswtv.coeff, 'MarkerStyle','o');
% ylim([0.4,1])
% grid on
% title('Median ACS with SWTV')
% % yline(0.64, 'k--')
% xlabel('Patient')
% ylabel('ACS [dB/cm/MHz]')


disp(mean(Tswtv.coeff))
disp(std(Tswtv.coeff))
%%
Trsld = T(T.method=='RSLD',:);

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(Trsld.patient), ...
    Trsld.coeff, 'MarkerStyle','o');
ylim([0.4,1])
grid on
title('Median ACS with RSLD')
% yline(0.64, 'k--')
xlabel('Patient')
ylabel('ACS [dB/cm/MHz]')


%%
Tbsc = T(T.method=='BSC',:);

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(Tbsc.patient), ...
    Tbsc.coeff, 'MarkerStyle','o');

set(gca,'YScale','log');
ylim([1e-5 1e-2])
grid on
title('BSC')
% yline(0.64, 'k--')
xlabel('Patient')
ylabel('BSC [1/cm-Sr]')

disp(mean(Tbsc.coeff))
disp(std(Tbsc.coeff))

%% Extra info
Tinfo = readtable('volunteers.xlsx');

%%


G = findgroups(Tswtv.patient);
acsPatientRsld = splitapply(@median,Trsld.coeff,G);
acsPatientSwtv = splitapply(@median,Tswtv.coeff,G);
acsPatientSwift = splitapply(@median,Tswift.coeff,G);
bscPatient = splitapply(@median,Tbsc.coeff,G);



iqrPatientRsld = splitapply(@iqr,Trsld.coeff,G);
iqrPatientSwtv = splitapply(@iqr,Tswtv.coeff,G);
iqrPatientSwift = splitapply(@iqr,Tswift.coeff,G);
iqrPatientbsc = splitapply(@iqr, Tbsc.coeff,G);


idsPatient = splitapply(@(x) x(1),Tswtv.patient,G);
sexPatient = Tinfo.sex(idsPatient);


figure('Units','centimeters', 'Position', [5 5 18 10]),
tiledlayout(1,3)
nexttile,
boxchart(categorical(sexPatient), ...
    acsPatientSwtv, 'MarkerStyle','o');
ylim([0.4,1])
grid on
title('ACS per patient, SWTV')
% yline(0.64, 'k--')
xlabel('Sex')
ylabel('ACS [dB/cm/MHz]')

nexttile,
boxchart(categorical(sexPatient), ...
    acsPatientSwift, 'MarkerStyle','o');
ylim([0.4,1])
grid on
title('ACS per patient')
% yline(0.64, 'k--')
xlabel('Sex')
ylabel('ACS [dB/cm/MHz]')

nexttile,
boxchart(categorical(sexPatient), ...
    bscPatient, 'MarkerStyle','o');
set(gca,'YScale','log');
grid on
title('BSC per patient')
% yline(0.64, 'k--')
xlabel('Sex')
ylabel('BSC [1/cm-Sr]')
ylim([1e-5 1e-2])



%%
% acsSex = acsPatientSwtv(sexPatient=="F");
% fprintf("%.2f %.2f\n\n", median(acsSex),iqr(acsSex))

% acsSex = acsPatientSwtv(sexPatient=="M");
% fprintf("%.2f %.2f\n\n", median(acsSex),iqr(acsSex))
%%
acsSex = acsPatientSwift(sexPatient=="F");
fprintf("ACS %.2f %.2f\n\n", median(acsSex),iqr(acsSex))

acsSex = acsPatientSwift(sexPatient=="M");
fprintf("ACS %.2f %.2f\n\n", median(acsSex),iqr(acsSex))


%%
bscSex = bscPatient(sexPatient=="F");
fprintf("BSC %.2f %.2f\n\n", median(bscSex)*1e4,iqr(bscSex)*1e4)

bscSex = bscPatient(sexPatient=="M");
fprintf("BSC %.2f %.2f\n\n", median(bscSex)*1e4,iqr(bscSex)*1e4)
%%
bmiPatient = Tinfo.bmi(idsPatient);

figure('Units','centimeters', 'Position', [5 5 20 10]),
tiledlayout(1,3)
nexttile,
scatter(bmiPatient,acsPatientSwtv)
grid on
title('Median ACS, SWTV')
xlabel('BMI')
ylabel('ACS [dB/cm/MHz]')
xlim([20 35])
ylim([0.5 0.9])

nexttile,
scatter(bmiPatient,acsPatientSwift)
grid on
title('Median ACS')
xlabel('BMI')
ylabel('ACS [dB/cm/MHz]')
xlim([20 35])
ylim([0.5 0.9])

nexttile,
scatter(bmiPatient,bscPatient)
grid on
title('Median BSC')
xlabel('BMI')
ylabel('BSC [1/cm-Sr]')
ylim([1e-5 1e-2])
set(gca,'YScale','log');

%%
% finalTable = array2table([acsPatientRsld, stdPatientRsld, ...
%     acsPatientSwtv, stdPatientSwtv, acsPatientSwift, stdPatientSwift], ...
%     'VariableNames', {'rsldMean','rsldStd','swtvMean','swtvStd', ...
%     'swiftMean','swiftStd'});
finalTable = array2table([acsPatientRsld, iqrPatientRsld, ...
    acsPatientSwtv, iqrPatientSwtv, acsPatientSwift, iqrPatientSwift, ...
    bscPatient*1e4, iqrPatientbsc*1e4], ...
    'VariableNames', {'rsldMedian','rsldIqr','swtvMedian','swtvIqr', ...
    'swiftMedian','swiftIqr', 'bscMedian', 'bscIqr'});
writetable(finalTable,'allResults.xlsx')
