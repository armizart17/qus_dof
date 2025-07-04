% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
setup,
% warning('off'); %% Important to turn warns in loop off

baseDir = '.';
resultsDir = fullfile(baseDir,'results');

resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles = length(resultFiles);
%%

varNames = {'coeff','method','patient','sample'};
varTypes = {'double','string','int16','int16'};
T = table('Size',[nFiles*4,4], 'VariableTypes',varTypes,'VariableNames',varNames);


resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles = length(resultFiles);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T(ii*4 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample};
    T(ii*4 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample};
    T(ii*4 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample};

    load(fullfile(resultsDir,fileName(1:6)+"_BSC.mat"));
    load(fullfile(baseDir,"band"));
    T(ii*4 + 4,:) = {mean(BSC(:,:,find(band>=3,1)),'all'),"BSC",patient,sample};
end
writetable(T,fullfile(resultsDir,'meanACS.xlsx'))

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
