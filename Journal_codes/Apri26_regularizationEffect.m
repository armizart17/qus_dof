list_mu_a  = 10.^[1.5, 2.5 3.5 4.5 5.5 6 8 9];

%% MIDDLE REGULARIZATION (10^3)
% NRMSE 
list_nrmse_a_3dof  = 100*[0.4441, 0.3055, 0.1854, 0.1042, 0.1886, 0.1088, 0.2756, 0.6066];
list_nrmse_a_2dofa = nan(size(list_mu_a)); % a is fixed
list_nrmse_a_2dofb = 100*[0.3788, 0.2894, 0.1784, 0.0956, 0.1963, 0.1014, 0.2717, 0.3673];
list_nrmse_a_2dofn = 100*[0.3789, 0.2895, 0.1785, 0.0957, 0.1316, 0.1008, 0.2627, 0.3333];

list_nrmse_b_3dof  = 100*[0.0122, 0.0122, 0.0122, 0.0122, 0.0379, 0.0122, 0.1695, 0.1793];
list_nrmse_b_2dofa = 100*[0.0199, 0.0199, 0.0199, 0.0199, 0.0331, 0.0199, 0.0199, 0.0199];
list_nrmse_b_2dofb = nan(size(list_mu_a)); % a is fixed
list_nrmse_b_2dofn = 100*[0.0062, 0.0062, 0.0062, 0.0062, 0.0393, 0.0062, 0.1938, 0.3150];

list_nrmse_n_3dof  = 100*[0.0078, 0.0078, 0.0078, 0.0078, 0.1243, 0.0083, 0.8029, 0.8885];
list_nrmse_n_2dofa = 100*[0.0278, 0.0278, 0.0278, 0.0278, 0.0522, 0.0278, 0.0278, 0.0278];
list_nrmse_n_2dofb = 100*[0.0299, 0.0299, 0.0299, 0.0299, 0.1288, 0.0301, 0.4152, 0.5443];
list_nrmse_n_2dofn = nan(size(list_mu_a)); % a is fixed

% MAE
list_mae_a_3dof  = 100*[0.3606, 0.2536, 0.1528, 0.0745, 0.1729, 0.0672, 0.2743, 0.4884];
list_mae_a_2dofa = nan(size(list_mu_a)); % a is fixed
list_mae_a_2dofb = 100*[0.3095, 0.2385, 0.1465, 0.0703, 0.1838, 0.0634, 0.2704, 0.3323];
list_mae_a_2dofn = 100*[0.3095, 0.2386, 0.1467, 0.0704, 0.0925, 0.0636, 0.2609, 0.3333];

list_mae_b_3dof  = 100*[0.0122, 0.0122, 0.0122, 0.0122, 0.0285, 0.0122, 0.1694, 0.1791];
list_mae_b_2dofa = 100*[0.0199, 0.0199, 0.0199, 0.0199, 0.0281, 0.0199, 0.0199, 0.0199];
list_mae_b_2dofb = nan(size(list_mu_a)); % b is fixed
list_mae_b_2dofn = 100*[0.0062, 0.0062, 0.0062, 0.0062, 0.0320, 0.0062, 0.1936, 0.3150];

list_mae_n_3dof  = 100*[0.0078, 0.0078, 0.0078, 0.0078, 0.0975, 0.0078, 0.8021, 0.8873];
list_mae_n_2dofa = 100*[0.0277, 0.0277, 0.0277, 0.0277, 0.0427, 0.0277, 0.0277, 0.0277];
list_mae_n_2dofb = 100*[0.0299, 0.0299, 0.0299, 0.0299, 0.1081, 0.0299, 0.4134, 0.5422];
list_mae_n_2dofn = nan(size(list_mu_a)); % n is fixed

%% LOW REGULARIZATION (10^0.5)

list_mu_a = 10.^[1.5, 2.5, 3.5, 4.5, 5.5, 6, 8, 9];

% MAE
% ---- Parameter a ----
list_mae_a_3dof  = 100*[1.3264, 0.7947, 0.3497, 0.1411, 0.1729, 0.2555, 0.3316, 0.3443];
list_mae_a_2dofa = nan(size(list_mu_a));
list_mae_a_2dofb = 100*[0.3708, 0.2519, 0.1496, 0.1212, 0.1838, 0.2544, 0.3319, 0.3254];
list_mae_a_2dofn = 100*[0.3097, 0.2213, 0.1166, 0.0782, 0.0925, 0.1297, 0.3260, 0.3292];
% ---- Parameter b ----
list_mae_b_3dof  = 100*[0.1088, 0.0901, 0.0609, 0.0385, 0.0285, 0.0266, 0.1275, 0.1749];
list_mae_b_2dofa = 100*[0.0281, 0.0281, 0.0281, 0.0281, 0.0281, 0.0281, 0.0281, 0.0281];
list_mae_b_2dofb = nan(size(list_mu_a));
list_mae_b_2dofn = 100*[0.0256, 0.0256, 0.0259, 0.0261, 0.0320, 0.0422, 0.0976, 0.3418];
% ---- Parameter n ----
list_mae_n_3dof  = 100*[0.4952, 0.3769, 0.2267, 0.1070, 0.0975, 0.1470, 0.6631, 0.8802];
list_mae_n_2dofa = 100*[0.0427, 0.0427, 0.0427, 0.0427, 0.0427, 0.0427, 0.0427, 0.0427];
list_mae_n_2dofb = 100*[0.1017, 0.0940, 0.0857, 0.0818, 0.1081, 0.1408, 0.5632, 0.5778];
list_mae_n_2dofn = nan(size(list_mu_a));

% NRMSE
% ---- Parameter a ----
list_nrmse_a_3dof  = 100*[2.1295, 1.1087, 0.4230, 0.1676, 0.1886, 0.2581, 0.3319, 0.4145];
list_nrmse_a_2dofa = nan(size(list_mu_a));
list_nrmse_a_2dofb = 100*[0.4793, 0.3138, 0.1895, 0.1583, 0.1963, 0.2570, 0.3320, 0.3320];
list_nrmse_a_2dofn = 100*[0.3836, 0.2733, 0.1437, 0.1080, 0.1316, 0.1562, 0.3260, 0.3923];
% ---- Parameter b ----
list_nrmse_b_3dof  = 100*[0.1397, 0.1195, 0.0857, 0.0557, 0.0379, 0.0359, 0.1327, 0.1775];
list_nrmse_b_2dofa = 100*[0.0331, 0.0331, 0.0331, 0.0331, 0.0331, 0.0331, 0.0331, 0.0331];
list_nrmse_b_2dofb = nan(size(list_mu_a));
list_nrmse_b_2dofn = 100*[0.0336, 0.0333, 0.0331, 0.0335, 0.0393, 0.0510, 0.1137, 0.3606];
% ---- Parameter n ----
list_nrmse_n_3dof  = 100*[0.6294, 0.4947, 0.2914, 0.1415, 0.1243, 0.1808, 0.7075, 0.9125];
list_nrmse_n_2dofa = 100*[0.0522, 0.0522, 0.0522, 0.0522, 0.0522, 0.0522, 0.0522, 0.0522];
list_nrmse_n_2dofb = 100*[0.1309, 0.1193, 0.1079, 0.1004, 0.1288, 0.1673, 0.5996, 0.6135];
list_nrmse_n_2dofn = nan(size(list_mu_a));


%% LOW-MIDDLE REGULARIZATION (10^2)
list_mu_a = 10.^[1.5, 2.5, 3.5, 4.5, 5.5, 6, 8, 9];
% MAE
list_mae_a_3dof  = [0.3595, 0.2516, 0.1494, 0.0713, 0.0573, 0.0742, 0.3191, 0.3329];
list_mae_a_2dofa = nan(size(list_mu_a));
list_mae_a_2dofb = [0.3088, 0.2371, 0.1441, 0.0680, 0.0525, 0.0670, 0.3188, 0.3305];
list_mae_a_2dofn = [0.3091, 0.2379, 0.1457, 0.0693, 0.0516, 0.0628, 0.3033, 0.3783];

list_mae_b_3dof  = [0.0122, 0.0122, 0.0122, 0.0122, 0.0122, 0.0122, 0.0587, 0.1987];
list_mae_b_2dofa = [0.0199, 0.0199, 0.0199, 0.0199, 0.0199, 0.0199, 0.0199, 0.0199];
list_mae_b_2dofb = nan(size(list_mu_a));
list_mae_b_2dofn = [0.0062, 0.0062, 0.0062, 0.0062, 0.0062, 0.0074, 0.1218, 0.3418];

list_mae_n_3dof  = [0.0078, 0.0078, 0.0078, 0.0079, 0.0128, 0.0222, 0.6553, 0.8714];
list_mae_n_2dofa = [0.0277, 0.0277, 0.0277, 0.0277, 0.0277, 0.0277, 0.0277, 0.0277];
list_mae_n_2dofb = [0.0299, 0.0299, 0.0299, 0.0299, 0.0299, 0.0353, 0.5023, 0.5852];
list_mae_n_2dofn = nan(size(list_mu_a));

% NRMSE (%)
list_nrmse_a_3dof  = 100*[0.4429, 0.3034, 0.1818, 0.1021, 0.1029, 0.1212, 0.3191, 0.3359];
list_nrmse_a_2dofa = nan(size(list_mu_a));
list_nrmse_a_2dofb = 100*[0.3781, 0.2880, 0.1757, 0.0940, 0.0934, 0.1120, 0.3189, 0.3305];
list_nrmse_a_2dofn = 100*[0.3785, 0.2889, 0.1774, 0.0949, 0.0899, 0.1036, 0.3035, 0.4652];

list_nrmse_b_3dof  = 100*[0.0123, 0.0123, 0.0123, 0.0123, 0.0125, 0.0129, 0.0602, 0.1991];
list_nrmse_b_2dofa = 100*[0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200, 0.0200];
list_nrmse_b_2dofb = nan(size(list_mu_a));
list_nrmse_b_2dofn = 100*[0.0063, 0.0063, 0.0063, 0.0064, 0.0068, 0.0082, 0.1435, 0.3505];

list_nrmse_n_3dof  = 100*[0.0082, 0.0082, 0.0083, 0.0087, 0.0138, 0.0245, 0.6714, 0.8850];
list_nrmse_n_2dofa = 100*[0.0286, 0.0286, 0.0286, 0.0286, 0.0286, 0.0286, 0.0286, 0.0286];
list_nrmse_n_2dofb = 100*[0.0300, 0.0300, 0.0301, 0.0302, 0.0323, 0.0392, 0.5297, 0.6047];
list_nrmse_n_2dofn = nan(size(list_mu_a));

%%
figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 900, 300]); % [x, y, width, height]


% --------- Plot 1: Parameter a ----------
subplot(1,3,1);
semilogx(list_mu_a, list_nrmse_a_3dof, '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_nrmse_a_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_a_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_a_2dofa, '-^', 'LineWidth', 1.5);
    
grid minor;
xlabel('\mu');
ylabel('NRMSE');
title('a');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Plot 2: Parameter b ----------
subplot(1,3,2);
semilogx(list_mu_a, list_nrmse_b_3dof, '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_nrmse_b_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_b_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_b_2dofa, '-^', 'LineWidth', 1.5);

grid minor;
xlabel('\mu');
ylabel('NRMSE');
title('b');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Plot 3: Parameter n ----------
subplot(1,3,3);
semilogx(list_mu_a, list_nrmse_n_3dof, '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_nrmse_n_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_n_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_nrmse_n_2dofa, '-^', 'LineWidth', 1.5);

grid minor;
xlabel('\mu');
ylabel('NRMSE');
title('n');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Global formatting ----------
set(gcf, 'Color', 'w'); % white background

%
figure;
set(gcf, 'Units', 'pixels', 'Position', [50, 100, 900, 300]); % [x, y, width, height]

% --------- Plot 1: Parameter a ----------
subplot(1,3,1);
semilogx(list_mu_a, list_mae_a_3dof,  '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_mae_a_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_a_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_a_2dofa, '-^', 'LineWidth', 1.5);

grid minor;
xlabel('\mu');
ylabel('MAE');
title('a');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Plot 2: Parameter b ----------
subplot(1,3,2);
semilogx(list_mu_a, list_mae_b_3dof,  '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_mae_b_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_b_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_b_2dofa, '-^', 'LineWidth', 1.5);

grid minor;
xlabel('\mu');
ylabel('MAE');
title('b');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Plot 3: Parameter n ----------
subplot(1,3,3);
semilogx(list_mu_a, list_mae_n_3dof,  '-o', 'LineWidth', 1.5); hold on;
semilogx(list_mu_a, list_mae_n_2dofb, '-s', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_n_2dofn, '-d', 'LineWidth', 1.5);
semilogx(list_mu_a, list_mae_n_2dofa, '-^', 'LineWidth', 1.5);

grid minor;
xlabel('\mu');
ylabel('MAE');
title('n');

legend({'3-DoF','2-DoF (b)','2-DoF (n)','2-DoF (a)'}, 'Location','best');

% --------- Global formatting ----------
set(gcf, 'Color', 'w'); % white background