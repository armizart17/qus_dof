% ====================================================================== %
% Script for clinical data.
% Created on June 19th, 2024
% ====================================================================== %
% addpath('D:\Lavarello\QUS_codes', 'D:\Lavarello\QUS_codes\BSC estimation\')
%% setup
clear, clc,
close all

baseDir = pwd;
% baseDir = 'LiverAcquisition_24_09_07';

addpath(fullfile(baseDir,'AttenuationUtils'))
addpath(fullfile(baseDir,'UltrasoundUtils'))
addpath(fullfile(baseDir,'BSC estimation'))

sampleDir = fullfile(baseDir,'samples');
refsDir = fullfile(baseDir,'refs','joined');
resultsDir = fullfile(baseDir,'results');
figsDir = fullfile(baseDir,'figures');

[~,~,~] = mkdir(resultsDir);
[~,~,~] = mkdir(figsDir);

sampleFiles = dir(fullfile(sampleDir,'*.mat'));

for iFile = 104:length(sampleFiles)
    % for iFile = 1:length(sampleFiles)
    % samName = "022-02";
    samName = sampleFiles(iFile).name(1:end-4);
    disp(samName)
    
    % %%
    % % Directory containing the .mat files
    % folderPath = 'D:\Lavarello\PROCIENCIA2023\PROCIENCIA2020 Samples\samples_prociencia';
    % 
    % % Get list of .mat files
    % filePattern = fullfile(folderPath, '*.mat');
    % theFiles = dir(filePattern);
    % 
    % % Display the list of files for reference
    % for k = 1:length(theFiles)
    %     fprintf(1, 'File #%d: %s\n', k, theFiles(k).name);
    % end
    % 
    % % Start renaming from the 12th file
    % startIndex = 12;
    % 
    % for k = 1:length(theFiles)
    %     % Get the old file name
    %     baseFileName = theFiles(k).name;
    % 
    %     % Split the file name to construct the new name
    %     [pathStr, name, ext] = fileparts(baseFileName);
    %     nameParts = strsplit(name, '-');
    % 
    %     if str2double(nameParts{1}) >= startIndex
    %         newBaseFileName = sprintf('%03d-%02d.mat', str2double(nameParts{1})-3, str2double(nameParts{2}));
    % 
    %         % Construct the full file path for old and new names
    %         oldFileName = fullfile(folderPath, baseFileName);
    %         newFileName = fullfile(folderPath, newBaseFileName);
    % 
    %         % Rename the file
    %         movefile(oldFileName, newFileName);
    % 
    %         fprintf('Renamed: %s to %s\n', oldFileName, newFileName);
    %     end
    % 
    % end
    % 
    % disp('File renaming completed.');
    
    %%
    
    % % Asumiendo que la tabla se llama volunteersS1
    % % La tabla debe tener las siguientes columnas: 'sex', 'BSCMean', 'BSCStd', 'ACSMean', 'ACSStd'
    % 
    % % Crear una figura
    % figure('Position',[50 50 400 600]);
    % 
    % % Crear el boxchart
    % boxchart(categorical(volunteersS1.sex), volunteersS1.ACSMean);
    % hold on;
    % yline(mean(volunteersS1.ACSMean), '--');
    % % Agregar título y etiquetas a los ejes
    % ylim([0.2 1])
    % % title('Mean ACS, SWIFT');
    % % xlabel('Genero');
    % ylabel('ACS [dB/cm/MHz]');
    % 
    % % Mostrar la grid
    % grid on;
    % set(gca, 'Fontsize', 18);
    % 
    % % Separar los datos en dos grupos según el sexo
    % group_M = volunteersS1.ACSMean(volunteersS1.sex == 'M');
    % group_F = volunteersS1.ACSMean(volunteersS1.sex == 'F');
    % 
    % % Realizar la prueba t de Student
    % [~, p_value] = ttest2(group_M, group_F);
    % 
    % 
    % % Guardar el gráfico si es necesario
    % % saveas(gcf, 'mean_acs_swift_plot.png');
    % % Crear una figura
    % figure('Position',[50 50 400 600]);
    % 
    % % Crear el boxchart
    % boxchart(categorical(volunteersS1.sex), volunteersS1.BSCMean*1e-4);
    % hold on;
    % yline(mean(volunteersS1.BSCMean*1e-4), '--');
    % 
    % set(gca, 'YScale', 'log');
    % % Agregar título y etiquetas a los ejes
    % ylim([1e-4 1e-2])
    % % title('Mean ACS, SWIFT');
    % % xlabel('Sex');
    % ylabel('BSC [1/cm-Sr]');
    % 
    % % Mostrar la grid
    % grid on;
    % set(gca, 'Fontsize', 18);
    % 
    % % Separar los datos en dos grupos según el sexo
    % group_M = volunteersS1.BSCMean(volunteersS1.sex == 'M');
    % group_F = volunteersS1.BSCMean(volunteersS1.sex == 'F');
    % 
    % % Realizar la prueba t de Student
    % [~, p_value] = ttest2(group_M, group_F);
    
    %% ACS estimation hyperparameters
    blocksize = 12;   % Axial block size in wavelengths
    blocklines = 8;   % Num of lines, lateral block size
    overlap_pc = 0.8;
    
    % Bandwidth
    fixedBW = true;
    ratio = db2mag(-30);
    freq_L = 1.5e6; freq_H = 4e6;
    
    % Weight parameters
    ratioCutOff = 10;
    order = 5;
    reject = 0.1;
    extension = 3;
    
    % SWTV
    aSNR = 5; bSNR = 0.09;
    desvMin = 15;
    
    % Reg weights (same as in thyroid cases)
    muBtv = 10^3; muCtv = 10^3;
    muBswtv = 10^2.5; muCswtv = 10^-0.5;
    muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
    muBwfr = 10^3; muCwfr = 10^0.5;
    
    % Plotting constants
    dynRange = [-60,0];
    attRange = [0,1.5];
    bsRange = [-15 15];
    NptodB = log10(exp(1))*20;
    
    % Loading file and variables
    load(fullfile(sampleDir,samName+".mat"));
    fprintf("Loading sample %s \n", samName)
    RcvData = cell2mat(RcvData);
    n_frame = size(RcvData,3); % Frame selector
    RcvData = RcvData(:, :, n_frame); % Select frame of RF Data
    
    % Additional variables
    central_freq = Receive(1).demodFrequency*1e6; % Central frequency of pulse
    fs = Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
    n_pulses = P.numRays; % number of pulses
    n_elements = Trans.numelements; % number of transducer elements
    num_samples = Receive(1).endSample - Receive(1).startSample +1; % samples per channel
    sound_speed = Resource.Parameters.speedOfSound; % [m/s]
    wvl = sound_speed/central_freq; % [m] wavelength in meters
    scalemm2wvl = 1/wvl;
    
    % Initialize variables
    rf_channel = zeros(num_samples , n_elements, n_pulses);
    rx_apods = zeros(1, n_elements, n_pulses);
    rf = zeros(num_samples, n_pulses);
    
    % Organize data
    for n = 1:n_pulses % Iterate through pulses
        rf_channel(:, :, n) = RcvData(Receive(n).startSample:Receive(n).endSample, :);
    end
    
    
    % Acommodate to time delays and rf signals
    focus = 20/1000;
    t = (0:(num_samples-1))/fs; % [sec.] time domain 0:T:(N_sample-1)*T
    [rx_delays] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl);
    
    % Dynamic Aperture
    f_num = 3;
    z = sound_speed*t/2;
    elem_pitch = Trans.spacingMm*1e-3;
    maxAprSz = 32;
    dyn_aperture = zeros(length(z), n_elements, n_pulses);
    for n = 1:n_pulses
        for z_i = 1:length(z)
            a = z(z_i)/(2*f_num);
            hlfAprSz = floor(a / elem_pitch);
            if (hlfAprSz > maxAprSz/2)
                hlfAprSz = floor(maxAprSz / 2);
            end
            a_i = -hlfAprSz: hlfAprSz;    % aperture indices
            fulAprSz = 2*hlfAprSz + 1;
            aper_center = n;
            aper = aper_center + a_i;
            aper = aper(aper>=1);
            aper = aper(aper<=128);
            dyn_aperture(z_i, aper, n) = 1;
        end
    end
    
    
    % Beamforming
    
    % Delay-and-sum
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            % rx_apods(e, n);
            rf_channel(:, e, n) = ...
                interp1(t, rf_channel(:, e, n), rx_delays(:,e, n), 'linear', 0);
        end
        rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
        % Summing
        rf(:, n) = sum(rf_channel(:, :, n), 2);
    end
    
    % B-Mode and coordinates
    b_mode = 20*log10(abs(hilbert(rf)));
    b_mode = b_mode-max(b_mode(:));
    
    param = getparam('C5-2v');
    siz = size(rf);
    z = sound_speed*t/2;
    zmax = z(end);
    R = param.radius;
    p = param.pitch;
    N = param.Nelements;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    
    th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
    r = linspace(R+p,-z0+zmax,siz(1));
    
    % To Polar Coordinates
    [xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);
    
    
    % Selecting ROI
    BmodeFull = db(hilbert(rf));
    BmodeFull = BmodeFull - max(BmodeFull(:));
    
    xFull = th; % [deg]
    r0 = r(1);
    zFull = (r-r0)*1e2; % [cm]
    
    % % % figure('Units','centimeters', 'Position',[5 5 15 15]),
    % % % imagesc(xFull,zFull,BmodeFull,dynRange);
    % % % ylim([0 10])
    % % % colormap gray; clim(dynRange);
    % % % hb2=colorbar; ylabel(hb2,'dB')
    % % % xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    % % % title('Liver')
    % % % 
    % % % confirmation = '';
    % % % while ~strcmp(confirmation,'Yes')
    % % %     rect = getrect;
    % % %     confirmation = questdlg('Sure?');
    % % %     if strcmp(confirmation,'Cancel')
    % % %         disp(rect)
    % % %         break
    % % %     end
    % % % end
    % % close,
    % load(fullfile(resultsDir,samName+".mat"),'rect');
    %% Setting up
    % x_inf = rect(1); x_sup = rect(1)+rect(3);
    % z_inf = rect(2); z_sup = rect(2)+rect(4);
    
    x_inf = xFull(1); x_sup = xFull(end);
    z_inf = zFull(1); z_sup = zFull(end);
    
    dz = (zFull(2) - zFull(1))/100;
    
    % Limits for ACS estimation
    ind_x = x_inf <= xFull & xFull <= x_sup;
    ind_z = z_inf <= zFull & zFull <= z_sup;
    x = xFull(ind_x);
    z = zFull(ind_z);
    sam1 = rf(ind_z,ind_x);
    Bmode = BmodeFull(ind_z,ind_x);
    Bmode = Bmode - max(Bmode(:));
    
    % Wavelength size
    c0 = 1540;
    wl = c0/mean([freq_L freq_H]);   % Wavelength (m)
    
    % Lateral samples
    wx = round(blocklines*(1-overlap_pc));  % Between windows
    nx = blocklines;                 % Window size
    x0 = 1:wx:length(x)-nx;
    x_ACS = x(1,x0+round(nx/2));
    n  = length(x0);
    
    % Axial samples
    wz = round(blocksize*wl*(1-overlap_pc)/dz ); % Between windows
    nz = 2*round(blocksize*wl/dz /2); % Window size
    % nz = 2*round(blocksize*wl/dz /2); % Window size
    L = (nz/2)*dz*100;   % (cm)
    z0p = 1:wz:length(z)-nz;
    z0d = z0p + nz/2;
    z_ACS = z(z0p+ nz/2);
    m  = length(z0p);
    
    % BW from spectrogram
    NFFT = 2^(nextpow2(nz/2)+2);
    band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
    rang = band > freq_L & band < freq_H ;   % useful frequency range
    f  = band(rang)*1e-6; % [MHz]
    p = length(f);
    %%
    % Plot
    [pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
    meanSpectrum = mean(pxx,2);
    meanSpectrum(1) = 0;
    % figure,
    % plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
    % xlim([0, fs/2e6])
    % hold on
    % xline(freq_L/1e6, 'k--')
    % xline(freq_H/1e6, 'k--')
    % hold off
    % xlabel('Frequency [MHz]')
    % ylabel('Magnitude [dB]')
    % saveas(gcf,fullfile(figsDir,"sample"+samName+"_Spectrum.png"))
    %close
    
    fprintf('Frequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
    fprintf('Blocksize in wavelengths: %i\n',blocksize)
    fprintf('Blocksize x: %.2f lines, z: %.2f mm\n',blocklines,nz*dz*1e3)
    fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
    fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);
    
    %% BSC
    % yLimits = [0 12];
    [X,Z] = meshgrid(xFull,zFull);
    roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
    % figure('Units','centimeters', 'Position',[5 5 20 18])
    % imagesc(xFull,zFull,BmodeFull,dynRange); % axis image;
    % title('B-mode')
    % % ylim(yLimits)
    % hold on
    % contour(xFull,zFull,roi,1,'w--')
    % hold off
    % xlabel('Lateral [deg]')
    % ylabel('Axial [cm]')
    % hBm = colorbar;
    % hBm.Label.String = 'dB';
    % hBm.Location = 'westoutside';
    % colormap gray;
    
    n_lay = 3;
    ACS_map = drawlayers(n_lay, BmodeFull, xFull, zFull, xFull(end), [-60 0]);
        save(fullfile(resultsDir,samName+"_ACS.mat"), ...
            'ACS_map')
    % if exist(fullfile(resultsDir,samName+"_ACS.mat"), 'file') == 0
    %     ACS_map = drawlayers(n_lay, BmodeFull, xFull, zFull, xFull(end), [-60 0]);
    %     save(fullfile(resultsDir,samName+"_ACS.mat"), ...
    %         'ACS_map')
    % end
    
    % title("Do not close the figure");
    
    % % Draw a polygon interactively on the image
    % h = drawpolygon('Color','b');
    % 
    % % Wait for the polygon to be finalized
    % wait(h);
    % 
    % % Create a binary mask from the polygon
    % paren_seg_map = createMask(h);
    
    % paren_seg_map(paren_seg_map ~= 1) = 0;
    % paren_seg_map(ACS_map_seg ~= 3) = 0;
    
    % figure;
    % imagesc(paren_seg_map)
    % save(fullfile(resultsDir,samName+"_paren.mat"), ...
    %     'paren_seg_map')
    
    % clear ACS_map;
    % clear paren_seg_map;
    close all;
    
    %%
    
    load(fullfile(resultsDir,samName+"_ACS.mat"),'ACS_map');
    % load(fullfile(resultsDir,samName+"_paren.mat"),'paren_seg_map');
    
    % For looping
    refFiles = dir([refsDir,'\*.mat']);
    Nref = length(refFiles);
    
    % Memory allocation
    samRefs = zeros(size(rf,1),size(rf,2),Nref);
    
    for iRef = 1:Nref
        out = load([refsDir,'\',refFiles(iRef).name]);
        samRefs = out.rf; % Cropping
    end
    
    % Bmoderef = 20*log10(abs(hilbert(samRefs(:,:,1))));
    % Bmoderef = Bmoderef - max(Bmoderef(:));
    % yLimits = [0 12];
    % [X,Z] = meshgrid(xFull,zFull);
    % roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
    % figure('Units','centimeters', 'Position',[5 5 20 18])
    % imagesc(xFull,zFull,Bmoderef,dynRange); % axis image;
    % title('B-mode')
    % ylim(yLimits)
    % hold on
    % contour(xFull,zFull,roi,1,'w--')
    % hold off
    % xlabel('Lateral [deg]')
    % ylabel('Axial [cm]')
    % hBm = colorbar;
    % hBm.Label.String = 'dB';
    % hBm.Location = 'westoutside';
    % colormap gray;
    
    att_ref_dB = 0.54;
    
    % Layers ACS
    n_lay = 3;
    alphaG = 1e-6;                     % Gel [dB/cm/MHz]
    alphaM = 1;                   % Muscle [dB/cm/MHz]
    alphaL = 0.65;
    
    REF.RF = samRefs;
    REF.acs = att_ref_dB;
    
    coeffvals = load("D:\Lavarello\mpUS_Verasonics\BSC_ref_coeff.mat");
    coeffvals = coeffvals.coeffvals;
    BSC_phan_fun = @(f) coeffvals(1)*f.^coeffvals(2);
    REF.BSC_ref = BSC_phan_fun;
    
    pars.fs = fs;
    pars.c = 1540;
    pars.BW = [freq_L freq_H]/1e6;
    pars.nb_lambda_axial = 12;
    pars.nb_lambda_lateral = 8;
    pars.overlap_axial = 0.8;
    pars.overlap_lateral = 0.8;
    pars.P = 2^10;
    
    
    load(fullfile(resultsDir,samName+".mat"),'rect');
    % Setting up
    
    % ROI = [xFull(1) xFull(end) zFull(1) zFull(end)];
    ROI = [rect(1) rect(1)+rect(3) rect(2) rect(2)+rect(4)];
    
    
    % Extra
    ACS_map_ind = ACS_map;
    
    ACS_map = zeros(size(ACS_map));
    
    alpha_layers = [alphaG, alphaM, alphaL];
    
    for i = 1:n_lay
        ACS_map(ACS_map_ind==i) = alpha_layers(i);    
    end
    
    DATA.RF = rf;
    DATA.acs = ACS_map;
    DATA.x = xFull*1e-2;
    DATA.z = zFull*1e-2;
    pars.z_ini = ROI(3)*1e-2;
    pars.z_end = ROI(4)*1e-2;
    pars.x_ini = ROI(1)*1e-2;
    pars.x_end = ROI(2)*1e-2;
    %
    [BSC, band,~,~,~, X_ROI, Z_ROI, SNR] = bsc_estimation(DATA, REF, pars);
    save(fullfile(resultsDir,samName+"_BSC.mat"), ...
        'BSC')
    %%
    % paren_ind = imresize(paren_seg_map, size(BSC, [1 2]), "nearest");
    
    % figure;
    % imagesc(paren_ind);
    % 
    
    % BSC_paren = BSC;
    % BSC_paren(repmat(paren_ind, [1 1 size(band)])==0) = NaN;
    % BSC_paren(71:end,:,:) = NaN;
    
    BSC_paren_single = BSC(:,:,find(band>=3,1)); %BSC_paren(:,:,find(band>=3,1));
    %
    % figure;
    % imagesc(ACS_map); colorbar;
    % figure;
    % imagesc(BSC_paren_single, [1e-5 1e-2]); colorbar;
    % set(gca,'colorscale','log');
    
    bsc_mean = mean(BSC_paren_single(:), 'omitnan');
    bsc_std = std(BSC_paren_single(:), 'omitnan');
    
    filename = 'BSC_data.xlsx';
    % writematrix(["samName" "mean" "std" squeeze(band)], filename, 'Range', sprintf('A%d:GS%d', 1, 1));
    
    % Specify the starting row
    starting_row = iFile+1;  % Example: writing to row i
    
    % Specify the range to write the data
    range = sprintf('A%d:GS%d', starting_row, starting_row);
    
    % Write the data to the specified range in the Excel file
    writematrix([string(samName)  bsc_mean*1e4 bsc_std*1e4], filename, 'Range', range);
    
    % [string(samName)  bsc_mean bsc_std squeeze(mean(BSC_paren,[1 2], 'omitnan'))']
    
    
    % Setting the color scale
    % clim([log10(min(BSC(:))), log10(max(BSC(:)))]);
    % 
    % cb = colorbar();
    % cb.Ruler.Scale = 'log';
    % cb.Ruler.MinorTick = 'on';
    
    % figure;
    % semilogy(band, squeeze(BSC(45,70,:[string(samName)  bsc_mean bsc_std]))')
    % %% Overlay
    % % yLimits = [0, 12];
    % 
    dynRange_BSC = [9e-6 1e-2];
    % 
    % [X,Z] = meshgrid(xFull,zFull);
    % roi = X >= xFull(1) & X <= xFull(end) & Z >= zFull(1) & Z <= zFull(end);
    % 
    % figure('Units','centimeters', 'Position',[5 5 20 12])
    % tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
    % t1 = nexttile();
    % imagesc(xFull,zFull,BmodeFull,dynRange); % axis image;
    % title('B-mode')
    % % ylim(yLimits)
    % hold on
    % % contour(xFull,zFull,roi,1,'w--')
    % hold off
    % xlabel('Lateral [deg]')
    % ylabel('Axial [cm]')
    % hBm = colorbar;
    % hBm.Label.String = 'dB';
    % hBm.Location = 'westoutside';
    % 
    % t2 = nexttile();
    % [~,~,hColor] = imOverlayInterp(BmodeFull,BSC_paren_single,dynRange,dynRange_BSC,0.7,...
    %     X_ROI,Z_ROI,roi,xFull,zFull);
    % title('BSC')
    % subtitle(['BSC at 3MHz: ',num2str(round(1e4*mean(BSC_paren_single,[1 2],'omitnan'),2)),'x10^4 (cm*sr)^{-1}'])
    % axis normal
    % % colorbar off
    % colorbar;
    % set(gca,'colorscale','log');
    % 
    % % ylim(yLimits)
    % hold on
    % % contour(xFull,zFull,roi,1,'w--')
    % hold off
    % % axis off
    % %xlabel('x [cm]')
    % xlabel('Lateral [deg]')
    % colormap(t1,'gray')
    % colormap(t2,'parula')
    % cbh = colorbar ; %Create Colorbar
    % % cbh.Ticks =  logspace(-6,-2,5); %Create 8 ticks from zero to 1
    % % tickLabels = arrayfun(@(x) ['10^{' num2str(log10(x)) '}'], cbh.Ticks, 'UniformOutput', false);
    % % cbh.TickLabels = tickLabels;
    % 
    % % exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_BSC.png"), ...
    % %     'Resolution','300')
    % % close all;
    % 
    figure;
    imagesc(BSC_paren_single, dynRange); colormap gray;
    hBm = colorbar('NorthOutside');
    hBm.Label.String = 'Bmode [dB]';
    set(gca,'Fontsize',16);
    
    figure;
    imagesc(BSC_paren_single, [0 1.5]); colormap turbo;
    hBm = colorbar('NorthOutside');
    hBm.Label.String = 'ACS [dB/cm/MHz]';
    set(gca,'Fontsize',16);
    
    
    figure;
    imagesc(BSC_paren_single, [9e-6 1e-2])
    hBm = colorbar('NorthOutside');
    hBm.Label.String = 'BSC [1/cm-Sr]';
    set(gca,'colorscale','log', 'Fontsize',16);
    
    
    % Plot in cartesian cords
    [TH_acs,R_acs] = meshgrid(-X_ROI*pi/180 + pi/2,Z_ROI/100 + r0);
    [xPolarACS,zPolarACS] = pol2cart(TH_acs,R_acs);
    zPolarACS = zPolarACS + z0Polar;
    
    
    figure('Units','centimeters', 'Position',[5 5 6 6]);
    [ax1,~] = imOverlayPolar(BmodeFull,ones(size(BSC_paren_single)),dynRange,dynRange_BSC,0, ...
        xPolar,zPolar,xPolarACS,zPolarACS);
    set(gca,'colorscale','log');
    yticks(ax1,[4 8 12 16])
    title(ax1,'B-mode')
    xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
    xlim([-9 9]), ylim([0 18])
    hold on
    % contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
    hold off
    exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_Bm.png"), ...
        'Resolution','300')
    
    figure('Units','centimeters', 'Position',[5 5 6 6]);
    [ax1,~] = imOverlayPolar(BmodeFull,BSC_paren_single,dynRange,dynRange_BSC,0.5, ...
        xPolar,zPolar,xPolarACS,zPolarACS);
    yticks(ax1,[4 8 12 16])
    title(ax1,['BSC at 3MHz: ',num2str(round(1e4*mean(BSC_paren_single,[1 2],'omitnan'),2)),'x10^4 (cm*sr)^{-1}'])
    xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
    xlim([-9 9]), ylim([0 18])
    hold on
    % contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
    hold off
    exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_BSCpolar.png"), ...
        'Resolution','300')
    
    pause(0.25)
    close all
end

% close all


% %% Generating Diffraction compensation
% 
% % Generating references
% att_ref = 0.54*f/NptodB; % From 20960001 _ID203458544
% att_ref_map = zeros(m,n,p);
% for jj=1:n
%     for ii=1:m
%         att_ref_map(ii,jj,:) = att_ref;
%     end
% end
% 
% % Windows for spectrum
% % windowing = tukeywin(nz/2,0.25);
% windowing = hamming(nz/2);
% windowing = windowing*ones(1,nx);
% 
% % For looping
% refFiles = dir([refsDir,'\*.mat']);
% Nref = length(refFiles);
% swrap = 0; % Correction factor for phantom data
% 
% % Memory allocation
% Sp_ref = zeros(m,n,p,Nref);
% Sd_ref = zeros(m,n,p,Nref);
% for iRef = 1:Nref
%     out = load([refsDir,'\',refFiles(iRef).name]);
%     samRef = out.rf;
%     samRef = samRef(ind_z,ind_x); % Cropping
%     % figure,imagesc(db(hilbert(samRef)))
%     for jj=1:n
%         for ii=1:m
%             xw = x0(jj) ;   % x window
%             zp = z0p(ii);
%             zd = z0d(ii);
% 
%             sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
%             sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
%             [tempSp,~] = spectra(sub_block_p,windowing,swrap,nz/2,NFFT);
%             [tempSd,~] = spectra(sub_block_d,windowing,swrap,nz/2,NFFT);
% 
%             Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
%             Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
%         end
%     end
% end
% 
% Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
% compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;
% 
% % Liberating memory to avoid killing my RAM
% clear Sp_ref Sd_ref
% 
% %% Spectral log difference
% 
% % Spectrum
% Sp = zeros(m,n,p);
% Sd = zeros(m,n,p);
% for jj=1:n
%     for ii=1:m
%         xw = x0(jj) ;   % x window
%         zp = z0p(ii);
%         zd = z0d(ii);
% 
%         sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
%         sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);
% 
%         [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
%         [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
%         Sp(ii,jj,:) = (tempSp(rang));
%         Sd(ii,jj,:) = (tempSd(rang));
%     end
% end
% 
% % System of eq
% A1 = kron( 4*L*f , speye(m*n) );
% A2 = kron( ones(size(f)) , speye(m*n) );
% b = (log(Sp) - log(Sd)) - (compensation);
% 
% % Optimization constants
% tol = 1e-3;
% clear mask
% mask = ones(m,n,p);
% 
% %% RSLD-TV
% [Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
% BR = reshape(Bn*NptodB,m,n);
% 
% %% SWTV-ACE
% % Calculating SNR
% envelope = abs(hilbert(sam1));
% SNR = zeros(m,n);
% for jj=1:n
%     for ii=1:m
%         xw = x0(jj) ;   % x window
%         zp = z0p(ii);
%         zd = z0d(ii);
% 
%         sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
%         sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
% 
%         temp = [sub_block_p(:);sub_block_d(:)];
%         SNR(ii,jj) = mean(temp)/std(temp);
%     end
% end
% 
% % Calculating weights
% SNRopt = sqrt(1/(4/pi - 1));
% desvSNR = abs(SNR-SNRopt)/SNRopt*100;
% wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));
% 
% % Method
% [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
%     m,n,tol,mask(:),wSNR);
% BSWTV = reshape(Bn*NptodB,m,n);
% CRSWTV = reshape(Cn*NptodB,m,n);
% 
% %% SWIFT
% 
% % First iteration
% [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBwfr,muCwfr,m,n,tol,mask(:));
% bscMap = reshape(Cn*NptodB,m,n);
% 
% % Weight map
% w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
% wExt = movmin(w,extension);
% 
% % Weight matrices and new system
% W = repmat(wExt,[1 1 p]);
% W = spdiags(W(:),0,m*n*p,m*n*p);
% bw = W*b(:);
% A1w = W*A1;
% A2w = W*A2;
% 
% % Second iteration
% [Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
% BSWIFT = reshape(Bn*NptodB,m,n);
% 
% % Weight plot
% % figure('Units','centimeters', 'Position',[5 5 18 4]);
% % tl = tiledlayout(1,3, 'TileSpacing','tight', 'Padding','compact');
% % t2 = nexttile;
% % imagesc(x_ACS,z_ACS,bscMap, [-20 20])
% % colormap(t2,parula)
% % title('BSC map')
% % c = colorbar;
% % c.Label.String = '\Delta BSC [db/cm]';
% % t2 = nexttile;
% % imagesc(x_ACS,z_ACS,w, [0 1])
% % colormap(t2,parula)
% % title('Weights')
% % c = colorbar;
% % c.Label.String = '[a.u.]';
% % t2 = nexttile;
% % imagesc(x_ACS,z_ACS,BSWIFT, attRange)
% % colormap(t2,turbo)
% % title('SWIFT')
% % c = colorbar;
% % c.Label.String = 'ACS [db/cm/MHz]';
% 
% %% SLD and SWIFT fit
% regionMaskAcs = ones(m,n);
% sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*NptodB/sum(regionMaskAcs(:));
% fit1 = f\sldLine;
% fit2 = [f ones(length(f),1)]\sldLine;
% 
% figure('Units','centimeters', 'Position',[5 5 20 10]),
% tiledlayout(1,2),
% nexttile
% plot(f,sldLine),
% hold on,
% plot(f,fit1*f, '--')
% plot(f,fit2(1)*f + fit2(2), '--')
% hold off,
% grid on,
% xlim([0,freq_H*1.1]/1e6),
% ylim([-1 5]),
% xlabel('Frequency [MHz]')
% ylabel('Attenuation [dB/cm]')
% title('Mean SLD')
% legend('Liver',sprintf('ACS ZC = %.2f\n',fit1), ...
%     sprintf('ACS NZC = %.2f\n',fit2(1)), 'Location','northwest')
% 
% sldLine = squeeze(sum(sum(b.*wExt,2),1))/4/L*NptodB/sum(wExt(:));
% fit1 = f\sldLine;
% fit2 = [f ones(length(f),1)]\sldLine;
% nexttile,
% plot(f,sldLine),
% hold on,
% plot(f,fit1*f, '--')
% plot(f,fit2(1)*f + fit2(2), '--')
% hold off,
% grid on,
% xlim([0,freq_H*1.1]/1e6),
% ylim([-1 5]),
% xlabel('Frequency [MHz]')
% ylabel('Attenuation [dB/cm]')
% title('Weighted mean SLD')
% legend('Liver',sprintf('ACS ZC = %.2f\n',fit1), ...
%     sprintf('ACS NZC = %.2f\n',fit2(1)), 'Location','northwest')
% 
% saveas(gcf,fullfile(figsDir,"sample"+samName+"_sldFit.png"))
% close
% %% Overlay
% yLimits = [0, 12];
% meanRsld = mean(BR,'all');
% meanSwtv = mean(BSWTV,'all');
% meanSwift = mean(BSWIFT,'al');
% 
% [X,Z] = meshgrid(xFull,zFull);
% roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
% 
% figure('Units','centimeters', 'Position',[5 5 20 6])
% tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
% t1 = nexttile();
% imagesc(xFull,zFull,BmodeFull,dynRange); % axis image;
% title('B-mode')
% ylim(yLimits)
% hold on
% contour(xFull,zFull,roi,1,'w--')
% hold off
% xlabel('Lateral [deg]')
% ylabel('Axial [cm]')
% hBm = colorbar;
% hBm.Label.String = 'dB';
% hBm.Location = 'westoutside';
% 
% nexttile,
% [~,~,hColor] = imOverlayInterp(BmodeFull,BR,dynRange,attRange,0.7,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('RSLD')
% subtitle(['ACS = ',num2str(meanRsld,2),' dB/cm/MHz'])
% axis normal
% colorbar off
% ylim(yLimits)
% hold on
% contour(xFull,zFull,roi,1,'w--')
% hold off
% % axis off
% %xlabel('x [cm]')
% xlabel('Lateral [deg]')
% 
% nexttile,
% [~,hB,hColor] = imOverlayInterp(BmodeFull,BSWTV,dynRange,attRange,0.7,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('SWTV-ACE')
% subtitle(['ACS = ',num2str(meanSwtv,2),' dB/cm/MHz'])
% colorbar off
% axis normal
% ylim(yLimits)
% hold on
% contour(xFull,zFull,roi,1,'w--')
% hold off
% % axis off
% %xlabel('x [cm]')
% xlabel('Lateral [deg]')
% 
% 
% nexttile,
% [~,hB,hColor] = imOverlayInterp(BmodeFull,BSWIFT,dynRange,attRange,0.7,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('SWIFT')
% subtitle(['ACS = ',num2str(meanSwift,2),' dB/cm/MHz'])
% axis normal
% ylim(yLimits)
% hold on
% contour(xFull,zFull,roi,1,'w--')
% hold off
% xlabel('Lateral [deg]')
% % hColor.Location = 'northoutside';
% % hColor.Layout.Tile = 'northoutside';
% hColor.Label.String = 'ACS [dB/cm/MHz]';
% colormap(t1,'gray')
% % fontsize(gcf,9,'points')
% 
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_Polar.png"), ...
%     'Resolution','300')
% 
% %% Plot in cartesian cords
% [TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS/100 + r0);
% [xPolarACS,zPolarACS] = pol2cart(TH_acs,R_acs);
% zPolarACS = zPolarACS + z0Polar;
% 
% 
% figure('Units','centimeters', 'Position',[5 5 6 6]);
% [ax1,~] = imOverlayPolar(BmodeFull,ones(m,n),dynRange,attRange,0, ...
%     xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
% title(ax1,'B-mode')
% xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
% xlim([-9 9]), ylim([0 18])
% hold on
% contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
% hold off
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_Bm.png"), ...
%     'Resolution','300')
% 
% figure('Units','centimeters', 'Position',[5 5 6 6]);
% [ax1,~] = imOverlayPolar(BmodeFull,BR,dynRange,attRange,0.5, ...
%     xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
% title(ax1,'RSLD')
% xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
% xlim([-9 9]), ylim([0 18])
% hold on
% contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
% hold off
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_rsld.png"), ...
%     'Resolution','300')
% 
% 
% figure('Units','centimeters', 'Position',[5 5 6 6]);
% [ax1,~] = imOverlayPolar(BmodeFull,BSWTV,dynRange,attRange,0.5, ...
%     xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
% title(ax1,'SWTV-ACE')
% xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
% xlim([-9 9]), ylim([0 18])
% hold on
% contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
% hold off
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_swtv.png"), ...
%     'Resolution','300')
% 
% figure('Units','centimeters', 'Position',[5 5 6 6]);
% [ax1,~] = imOverlayPolar(BmodeFull,BSWIFT,dynRange,attRange,0.5, ...
%     xPolar,zPolar,xPolarACS,zPolarACS);
% yticks(ax1,[4 8 12 16])
% title(ax1,'SWIFT')
% xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
% xlim([-9 9]), ylim([0 18])
% hold on
% contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
% hold off
% exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_swift.png"), ...
%     'Resolution','300')
% pause(0.5)
% close all
% 
% 
% %% Saving ACS maps and ROI
% save(fullfile(resultsDir,samName+".mat"), ...
%     'BR','BSWTV','BSWIFT','rect')
% 
% % end



%% Auxiliary functions

% Get delays
function [t_delay] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)

t_delay = zeros(length(t), n_elements, n_pulses);
% (x, z) [m] Obtain positions of center of every element
element_pos_x = Trans.ElementPos(:, 1)*wvl;
element_pos_z = Trans.ElementPos(:, 3)*wvl;
phi = Trans.ElementPos(:, 4);

for n = 1:n_pulses
    for e = 1:n_elements
        focus = sound_speed*t(:)/2;
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + ...
            (xfocus - element_pos_x(e)).^2))/sound_speed;
    end
end

end