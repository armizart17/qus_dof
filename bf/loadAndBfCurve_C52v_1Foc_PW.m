function [MONOFOC, PW] = loadAndBfCurve_C52v_1Foc_PW(fileName, presetName, ax)
% Based on CS code for MonoFoc and PW acquisition BF
% Major changes:
% - Monofocus getBModeC52v change to end both bMode and rf_data
% - Omit zp = zp - Trans.radiusMm*1e-3*(1-cos(phi(1))); (later to do correct
% superposition QUS image+Bmode
% Returns MONOFOC AND PW ACQ (inside rf,bMode,xp,zp,z0p,xr,zr,r0,fs]

% ax = axes(figure);

[~, ~, ext] = fileparts(fileName);
% Whether is a '.vrs' or '.mat' recover the RF Channel Data
if (ext == ".vrs")
    bufferData = vsv.file.readVRSFile(fileName);
    
elseif (ext == ".mat")
    bufferData = load(fileName);
    bufferData = cell2mat(bufferData.RcvData);
end

presetData = load(presetName);
presetData = presetData.preSet;

Trans = presetData.Trans;
Receive = presetData.Receive;
P = presetData.P;
txNum = 1;
RcvZone = P.endDepth;

% Check if the field exists before assignment
if isfield(presetData, 'Refocus') && isfield(presetData.Refocus, 'TXNum')
    txNum = presetData.Refocus.TXNum; % Safe assignment
elseif isfield(presetData, 'TXNum')
    txNum = presetData.TXNum; % Safe assignment
end

centralFreq = Receive(1).demodFrequency*1e6; % Central freq. [Hz]
sampleFreq = Receive(1).decimSampleRate*1e6; % Sample freq. [Hz]

nPulses = P.numRays;
endDepth = P.endDepth;
nElements = Trans.numelements;
nSamples = Receive(1).endSample - Receive(1).startSample + 1;
soundSpeed = 1540; % [m/s]
wvl = soundSpeed/centralFreq;

% Aux variables
t = (0:(nSamples-1))/sampleFreq; % [sec.] time domain 0:T:(N_sample-1)*T
z = soundSpeed*t/2*1e3;
elementPosX = Trans.ElementPos(:, 1)*wvl; % x [m] X center of Trans.element
elementPosZ = Trans.ElementPos(:, 3)*wvl; % z [m] Z center of Trans.element
phi = Trans.ElementPos(:, 4); % phi [rad] angle for every Trans.element
focus = soundSpeed*t(:)/2;
tgcDbCm = 0.0 * (centralFreq * 1e-6) ^ 1;
tgcNpM = tgcDbCm / 8.686 * 100;
r = exp(tgcNpM * z*1e-3); % [m]

% Get reception delays
rxDelays = getRXDelays(elementPosX, elementPosZ, phi, focus, nElements, nPulses, soundSpeed);
rxSamples = round(rxDelays*sampleFreq);
rxSamples(rxSamples<1) = 1;
rxSamples(rxSamples>nSamples) = nSamples;

% Apodization
maxAperSize = 65;
apodAperture = getApodCoeff(maxAperSize, nElements, nPulses);

% Organize and BMode
rfChannel = zeros(nSamples, nElements, nPulses);
plotDepth = floor(endDepth*2*4);
bModeImage = zeros(plotDepth, nPulses); 
rfChannel(:, :, :) = getOrderedData(bufferData, Receive, nSamples, nElements, nPulses, 1);

try
    % Intentar con la asignación de tres dimensiones
    for i = 1:txNum
        rfChannel(:, :, :, i) = getOrderedData(bufferData, Receive, nSamples, nElements, nPulses, i);
    end

    samples2Recon = [nSamples, RcvZone*2*4, 1];
    samples2Recon(2) = [];

    for i = 1:txNum
        rfChannel(samples2Recon(i+1):samples2Recon(i), :, :, end) = rfChannel(samples2Recon(i+1):samples2Recon(i), :, :, i);
    end

    rfChannel = rfChannel(:, :, :, end);
%     rfChannel = getOrderedData(bufferData, Receive, nSamples, nElements, nPulses, 1);
catch
    % Si hay un error, probar con la otra asignación
    rfChannel = getOrderedDataAnt(bufferData, Receive, nSamples, nElements, nPulses, txNum);
end


%%%%%%%%%% Monofoc %%%%%%%%%%
[bModeImage, rf_1foc] = getBModeC52v(rfChannel, rxSamples, apodAperture, nSamples, nPulses, nElements, plotDepth, r);

% For polar grid
param = getparamC52v();
param.fc = centralFreq;
param.fs = sampleFreq;
z2 = (0:size(bModeImage, 1) - 1)/sampleFreq*soundSpeed/2; % [m]

% [xp, zp] = toPolarGrid(size(bModeImage, [1, 2]), z2(end),param);
% vLast
[xp,zp,z0p,xr,zr,r0] = toPolarGridRec(size(bModeImage,[1,2]), z2(end), param);  

% zp = zp - Trans.radiusMm*1e-3*(1-cos(phi(1)));

% [plane, patientCode, frameN] = getAcqData(fileName);
% plotBMode(bModeImage(:, :, end), ax, xp, zp, plane, patientCode, frameN);

rectangleROI = gobjects(0);
rectangleCount = 0;
assignin('base', 'rectangleROI', rectangleROI);
assignin('base', 'rectangleCount', rectangleCount);

% PACKAGING
MONOFOC.rf      = rf_1foc;
MONOFOC.bMode   = bModeImage;
MONOFOC.xp      = xp;
MONOFOC.zp      = zp;
MONOFOC.z0p     = z0p;
MONOFOC.xr      = xr;
MONOFOC.zr      = zr;
MONOFOC.r0      = r0;
MONOFOC.fs      = sampleFreq;


%%%%%%%%%% PWI %%%%%%%%%%

% [rfChannelHP, rfChannelLP] = getOrderedDataPW(bufferData, Receive, nSamples, nElements);
rfChannelPW   = getOrderedDataPW(bufferData, Receive);
rfChannelPW   = double(rfChannelPW);

fNum = 2.5;
maxAperSize = 32;
elementPitch = elementPosX(2) - elementPosX(1);

% Beamforming PW
% Param
param.x_pos = elementPosX';
param.z_pos = elementPosZ' + -1*min(elementPosZ);
param.c = soundSpeed;
param.t0 = 0;
param.fnumber = 2;
param.theta = phi(end) - phi(1);
param.origin = -0.0394;
param.plotear = false;
param.save = false;
param.f_order = 50;

beamPW      = bf2(rfChannelPW(1:plotDepth, :),param,'');
bModePW     = db(abs(hilbert(beamPW)));
bModePW     = bModePW - max(bModePW(:));

% PACKAGING
PW.rf      = beamPW;
PW.bMode   = bModePW;
PW.xp      = xp;
PW.zp      = zp;
PW.z0p     = z0p;
PW.xr      = xr;
PW.zr      = zr;
PW.r0      = r0;
PW.fs      = sampleFreq;

% figure;
% tiledlayout(1, 1, "Padding","tight", 'TileSpacing','compact');
% dbRange = [-55, 0];
% 
% % WideBand
% t1 = nexttile;
% pcolor(xp*1e3, zp*1e3, bModePW);
% clim(dbRange);
% title('Plane Wave - Wideband');
% ylabel('Axial [cm]')
% xlabel('Lateral [cm]')
% shading interp
% axis equal ij tight
% 
% % Colormaps
% colormap(t1,gray);

end
%% Aux Functions
function plotBMode(Image, axes, xp, zp, plane, patientCode, frameN)
        pcolor(axes, xp*1000, zp*1000, Image);
        axes.Color = [0 0 0];
        axes.FontSize = 13;
        cbar = colorbar;
        ylabel(cbar, '[dB]');
        %cbar.Ticks = [];
        clim([-55, 0]);
        colormap gray
        fontSize = 18;
        caption = append('Paciente: ', patientCode, ' - Plano: ', plane, ...
            ' - Fotograma: ', frameN);
        title('B-Mode Multifocal', 'FontSize', fontSize);
        ylabel('[mm]', 'FontSize', 12);
        xlabel('[mm]', 'FontSize', 12);
        shading interp
        axis equal ij tight
end

% Dynamic Focusing Delays
function rx_delays = getRXDelays(element_pos_x, element_pos_z, phi, focus, n_elements, n_pulses, sound_speed)
    rx_delays = zeros(length(focus), n_elements, n_pulses); % length(focus) is same as length(t)
    for n = 1:n_pulses
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        for e = 1:n_elements
            rx_delays(:,e,n) = (focus + ...
                sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
        end
    end
end

% Dynamic aperture matrix
function  apodAperture = getApodCoeff(maxAperSize, nElements, nPulses)
    apodAperture = zeros(nElements, nPulses);
    for i = 1:nPulses
        aperCenter = i;
        halfAperSize = floor(maxAperSize/2);
        aIndex = -halfAperSize:halfAperSize;
        aperture = aperCenter + aIndex;
        aperture = aperture(aperture>=1);
        aperture = aperture(aperture<=nElements);
        apodAperture(aperture, i) = 1;
    end
end

% Organize data
function [rfChannelPW] = getOrderedDataPW(BuffData, Receive)
    rfChannelPW = BuffData(Receive(129).startSample:Receive(129).endSample, :);
end


function rf_channel = getOrderedData(BuffData, Receive, n_samples, n_elements, n_pulses, k)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(n + n_pulses*(k-1)).startSample:Receive(n + n_pulses*(k-1)).endSample, :); % RF Data from Buffer
    end
end

% Organize data Ant
function rf_channel = getOrderedDataAnt(BuffData, Receive, n_samples, n_elements, n_pulses, k)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(k*(n-1)+1).startSample:Receive(k*(n-1)+1).endSample, :); % RF Data from Buffer
    end
end

% Get BMode
function [bModeImage, rf_data] = getBModeC52v(rf_channel, rx_samples, apodAperture, n_samples, n_pulses, n_elements, plotDepth, tgcr)
    rf_data = zeros(n_samples, n_pulses);
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            rf_channel(:, e, n) = rf_channel(rx_samples(:, e, n), e, n);
        end
        % Apodization
        rf_channel(:, :, n) = rf_channel(:, :, n) .* apodAperture(n, :);
        % Summing
        rf_data(:, n) = sum(rf_channel(:, :, n), 2);
    end

    % Depth segmentation
    rf_data = rf_data(1:plotDepth, :);

    % TGC
    tgcr = tgcr(1:plotDepth);
    rf_data = bsxfun(@times, tgcr', rf_data);

    rf_data_padded = zeros(floor(size(rf_data, 1)*1.5), size(rf_data, 2));
    rf_data_padded(1:size(rf_data, 1), :) = rf_data;

    % BMode
    bModeImage = abs(hilbert(rf_data_padded));
    bModeImage = 20*log10(bModeImage);
    bModeImage = bModeImage(1:size(rf_data ,1), :);
    bModeImage = bModeImage-max(bModeImage(:));
end

% Get Acquisition Plane
function [plane, patientCode, frameN] = getAcqData(fileName)
    [~, fName, ~] = fileparts(fileName);
    fParts = strsplit(fName, '_');
    frameN = '1';
    patientCode = fParts{1};
    plane = getPlane(fParts{2});
    if (length(fParts) > 3)
        frameN = fParts{end};
    end       
end


function acqPlane = getPlane(ac_plane)
    switch ac_plane
        case 'LCM'
            acqPlane = 'Línea Clavicular Media';
        case 'IOLAI'
            acqPlane = 'Intercostal Oblicuo Línea Axilar Inferior';
        case 'IHR'
            acqPlane = 'Interfaz Hepatorrenal';
        case 'LHI'
            acqPlane = 'Lóbulo Hepático Izquierdo';
        case 'PL'
            acqPlane = 'Plano Libre';
        otherwise
           acqPlane = 'Plano Legado (antiguo)';
    end
end

function rf=bf2(data, param, file_out)

% Offset
data(1:10,:,:)=0;

fs = param.fs;
c0 = param.c;

freq_vector = [param.fc];
for ii = 1:length(freq_vector)

    fnumber = param.fnumber;

    % param.z_pos = (1:length(data))*(1/fs)*c0/2;
    for rr = 1:size(data,3)
        [rf(:,:,rr),z,x] = bf_planewave_pwc((data(:,:,rr))',fs,fnumber,0,param);
    end

    z=z-param.radius;

    if param.plotear

        % env_ps_data = abs(hilbert(rf(500:end,:,1)));
        env_ps_data = abs(hilbert(rf));
        env_ps_data_norm = env_ps_data./max(env_ps_data(:));
        env_ps_data_norm_comp = 20*log10(env_ps_data_norm);

        % ACTUALIZAR A COORDENADAS POLARES
        % z=z-param.radius; %arriba
        figure;
        % imagesc(x*1e3,z(500:end)*1e3,env_ps_data_norm_comp);
        % imagesc(env_ps_data_norm_comp);
        % s=pcolor(x(500:end,:)*1e3,z(500:end,:)*1e3,env_ps_data_norm_comp);
        s=pcolor(x*1e3,z*1e3,env_ps_data_norm_comp(:,:,1));
        set(s, 'EdgeColor', 'none');
        set(gca,'YDir','reverse')
        colormap gray; colorbar;clim([-50 0]); axis image;
        xlabel('Lateral distance (mm)');ylabel('Depth (mm)');

        parts = split(file_out, '\');
        title(parts{end});
    end

    if param.save
        save(file_out,'c0','fs','rf','x','z');
    end
end
end

function filtrado=filtrar(rf,param)
puntos=[
    0,                  -40;
    0.23229946524064204,-32.188612099644125;
    0.49411764705882355,-39.9288256227758;
    1.058823529411765,  -40.088967971530245;
    1.239786096256685,  -29.89323843416368;
    1.5016042780748666, -23.16725978647687;
    1.7467379679144384, -12.608540925266901;
    1.9828877005347594, -7.206405693950178;
    2.230588235294117,  -5.08185053380783;
    2.4641711229946526, -2.7758007117437717;
    2.7272727272727275, -0.9608540925266906;
    2.9967914438502676, -0.6939501779359434;
    3.2534759358288774, -1.0676156583629899;
    3.497326203208557,  -0.907473309608541;
    3.7347593582887706, -1.8149466192170811;
    3.9914438502673804, -1.6548042704626331;
    4.241711229946525,  -0.10676156583629925;
    4.4727272727272736, -0.8007117437722426;
    4.755080213903744,  -2.9893238434163703;
    5.005347593582888,  -6.245551601423489;
    5.262032085561498,  -12.758007117437725;
    5.486631016042781,  -20.391459074733095
    5.762566844919787,  -28.665480427046266;
    5.993582887700535,  -33.629893238434164;
    6.269518716577541,  -33.629893238434164;
    6.50053475935829,   -33.896797153024906;
    6.7572192513369,    -29.7864768683274;
    6.994652406417113,  -31.654804270462638;
    ];
% figure;plot(puntos(:,1),10.^(puntos(:,2)/20));
f=puntos(:,1)*1e6*2/param.fs;
f(end)=1;
m=10.^(puntos(:,2)/20);
filtro=fir2(param.f_order,f,m);
% figure;plot(filtro);
% figure;freqz(filtro,1,[],param.fs);
filtrado=filter(filtro,1,rf);
end

function [ps_data,z,x] = bf_planewave_pwc(sensor_data_kwave,fs,fnumber,steering_angle,param)

data = sensor_data_kwave';
c=param.c;

origin=param.origin;
z_pos=param.z_pos;
x_pos=param.x_pos;
radius=param.radius;

lenz=size(data,1);
% lenz=256;
lenx=size(data,2);


theta=linspace(0,param.theta,lenx)-param.theta/2;
% rho=c*((1:size(data,1)) /fs)/2;
rho=c*((1:size(data,1)) /fs)/2 +radius;
% rho=linspace(0,84.96e-3,lenz) +radius;
% rho=linspace(0,70e-3,lenz) +radius;

x_arc=radius*sin(theta);

% ayuda
% x=x_arc;
% z=rho-radius;
% save('verasonics\xzayuda','x','z');

[theta,rho] = meshgrid(theta,rho);
[z,x] = pol2cart(theta,rho);



for j=1:lenx %% Iterar en x

    %h1 = waitbar(j/length(x),h1,['Generando líneas... ' num2str(j) '/' num2str(length(x))]);

    for i=1:lenz %% Iterar en z

        %%% Posición del punto en x,z

        % x_point=x(j);
        % z_point=c*(i/fs)/2;
        x_point=x(i,j);
        z_point=z(i,j);
        

        %%% Calcular retardos
        % ang = steering_angle*pi/180;
        % d_ec=z_point*cos(ang)+x_point*sin(ang);
        % tau=(hypot(x_point,z_point)-radius +hypot(z_point+z0+origin-z_pos,x_point-x_pos))/c;
        tx=hypot(x_point,z_point)-radius;
        rx=hypot(z_point-(z_pos-origin),x_point-x_pos);
        tau=(tx + rx)/c;

        %%% Calcular muestras (no se usa interpolación)

        % muestra=round((tau)*fs);
        % muestra(muestra<=0)=1;
        % muestra(muestra>size(data,1))=1;
        % ind=muestra+lenz*(0:size(muestra,2)-1); %Indexación linear
        % values=data(ind);

        muestra=round(tau*fs);
        muestra(muestra<=0)=1;
        muestra(muestra>size(data,1))=1;
        ind=muestra+size(data,1)*(0:size(muestra,2)-1); % indexación linear
        % figure;plot(floor(ind)-size(data,1)*(0:size(muestra,2)-1));
        values=data(ind);
        % values=interp1(data(:),ind,'linear');
        % v1=data(floor(ind));
        % v2=data(ceil(ind));
        % values=v1+ (v2-v1).*(ind-floor(ind));

        %%% Sumar valores solo de la apertura
        %F = 3;
        if isinf(fnumber)
            apodization=1;
        else
            a=tx/(2*fnumber);
            % liminf=x_point-a;
            % limsup=x_point+a;
            % apodization=(tx_positions<limsup & tx_positions>liminf);

            % b=sqrt(rx.^2-tx^2);
            % b=hypot(x_pos(j)-x_pos,z_pos(j)-z_pos);
            b=abs(x_arc(j)-x_arc);
            apodization=b<a;
            % keyboard;

            
        end

        ps_data(i,j,1)=sum(apodization*values');

    end
end
end

% C5-2v Parameters (Verasonics)
function param = getparamC52v()
        param.fc = 3.57e6; % Transducer center frequency [Hz]
        param.kerf = 48e-6; % Kerf [m]
        param.width = 460e-6; % Width [m]
        param.pitch = 508e-6; % Pitch [m]
        param.Nelements = 128;
        param.bandwidth = 79; % Fractional bandwidth [%]
        param.radius = 49.57e-3; % Array radius [m]
        param.height = 13.5e-3; % Elevation height [m]
        param.focus = 60e-3; % Elevation focus [m]
end

function [xPolar, zPolar] = toPolarGrid(siz,zmax,param)
    N = param.Nelements;
    R = param.radius;
    p = param.pitch;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+p,-z0+zmax,siz(1)));
    [xPolar,zPolar] = pol2cart(th,r);
    zPolar = zPolar+z0;
end

function [xPolar,zPolar,z0,xFull,zFull,r0] = toPolarGridRec(siz,zmax,param)                     
    N = param.Nelements;
    R = param.radius;
    p = param.pitch;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+p,-z0+zmax,siz(1)));
    [xPolar,zPolar] = pol2cart(th,r);

    zPolar = zPolar+z0;

    % NEW EMZ
    th_r = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
    r_r = linspace(R+p,-z0+zmax,siz(1));
    
    xFull = th_r; % [deg]
    r0 = r_r(1);
    zFull = (r_r-r0); % [m]

end

