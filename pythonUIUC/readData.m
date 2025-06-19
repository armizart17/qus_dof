% p1 = load('pythonUIUC\phantom1.mat');
% p2 = load('pythonUIUC\phantom2.mat');
p3 = load('pythonUIUC\phantom3.mat');

elem_pitch = 0.30e-3/2; % because virtual channels so instead of 128 there is 256

ini_frame = 100;
end_frame = 500;
%% PHANTOM 1

% rf = p1.arr_0(:,:,ini_frame:end_frame);
rf = p1.arr_0(:,:,1:400); % better no line at half medium

% 
% bmode_sam = db(abs(hilbert(SAM.rf)));
% bmode_sam = bmode_sam - max(bmode_sam(:));

bmode = my_RF2Bmode(rf);
bmode_range = [-60 0];


% axAxis  = linspace(0, 4E-2, size(rf,1));

maxDepth = 4e-2; c0 = 1540;
fs =  (size(rf,1)-1)*c0/2 / maxDepth;
axAxis = 0:size(rf,1)-1; axAxis = axAxis*1/fs*c0/2;
latAxis = 0:size(rf,2)-1;  latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;

frame = 2; 

figure, 
imagesc(latAxis*1e2, axAxis*1e2, bmode(:,:,frame), bmode_range)
colormap("gray")
h = colorbar;
ylabel(h,'dB')
xlabel('Lateral [cm]'), ylabel('Depth [mm]');
title('Bmode Phantom 1')
set(gca, 'FontSize', 12)

x   = latAxis;
z   = axAxis;
save('./pythonUIUC/rf_phantom1.mat', 'rf', 'x', 'z', 'fs');

%% CREATE VIDEO

rf = p1.arr_0(:,:,1:400);
% rf = p1.arr_0(:,:,400:800);
% rf = p1.arr_0(:,:,800:end);

% 
% bmode_sam = db(abs(hilbert(SAM.rf)));
% bmode_sam = bmode_sam - max(bmode_sam(:));

bmode = my_RF2Bmode(rf);

% --- Parameters ---
video_filename = fullfile('./pythonUIUC/', 'bmode_phantom1.avi');
frame_rate = 20;   % frames per second, adjust as desired
frame_range = 1:size(bmode, 3);

% --- Video Writer Setup ---
v = VideoWriter(video_filename, 'Motion JPEG AVI');
v.FrameRate = frame_rate;
open(v);

figure('Color','w');
for frame = frame_range
    imagesc(latAxis*1e2, axAxis*1e2, bmode(:,:,frame), bmode_range);
    colormap('gray');
    h = colorbar;
    ylabel(h, 'dB');
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    title(['B-mode Phantom 1 - Frame ', num2str(frame)]);
    set(gca, 'FontSize', 12);
    axis image;
    drawnow;
    
    % % Capture frame and write to video
    % frame_image = getframe(gcf);
    % writeVideo(v, frame_image);
end

close(v);
disp(['Video saved as ', video_filename]);

%% PHANTOM 2
rf = p2.arr_0(:,:,ini_frame:end_frame);

bmode = my_RF2Bmode(rf);
bmode_range = [-60 0];

% axAxis = 0:size(rf,1)-1; axAxis = axAxis*1/fs*c0/2;

% axAxis  = linspace(0, 4E-2, size(rf,1));
% axAxis2 =  (0:size(rf,1)-1)*1/fs*c0/2;

maxDepth = 4e-2; c0 = 1540;
fs =  (size(rf,1)-1)*c0/2 / maxDepth;
axAxis = 0:size(rf,1)-1; axAxis = axAxis*1/fs*c0/2;
latAxis = 0:size(rf,2)-1;  latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;
frame = 2;
%
figure, 
imagesc(latAxis*1e2, axAxis*1e2, bmode(:,:,frame), bmode_range)
colormap("gray")
h = colorbar;
ylabel(h,'dB')
xlabel('Lateral [cm]'), ylabel('Depth [mm]');
title('Bmode Phantom 2')
set(gca, 'FontSize', 12)

x   = latAxis;
z   = axAxis;
save('./pythonUIUC/rf_phantom2.mat', 'rf', 'x', 'z', 'fs');

%% CREATE VIDEO

% --- Parameters ---
video_filename = fullfile('./pythonUIUC/', 'bmode_phantom2.avi');
frame_rate = 10;   % frames per second, adjust as desired
frame_range = 1:size(bmode, 3);

% --- Video Writer Setup ---
v = VideoWriter(video_filename, 'Motion JPEG AVI');
v.FrameRate = frame_rate;
open(v);

figure('Color','w');
for frame = frame_range
    imagesc(latAxis*1e2, axAxis*1e2, bmode(:,:,frame), bmode_range);
    colormap('gray');
    h = colorbar;
    ylabel(h, 'dB');
    xlabel('Lateral [cm]');
    ylabel('Depth [cm]');
    title(['B-mode Phantom 2 - Frame ', num2str(frame)]);
    set(gca, 'FontSize', 12);
    axis image;
    drawnow;
    
    % Capture frame and write to video
    frame_image = getframe(gcf);
    writeVideo(v, frame_image);
end

close(v);
disp(['Video saved as ', video_filename]);


%% PHANTOM 3
rf = p3.arr_0(:,:,ini_frame:end_frame);

bmode = my_RF2Bmode(rf);
bmode_range = [-60 0];

% axAxis = 0:size(rf,1)-1; axAxis = axAxis*1/fs*c0/2;

maxDepth = 4e-2; c0 = 1540;
fs =  (size(rf,1)-1)*c0/2 / maxDepth;
axAxis = 0:size(rf,1)-1; axAxis = axAxis*1/fs*c0/2;
latAxis = 0:size(rf,2)-1;  latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;
frame = 2;


figure, 
imagesc(latAxis*1e2, axAxis*1e2, bmode(:,:,frame), bmode_range)
colormap("gray")
h = colorbar;
ylabel(h,'dB')
xlabel('Lateral [cm]'), ylabel('Depth [mm]');
title('Bmode Phantom 3')
set(gca, 'FontSize', 12)

x   = latAxis;
z   = axAxis;
save('./pythonUIUC/rf_phantom3.mat', 'rf', 'x', 'z', 'fs');

%%
dirFigout = './pythonUIUC/';
if (~exist(dirFigout)); mkdir (dirFigout); end
titleFigout = 'Fig';
save_all_figures_to_directory(dirFigout, titleFigout)