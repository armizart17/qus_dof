%%
clear; clc;
close all;

%% Extract structures from .mat
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation\' ...
%     'Liver_24_06_28\set2'];
baseDir = pwd;
%%
for iRef=1:4
name = "r"+iRef+".mat";

input_f = fullfile(baseDir,'refs','left',name);
load(input_f)
rfLeft = rf;

input_f = fullfile(baseDir,'refs','right',name);
load(input_f)
rfRight = rf;

input_f = fullfile(baseDir,'refs','center',name);
load(input_f)
rfCenter = rf;

output_f = fullfile(baseDir,'refs','joined',name);
%%
rf = [rfLeft(:,1:42), rfCenter(:,43:86), rfRight(:,87:end)];
b_mode = db(hilbert(rf));
b_mode = b_mode - max(b_mode(:));

clim_min = -60;
figure();
x = linspace(1, 128, 128); % Same as verasonics image
imagesc(th, r*100-r(1)*100, b_mode);
xlabel('\theta [deg]')
ylabel('\DeltaR [cm]')
%ylim([0 80])
colormap("gray");
colorbar;
clim([clim_min, 0]);

%% To Polar Coordinates
figure();
pcolor(xPolar*1e2,zPolar*1e2,b_mode)
colorbar;
clim([clim_min, 0]);
colormap gray
title('Bmode image')
ylabel('[cm]')
shading interp
axis equal ij tight
%ylim([0 80])
%xlim([-50 50])
%set(gca,'XColor','none','box','off')

%% Saving
save(output_f,'rf','th','r','xPolar','zPolar',"z0Polar",'fs');

%% Auxiliary functions
end
% Get delays

function [t_delay] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)
    
    t_delay = zeros(length(t), n_elements, n_pulses);

    element_pos_x = Trans.ElementPos(:, 1)*wvl; % (x, z) [m] Obtain positions of center of every element
    element_pos_z = Trans.ElementPos(:, 3)*wvl;
    phi = Trans.ElementPos(:, 4);

    for n = 1:n_pulses
        for e = 1:n_elements
            focus = sound_speed*t(:)/2;
            xfocus = element_pos_x(n) + sin(phi(n)) * focus;
            zfocus = element_pos_z(n) + cos(phi(n)) * focus;
          %  keyboard 
            %focus = (20/1000)*ones(size(t));
            %t_delay(:,e,n) = (focus + sqrt((focus - element_pos_z(e)).^2 + (element_pos_x(e)-element_pos_x(n))^2))/sound_speed; 
            t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
            
        end
    end

end