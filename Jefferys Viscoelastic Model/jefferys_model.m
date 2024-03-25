%% Fitting Jeffrey's Model to Single Bead Pulling Data Using fminsearch Optimization Method

% Code was developed by Javad Najafi in Minc's lab at Institut Jacques Monod 
% Last Modified: 2024/03/05

% Description: This code analyzes trajectories of magnetic beads
% pulled through the cytoplasm of sea urchin spherical oocytes from XML files exported by TrackMate.
% The desired trajectory can be selected interactively by clicking on one end of the track.
% Geometrical features, including the bead's initial distance to the egg center and
% the pulling angle between the diagonal line crossing the bead's initial point and the displacement vector, are calculated.
% Jeffrey's viscoelastic model is then fitted to the pulling and release curve, from which viscoelastic parameters are derived.
% fminsearch optimization is very less sensitive to initial values and there is no need to change parameters for evey fit.
%
% Output: 
% - Pulling curve that is scaled displacement to force.
% - Viscoelastic parameters (k, gamma1, and gamma2) obtained from fitting Jeffrey's model.
% - Normalized release curve.
% - Relaxation offset and characteristic time parameters (a and tau) obtained from exponential fit.
% - Distance between egg center and bead initial point
% - Angle of diagonal line crossing the bead's initial position with respect to the x-axis (ang.egg_bead).
% - Angle of bead's displacement vector with respect to the x-axis (ang.bead_displace).
% - Pulling angle that is considered as the difference between the two previous angles (ang.pull).
%
% Additionally, figures of selected track, fits, angles, and workspace variables are exported for future reference.
% Exporting results properly as an Excel file is not compatible with linux and mac systems.


%% Clear Workspace
clearvars;
close all;
clc;


%% Setting the Parameters 
directory = 'C:\Users\Javad\Codes\Bead Model Fit';
name = 'Fluorescence_Cropped_Tracks.xml';                                  % Name of trajectory file
output_name = 'p1_distances.xls';                                          % Name of excel output file
output_folder = 'p1';                                                      % Name of output folder
cd(directory);

start_pull_fr = 22;                                                        % Initial pulling frame 
end_pull_fr = 37;                                                          % Final pulling frame
x_magnet = 329;                                                            % x magnet tip center in pixel                                                      
y_magnet = 1111.5;                                                         % y magnet tip center in pixel

x_egg = 438.7;                                                             % x position of the egg in pixel
y_egg = 147.77;                                                            % x position of the egg in pixel


pixel_size = 0.451;                                                        % um
time_stp = 1;                                                              % second
magnet_radius = 30;                                                        % um
bead_radius = 0.5;                                                         % um

v = @(x) 1286.8 * exp(-x / 4.4) + 88.9 * exp(-x / 135.1);                  % Magnet calibration function / velocity-distance function
viscosity_glycerol = 0.0857;                                               % 80% glycerol viscosity in Pa.s

start1 = [2, 50, 700];                                                     % Starting values of k gamma1 gamma2
start2 = [0.9, 25];                                                        % Starting values of a and tau


%% Loading and Selection of Trajectories
[tracks, ~] = importTrackMateTracks(fullfile(directory, name));
[track_id, f1] = trackSelection(tracks);                                   % Index and replotted figure of selected tracks 
selected_track = tracks{track_id, 1};  


%% Preliminary Trajectory Check
if selected_track(1, 1) > start_pull_fr || selected_track(end, 1) < end_pull_fr + 3
    msg_box = msgbox({'Selected track starts after pulling!',...
        'Not enough points to fit relaxation curve!'}, 'Error', 'error');
    return
end


%% Rotating Selected Track and Visual Check
x = selected_track(:, 2);
y = selected_track(:, 3);

theta = atan2((y(start_pull_fr) - y(end_pull_fr)), (x(start_pull_fr) - x(end_pull_fr)));

rotation_mat = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];

x_rot = x * rotation_mat(1, 1) + y * rotation_mat(1, 2);
y_rot = x * rotation_mat(2, 1) + y * rotation_mat(2, 2);

pull_index = (start_pull_fr:end_pull_fr);
xp = x_rot(pull_index);
yp = y_rot(pull_index);

f2 = figure(2);
subplot(1, 2, 1)
plot(x, y, '.-', 'linewidth', 1);
title('Original Track'); axis equal;
set(gca, 'ydir', 'reverse');
subplot(1, 2, 2)
plot(x_rot, y_rot, '.-', 'linewidth', 1); hold all;
p1 = plot(xp, yp, 'k.-', 'linewidth', 2);
title('Rotated Track'); axis equal; 
legend(p1, 'Pulled Segment');


%% Distance & Angle Calculation
r_egg = [x_egg, y_egg];                                                      
ri_bead = [x(start_pull_fr), y(start_pull_fr)];                            % Bead initial position                                   
rf_bead = [x(end_pull_fr), y(end_pull_fr)];                                % Bead final position 

distance = norm(ri_bead - r_egg) * pixel_size;                             % Distance between egg center and initial bead point

vec.egg_bead = ri_bead - r_egg;                                            % Vector connecting egg center to bead initial point 
ang.egg_bead = -atan2d(vec.egg_bead(2), vec.egg_bead(1));                  % Minus sign is due to flipping y axis
vec.bead_displace = rf_bead - ri_bead;                                     % Bead displacement vector
ang.bead_displace = -atan2d(vec.bead_displace(2), vec.bead_displace(1)); 
  
% Visual checking
f4 = figure(4);
plot(x, y, 'linewidth', 1); hold all;
p2 = plot(x_egg, y_egg, 'ko', 'markersize', 10);
plot([r_egg(1), ri_bead(1)], [r_egg(2), ri_bead(2)], 'r', 'linewidth', 2);
plot([ri_bead(1), rf_bead(1)], [ri_bead(2), rf_bead(2)], 'k', 'linewidth', 2);
title('Pulling Angle Measurement'); axis equal;
legend(p2, 'Egg Center');
set(gca, 'ydir', 'reverse');
axis equal;

% Transforming pulling angle to be between pi and -pi 
ang.pull = ang.egg_bead - ang.bead_displace;                               % Pulling angle of bead's displacement vector
if ang.pull > 180
    ang.pull = ang.pull - 360;
elseif ang.pull < -180
    ang.pull = ang.pull + 360;
end

output_geometry = [distance, ang.egg_bead, ang.bead_displace, ang.pull];   % To be exported


%% Pulling Phase
d = [(x_magnet - x(pull_index)) * pixel_size, (y_magnet - y(pull_index)) * pixel_size];
dist = sqrt(d(:, 1).^2 + d(:, 2).^2) - magnet_radius;                                            % Bead distance to magnet surface

pull_force = 6 * pi * viscosity_glycerol * v(dist) * bead_radius;                                % Magnetic pulling force in pN

x_shift = (-x_rot(start_pull_fr:end) + max(x_rot(start_pull_fr:end))) * pixel_size;              % Flipped and shifted x coordinate

pull_length = end_pull_fr - start_pull_fr + 1;
dx.pulling = x_shift(1:pull_length);
dx.pulling_n = dx.pulling ./ pull_force;                                                         % Scaled displacement 
t.pulling = (0:length(dx.pulling_n) - 1)' * time_stp;                                        

% Fitting Jeffrey's model to the pulling phase
jeffery_model = @(k, gamma1, gamma2, t) (1 - exp(-k * t / gamma1)) / k + t / gamma2;

equ1 = @(start1) norm(jeffery_model(start1(1), start1(2), start1(3), t.pulling) - dx.pulling_n); % Sum of squared errors

options = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 500);
fit1 = fminsearch(equ1, start1, options);

[k, gamma1, gamma2] = deal(fit1(1), fit1(2), fit1(3));

f3 = figure(3);
subplot(1, 2, 1);
plot(t.pulling, dx.pulling_n, 's');
hold all; axis square;
time = linspace(0, t.pulling(end), 1000);
y1 = jeffery_model(k, gamma1, gamma2, time);
plot(time, y1, 'linewidth', 1.5, 'color', 'r');
title('Pulling Phase');
xlabel('t [s]');
ylabel('dx/f [\mum/pN]');

pulling.curves = [t.pulling, dx.pulling_n];
pulling.fitting = [k, gamma1, gamma2];


%% Release Phase
% Fitting exponential model with an offset to the release curve
exp_fit = @(a, tau, t) (1 - a) * exp(-t / tau) + a;

dx.release_n = x_shift(pull_length:end) / x_shift(pull_length);                 % Normalized release displacement

vx_release = diff(x_rot(end_pull_fr:end));
vy_release = diff(y_rot(end_pull_fr:end));
alpha = atan2(vy_release, vx_release);                                          % Angle between two successive steps
correlation = cos(diff(alpha));
ind = find(correlation < 0, 1) + 1;                                             % Last relaxation step when particle moves backward

dx.release_n = dx.release_n(1:ind);
t.release = (0:ind - 1)' * time_stp;

equ2 = @(start2) norm(exp_fit(start2(1), start2(2), t.release) - dx.release_n); % Sum of squared errors

options = optimset('TolX', 1e-7, 'TolFun', 1e-7, 'MaxFunEvals', 1000);
fit2 = fminsearch(equ2, start2, options);

[a, tau] = deal(fit2(1), fit2(2));

subplot(1, 2, 2);
plot(t.release, dx.release_n, 's');
hold all; axis square;
time = linspace(0, t.release(end), 1000);
y2 = exp_fit(a, tau, time);
plot(time, y2, 'linewidth', 1.5, 'color', 'r');
title('Release Phase');
xlabel('t [s]');
ylabel('Normalizd displacement');

release.curves = [t.release, dx.release_n];
release.fitting = [a, tau];


%% Saving Results in Excel
warning('off');
mkdir(output_folder);
cd(fullfile(directory, output_folder));

% First sheet for the pulling phase
headers.pulling = {'t[s]', 'dx/f[um/pN]', 'k[pN/nm]', 'gamma1[pN.s/um]', 'gamma2[pN.s/um]'};
xlswrite(output_name, headers.pulling, 1, 'B2');
xlswrite(output_name, pulling.curves, 1, 'B3');
xlswrite(output_name, pulling.fitting, 1, 'D3');

% Second sheet for the release phase
headers.release = {'t[s]', 'dx/dx(0)', 'a', 'tau[s]'};
xlswrite(output_name, headers.release, 2, 'B2');
xlswrite(output_name, release.curves, 2, 'B3');
xlswrite(output_name, release.fitting, 2, 'D3');

% Third sheet distance and angles 
headers.geometry = {'distance[um]', 'ang1[deg]', 'ang2[deg]', 'pull_ang[deg]'};
xlswrite(output_name, headers.geometry, 3, 'B2');
xlswrite(output_name, output_geometry, 3, 'B3');


%% Saving Workspace and Figures
save(['workspace_variables', '.mat']);
saveas(f1, 'selected_track.jpg');
saveas(f2, 'trajectories.jpg');
saveas(f3, 'fits.jpg');
saveas(f4, 'pulling_angle.jpg');
