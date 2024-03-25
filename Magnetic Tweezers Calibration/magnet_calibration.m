%% Obtaining Main Calibration Curve of Magnet by Magnetic Bead Pulling 

% Code was developed by Javad Najafi in Minc's lab at Institut Jacques Monod 
% Last Modified: 2024/01/19 

% Description: This code analyzes trajectories of magnetic beads pulled through a test viscous fluid from XML files
% exported by TrackMate. The desired trajectory can be selected interactively by clicking on one end of the track.
% A double exponential function is fitted to the data of all trajectories to get the calibration coefficients.
%
% Output: 
% - Pulling curves for individual tracks (dist_speed)
% - Coefficient of calibration curve (a1, b1, a2, b2)


%% Clearing workspace
clearvars;
close all;
clc;


%% Input Parameters
directory = 'C:\Users\Javad\Codes\Magnetic Tweezers Calibration';
magnetic_tracks = 'pulled beads track inside test viscous fluid.xml';      % Name of trajectory file
output_name = 'workspace_variables';                                       % Name of output file
output_folder = 'output';                                                  % Name of output folder
cd(directory);

magnet_radius = 30;                                                        % Magnet tip radius in micrometers
pixel_size = 1;                                                            % Pixel size in micrometers
frame_rate = 1/3;                                                          % Frame rate as frames per second
x_magnet = 86.9;                                                           % x position of magnet tip center (pixel)
y_magnet = 333.6;                                                          % y position of magnet tip center (pixel)
drift = [0 0];                                                             % Drift speed of tracers in micrometers per second


%% Plotting Tracks and Selecting Pulled Ones
[tracks, ~] = importTrackMateTracks(fullfile(directory, magnetic_tracks));

f1 = figure(1);
for ii = 1:size(tracks, 1)
    plot(tracks{ii, 1}(:, 2), tracks{ii, 1}(:, 3));
    hold on;
    plot(x_magnet, y_magnet, 'x');
end
title('Select Pulled Tracks and Press Enter...', 'FontSize', 12);
set(gca, 'YDir', 'reverse');
axis equal; 

[xs, ys] = getpts;   

index = nan(size(xs, 1));
for jj = 1:size(xs, 1)
    distance = nan(size(tracks, 1));
    for ii = 1:size(tracks, 1)
        distance(ii, 1) = sqrt((ys(jj) - tracks{ii}(end, 3))^2 + (xs(jj) - tracks{ii}(end, 2))^2); % Distance from the end
        distance(ii, 2) = sqrt((ys(jj) - tracks{ii}(1, 3))^2 + (xs(jj) - tracks{ii}(1, 2))^2);     % Distance from the beginning
    end
    [index(jj, 1), ~] = find(distance == min(distance(:)));                                        % Index of selected tracks
end

dist_speed = cell(length(index), 1);
for ii = 1:length(index)
    n = index(ii);
    x = tracks{n, 1}(:, 2);
    y = tracks{n, 1}(:, 3);
    
    vx = diff(x) * pixel_size / frame_rate - drift(1);
    vy = diff(y) * pixel_size / frame_rate - drift(2);
    d = sqrt((x - x_magnet).^2 + (y - y_magnet).^2) - magnet_radius;                                % Distance from magnet surface 
    v = sqrt(vx.^2 + vy.^2);
    
    dist_speed{ii, 1} = [d(1:end-1), v];
    
    plot(tracks{n, 1}(:, 2), tracks{n, 1}(:, 3), 'k', 'LineWidth', 2);                              % Replot tracks to ensure from the selection
    hold on;
end   
pause(4);

f2 = figure(2);
for ii = 1:size(dist_speed, 1)
    plot(dist_speed{ii, 1}(:, 1), dist_speed{ii, 1}(:, 2));
    hold on;
end
title('Curves for Individual Tracks');
xlabel('Distance [\mum]');
ylabel('Speed [\mum/s]');


%% Fitting Double Exponential Function to All Tracks
dist_speed_mat = cell2mat(dist_speed(:));
dist = dist_speed_mat(:, 1);
speed = dist_speed_mat(:, 2);

% Double exponential fit function and initial values
double_exp = @(a1, b1, a2, b2, x) a1 * exp(-dist / b1) + a2 * exp(-dist / b2);
start = [500, 100, 100, 500]; 

equ = @(start) norm(double_exp(start(1), start(2), start(3), start(4), dist) - speed);              % Sum of squared errors
fit1 = fminsearch(equ, start);

[a1, b1, a2, b2] = deal(fit1(1), fit1(2), fit1(3), fit1(4));

f3 = figure(3);
scatter(dist, speed, 50, '.');
hold on;
fplot(@(x) a1 * exp(-x / b1) + a2 * exp(-x / b2), [50 900], 'r', 'LineWidth', 1.5);
title('Double Exponential Fit to Bead Pulling Data');
xlabel('Distance [\mum]');
ylabel('Speed [\mum/s]');
xlim([50 900]);
box on;


%% Uncomment to Save all Workspace and Fit Results 
% mkdir(output_folder);
% cd(fullfile(directory, output_folder));

% save([output_name, '.mat']);
% saveas(f3, 'double_expo_fit.jpg');

