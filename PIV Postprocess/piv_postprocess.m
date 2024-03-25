%% Plotting and Analysis of PIV Results inside ROI

% Code was developed by Javad Najafi in Minc's lab at Institut Jacques Monod 
% Last Modified: 2024/02/27

% Description:
% This code generates high-quality graphs depicting the Particle Image Velocimetry (PIV) analysis of the sea
% urchin cytoplasm when microtubule aster is manipulated by magnetic tweezers and positioned near the cell cortex.
% 
% It loads the PIV analysis data exported from PIVlab, performs artifact correction to fill the voids
% present in the original vector field after removing erroneous vectors, and generates a proper mask. 
% These masked regions of vector fields are excluded from the plotted heatmaps and velocity analysis.
% 
% The code also calculates the average velocity components over a defined interval by applying the mask
% to each frame of the smoothed vector field. Additionally, it computes the average x-component of the 
% velocity in an interactively drawn region of interest (ROI) after the aster reaches the cortex.
% 
% Outputs:
% - heatmap of average x-component of velocity
% - Streamlines of the average velocity
% - Temporal evolution of x-component of velocity in ROI 


%% Clear Workspace
clearvars;
close all;
clc;

%% Loading mat File of PIVLab Output and Setting the Parameters
directory = 'F:\Aster Pulling Inward Flow\20191009\1\DIC images';
name = 'piv_analysis';
load(fullfile(directory, [name, '.mat']));
cd(directory)

fps = 1/3;                                                                 % Frame per second
pixel = 0.216;                                                             % Pixel size in micron
f_frame = 180;                                                             % First frame number for averaging interval
l_frame = 200;                                                             % Last frame number for averaging interval
last_pull_frame = 80;

%% Generating Masks for Each Frame Using Original Vector Field Before Smoothing
mask = cell(size(u_original));
mask_final = cell(size(u_original));

for ii = 1:length(u_original)
    mask{ii,1} = u_original{ii,1} ./ u_original{ii,1};
    mask{ii,1}(isnan(mask{ii,1})) = 0;                                     % Mask to binary matrix
    mask_final{ii,1} = imclose(mask{ii,1}, strel('disk', 1));              % Closing the holes connected to the boundary
    mask_final{ii,1} = bwareaopen(mask_final{ii,1}, 30, 8);                % Removing holes outside of egg
    mask_final{ii,1} = double(mask_final{ii,1});                           % Logical to double matrix
    mask_final{ii,1}(~mask_final{ii,1}) = nan;                             % Replace with nan to skip in the plots
end

%% Averaging Velocity Components, Speed, and Vorticity
masked.vx = nan([(size(mask{1,1})), length(mask)]);
masked.vy = nan([(size(mask{1,1})), length(mask)]);
masked.vmag = nan([(size(mask{1,1})), length(mask)]);

% Applying masks and converting to 3d matrix
for ii = 1:length(v_smoothed)
    masked.vx(:,:,ii) = u_smoothed{ii,1} .* mask_final{ii,1} * pixel * fps;
    masked.vy(:,:,ii) = v_smoothed{ii,1} .* mask_final{ii,1} * pixel * fps;
    masked.vmag(:,:,ii) = velocity_magnitude{ii,1} .* mask_final{ii,1} * pixel * fps;
end

%% Plotting Heatmap of Speed 
% Averaging in defined temporal interval 
mean.vx = nanmean(masked.vx(:,:,f_frame:l_frame), 3);        
mean.vy = nanmean(masked.vy(:,:,f_frame:l_frame), 3); 
mean.vmag = nanmean(masked.vmag(:,:,f_frame:l_frame), 3);  

f1 = figure(1);
% Heatmap of Velocity Component along x-axis
pcolor(x{1,1}, y{1,1}, mean.vx); hold on;
title('Heatmap and Streamlines of Average Velocity');
colormap jet;
shading interp;
q = colorbar;
q.Ruler.Exponent = -3;  
ylabel(q, 'Velocity component x (\mum/s)','fontsize',16);

% Streamlines
g = streamslice(x{1,1}, y{1,1}, mean.vx, mean.vy, 2);
set(g, 'color', 'k', 'linewidth', 1.1);

set(gcf, 'color', 'w');
set(gca, 'ydir', 'reverse', 'xtick', [], 'ytick', [], 'color', 'w', 'fontsize', 11);
axis equal tight off
zoom(1.2);

%% X-Component Average in ROI 
f2 = figure(2); hold on; 
% Drawing vector field
quiver(x{1,1}, y{1,1}, mean.vx, mean.vy, 1.25, 'linewidth', .5, 'color', 'r');
title('Draw ROI you want to average speed...');
axis equal;

% Drawing ROI and Selecting Data
roi = impoly();                                                            % Draw ROI                                          
nodes = getPosition(roi);                                                  % Position of the nodes
inside = inpolygon(x{1,1}, y{1,1}, nodes(:,1), nodes(:,2));                % Mask of the points inside

% Generate Mask of ROI
inside = double(inside);
inside(~inside) = nan;

% Applying Mask and Replotting Vector Field in ROI
mean.vx_roi = mean.vx .* inside;
mean.vy_roi = mean.vy .* inside;
quiver(x{1,1}, y{1,1}, mean.vx_roi, mean.vy_roi, 1.25, 'linewidth', .3, 'color', 'k');

% Applying Mask to Each Frame
vx_roi_mean = nan(size(u_smoothed,1) - last_pull_frame + 1, 2);
for ii = last_pull_frame:size(u_smoothed,1)
    vx_roi = masked.vx(:,:,ii) .* inside;
    vx_roi_mean(ii - last_pull_frame + 1, :) = [ii, nanmean(vx_roi(:))];
end

f3 = figure(3);
plot(vx_roi_mean(:,1), vx_roi_mean(:,2));
title('Temporal Evolution of v_{x} inside ROI');
xlabel('Frame after pulling');
ylabel('Average v_{x} (\mum/s)');
