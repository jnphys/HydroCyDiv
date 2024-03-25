%% Vesicles Force Study during Aster Centration Using Direct Stoke's Force Balance

% Code was developed by Javad Najafi in Minc's lab at Institut Jacques Monod 
% Last Modified: 2024/02/22  

% Description:
% This code analyzes the trafficking of endomembranes during aster centration to primarily estimate 
% the exerted force by their movement on the aster. It takes as inputs tracks of endomembranes in XML 
% format exported by TrackMate and nucleus/aster track in CSV format. The interpolated aster track,
% which is the longest, is used along with its overlapping segment with each endomembrane track for the analysis.
% 
% The coordinate system is rotated to align the aster centeration along the x-axis, simplifying calculation.
% The entrance point of the aster on the cell periphery should be determined by clicking on the figure to set
% a reference for aster size estimation at each time point. The aster/nucleus coordinate system is used
% to distinguish whether the vesicles are in front or back of the aster.
% 
% Radius of gyration and confinement ratio are calculated as measures of persistence.
% The weighted average force of all vesicles is calculated using Stoke's drag,
% where the weight is the net temporal length of the vesicle's trajectory moving to the right or left.
% Aster size and the average distance between the aster and vesicle are used to determine if a vesicle is inside the aster.
% Displacement speed of vesicles and confinement ratio inside and outside the aster are also calculated. 
% The percentage of vesicles in the front moving inward/backward or in the back moving inward/backward is calculated.
% 
% Output:
% - Weighted force of each vesicle
% - Vesicle displacement speed
% - Vesicle distance to the aster center
% - Tracks confinement ratio
% - Tracks radius of gyration
% - Vesicle position (front/back of aster)
% - Vesicle position (inside/outside of aster)
% - Inward-moving vesicles in the front/back of aster
% - Vesicle percentage in aster’s front/back moving inward/outward
% - Average weighted force exerted by all vesicles


%% Clear Workspace
clearvars;
close all;
clc;

%% Setting Script Parameters 
micronppixel = 0.18;                                                       % Pixel size
fps = 1/1;                                                                 % Frames per second
directory = 'C:\Users\javad\Desktop\endosome';                             % Path of vesicles tracks file
file_name = 'Tracks.xml';                                                  % File name of TrackMate vesicles trajectories in XML format
lag = 30;                                                                  % Skipped frames in aster track to generate time points 
min_overlap_frame = 3;                                                     % Minimum overlapping frames between aster and vesicle
directory_n = directory;                                                   % Path of aster track file
file_name_n = 'Nucleus.csv';                                               % File name of aster track in CSV format
export_dir = 'Weighted Force';                                             % Export folder
cd(directory_n);

R = 0.25;                                                                   % Radius of vesicle in um
eta = 1;                                                                   % Viscosity of cytoplasm in Pa.s.
drag = 6*pi*eta*R;                                                         % Stoke's drag in um.Pa.s.

file_name_exp = 'min_overlap_3frames';                                     % Output file name

%% Loading Vesicle Trajectories
disp('Loading trajectories...');
[p.tracks, md] = importTrackMateTracks(fullfile(directory,file_name));                         % Import TrackMate tracks from XML file
clc; disp('Trajectories loaded!');

for id = 1:length(p.tracks)
    p.tracks_smooth{id,1}(:,1) = p.tracks{id,1}(:,1);                                          % Frame numbers
    p.tracks_smooth{id,1}(:,2) = smooth(p.tracks{id,1}(:,1),p.tracks{id,1}(:,2),1);            % Smooth tracks if desired
    p.tracks_smooth{id,1}(:,3) = smooth(p.tracks{id,1}(:,1),p.tracks{id,1}(:,3),1);
end

%% Loading Track of Nucleus/Aster Center
p.overlap_track = [];

trackload = csvread(fullfile(directory_n,file_name_n),1,4);                                    % Only X and Y coordinates are read

n.track = [(1:lag:(length(trackload)-1)*lag+1)' trackload];                                    % Frames, X, Y
n.track_interp(:,1) = (1:n.track(end,1))';  % Time points
n.track_interp(:,[2,3]) = interp1(n.track(:,1), n.track(:,2:3),n.track_interp(:,1),'spline');  % Interpolated nucleus trajectory

n.track_interp(:,2) = smooth(n.track_interp(:,1),n.track_interp(:,2),1);                       % Smooth tracks if you wish
n.track_interp(:,3) = smooth(n.track_interp(:,1),n.track_interp(:,3),1);

% Finding overlapping frames between vesicle and aster trajectories
track_counter = 1;
for id = 1:length(p.tracks)
    [overlap_inds,~] = find(ismember(p.tracks_smooth{id,1}(:,1), n.track_interp(:,1)) > 0);
    if length(overlap_inds) > min_overlap_frame
        p.overlap_track{1,track_counter} = p.tracks_smooth{id,1}(overlap_inds,(1:3));
        p.overlap_track{2,track_counter} = id;  % Index of overlapping tracks
        track_counter = track_counter+1;
    end
end

%% Original Trajectories
f1 = figure(1); 
hold all;
% Only a fraction of tracks are plotted to keep it clear
for id = 1:size(p.overlap_track,2)/5
    plot(p.overlap_track{1,id}(:,2),p.overlap_track{1,id}(:,3),'linewidth',1.5);
end
plot(n.track_interp(:,2),n.track_interp(:,3),'r','linewidth',2); 

title({'Original Trajectories', 'Click on the beginning of aster track', 'over the cell boundary'});
xlabel('x [pixel]');
ylabel('y [pixel]');
set(gca,'fontsize',9,'ydir','reverse');
daspect([1 1 1]);


%% Rotation of Trajectories and Analysis
% Nucleus entrance point should be determined to estimate the size of aster
particle_num = size(p.overlap_track,2);
p.overlap_rot  = cell(particle_num,1);
p.rot_relative = cell(particle_num,1);
p.vel_relative = cell(particle_num,1);
force_time = nan(particle_num,1);
temp_length = nan(particle_num,1);
weight = nan(particle_num,1);

[n.x0, n.y0] = ginput(1);                                                                                         % Approximate entrance point of aster

n.rot_ang = atan2(n.track_interp(end,3)-n.track_interp(1,3), n.track_interp(end,2)-n.track_interp(1,2));
rot_matrix = [cos(-n.rot_ang), -sin(-n.rot_ang); sin(-n.rot_ang), cos(-n.rot_ang)];
                   
n.track_rot = [n.track_interp(:,1), (rot_matrix*n.track_interp(:,(2:3))')'];

for id = 1:particle_num
    xp = p.overlap_track{1,id}(:,2);
    yp = p.overlap_track{1,id}(:,3);
    p.overlap_rot{1,id} = [p.overlap_track{1,id}(:,1), (rot_matrix*p.overlap_track{1,id}(:,(2:3))')'];            % Pixel
    ind = p.overlap_rot{1,id}(:,1);                                                                               % To find overlapping index in nucleus
    [overlap_inds,~] = find(ismember(n.track_rot(:,1), ind) > 0);
    p.vel_rot_labframe{1,id} = [diff(p.overlap_rot{1,id}(:,2)) diff(p.overlap_rot{1,id}(:,3))]*micronppixel*fps;  % Particle velocity in lab frame um/s 
    p.rot_relative{1,id} = [ind p.overlap_rot{1,id}(:,(2:3))-n.track_rot(overlap_inds,(2:3))]*micronppixel;       % Distance relative to the nucleus in um
    p.vel_relative{1,id} = [diff(p.rot_relative{1,id}(:,2)) diff(p.rot_relative{1,id}(:,3))]*fps;                 % Particle velocity in nucleus frame um/s 

    % Particle is in front or back
    if mean(p.rot_relative{1,id}(:,2)) <= 0      
        stat.front(id,1) = 0;               
    else
        stat.front(id,1) = 1;
    end
    
    % Particle movement direction along x axis
    if p.overlap_rot{1,id}(1,2) < p.overlap_rot{1,id}(end,2)
       stat.vx_sign(id,1) = 1;
    else
       stat.vx_sign(id,1) = -1;
    end    
    
    xmean = mean(xp);  % Pixel
    ymean = mean(yp);  % Pixel
    xmean_n = mean(n.track_interp(ind,2));
    ymean_n = mean(n.track_interp(ind,3));    
    
    stat.gyration_radius(id,1) = sqrt(mean((xp-xmean).^2+(yp-ymean).^2))*micronppixel;                            % Gyration radius in um                    
    stat.particle_nucleus_distance(id,1) = sqrt((xmean-xmean_n)^2+(ymean-ymean_n)^2)*micronppixel;                % Particle-nucleus distance in um
    track_length = sum(sqrt(diff(xp).^2+diff(yp).^2));                                                            % Pixel 
    displacement = sqrt( (xp(1)-xp(end))^2 + (yp(1)-yp(end))^2);                                                  % Pixel
    stat.confinement_ratio(id,1) = displacement/track_length;                                                     % Confinement ratio
    stat.particle_displacement(id,1) = displacement*micronppixel;                                                 % Particle displacement in um    
    
    % Speed summation along x axis for each particle
    p.vx_sum(id,1) = sum(p.vel_rot_labframe{1,id}(:,1)); 
    p.vy_sum(id,2) = sum(p.vel_rot_labframe{1,id}(:,2)); 
        
    % Average speed of particle displacement
    temp_length(id, 1) = (length(xp)-1)/fps;
    stat.mean_vel_labframe(id,1) = displacement*micronppixel/temp_length(id, 1);
    
    % Direct estimation of weighted force for each particle along x axis
    weight(id,1) = abs(sum(sign(p.vel_rot_labframe{1,id}(:,1))));
    force_time(id,1) = drag * (p.vx_sum(id,1)*weight(id, 1));  % pN*s
    
    % Track on average is placed inside/outside of the aster
    aster_size = sqrt((n.x0-n.track_interp(ind(end),2))^2+(n.y0-n.track_interp(ind(end),3))^2)*micronppixel;
    if stat.particle_nucleus_distance(id,1) <= aster_size
        stat.inside_aster(id,1) = 1;         
    else
        stat.inside_aster(id,1) = 0;
    end
end

%% Speed Component Summation Inside/Outside Aster and Weighted Force
vx_inside_aster  = p.vx_sum(stat.inside_aster==1);
vx_outside_aster = p.vx_sum(stat.inside_aster==0);

force_weighted = sum(force_time) / sum(weight);  % pN

%% Percentage Calculation
no_analysed = size(p.overlap_track,2);

stat.front_inward = [];
stat.back_inward = [];
for id=1:no_analysed
    % Front
    if stat.front(id)==1 && stat.vx_sign(id)==-1
        stat.front_inward = [stat.front_inward; 1];
    elseif stat.front(id,1)==1 && stat.vx_sign(id)==1
        stat.front_inward = [stat.front_inward; 0];
    % Back
    elseif stat.front(id)==0 && stat.vx_sign(id)==1
        stat.back_inward = [stat.back_inward; 1];
    elseif stat.front(id)==0 && stat.vx_sign(id)==-1
        stat.back_inward = [stat.back_inward; 0];
    end    
end
percent.front_inward = 100*[sum(stat.front_inward), length(stat.front_inward)-sum(stat.front_inward)]/no_analysed;  % Front inward 1 and front backward 0
percent.back_inward  = 100*[sum(stat.back_inward), length(stat.back_inward)-sum(stat.back_inward)]/no_analysed;     % Back inward 1 and back backward 0
percent.inside_aster = 100*[sum(stat.inside_aster), no_analysed-sum(stat.inside_aster)]/no_analysed;

%% Speed and Confinement Ratio Inside/Outside Aster
v_inside_aster  = stat.mean_vel_labframe(stat.inside_aster==1);
v_outside_aster = stat.mean_vel_labframe(stat.inside_aster==0);

cr_inside_aster  = stat.confinement_ratio(stat.inside_aster==1);
cr_outside_aster = stat.confinement_ratio(stat.inside_aster==0);

%% Plotting the Graphs

%% Rotated & Translated Trajectories
n.track_rot_shift = [n.track_rot(:,2)-n.track_rot(1,2) n.track_rot(:,3)-n.track_rot(1,3)];   

f2 = figure(2);
hold all;
for id = 1:size(p.rot_relative,2)/5
    plot((p.overlap_rot{1,id}(:,2)-n.track_rot(1,2))*micronppixel,...
        (p.overlap_rot{1,id}(:,3)-n.track_rot(1,3))*micronppixel,'linewidth',1.5);
end
plot(n.track_rot_shift(:,1)*micronppixel, n.track_rot_shift(:,2)*micronppixel,'r','linewidth',2);

title('Rotated Trajectories');
xlabel('x [\mum]');
ylabel('y [\mum]');
set(gca,'fontsize',12,'ydir','reverse');
daspect([1 1 1]);
set(gcf, 'color','w');


%% Histogram of Displacement Speed of Vesicles & Weighted Force of Individual Vesicles
particle_speed = mean(stat.mean_vel_labframe)*60;  % um/min

f3 = figure(3);
histogram(stat.mean_vel_labframe*60);
title('Histogram of Deisplacement Speed');
xlabel('displacement speed [\mum/min]','fontsize',14);
ylabel('count','fontsize',14);
set(gca,'fontsize',12);
set(gcf,'color','w');

f4 = figure(4);
histogram(force_time,75,'Normalization','probability'); 
title('Normalized Histogram of Vesicles Weighted Force');
xlabel('Fn*weight [pN]','fontsize',14);
ylabel('probability','fontsize',14);
set(gca,'fontsize',12,'yscale','log');
set(gcf,'color','w');


%% Bar Plots of Front/Back & Inwards/Outwards
f5 = figure(5); 
h = bar([percent.front_inward; percent.back_inward], 0.3,'stacked', 'linewidth', 1.2);

title('Percentage of Vesicles in Front/Back Moving Inwards/Backwards');
set(gca,'xticklabel',{'inwards','outwards'},'linewidth',1);
h(1).FaceColor = 'c';
h(2).FaceColor = 'g';
legend('front', 'back');
legend boxoff
ylabel('frequency%');
ylim([0 100]);
set(gcf,'color','w');
axis square

%% Export Results as Excel File in Export Directory
% Check if export directory exists, if not, create it
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
end
cd(export_dir);

warning('off');
% First sheet for weighted force and speeds
headers.sh1 = {'F_{n}t_{n}[pN.s]', 'speed[um/s]'};
output.sh1 = [force_time, stat.mean_vel_labframe];
xlswrite(file_name_exp, headers.sh1, 'B2');
xlswrite(file_name_exp, output.sh1, 1, 'B3');

% Second sheet for statistics and possible correlations
headers.sh2 = {'distance[um]', 'speed[um/s]', 'cr', 'gyr[um]', 'front', 'inside', 'front_in', 'back_in'};
output.sh2 = [stat.particle_nucleus_distance,...
    stat.mean_vel_labframe,...
    stat.confinement_ratio,...
    stat.gyration_radius,...
    stat.front,...
    stat.inside_aster];
xlswrite(file_name_exp, headers.sh2 ,2 ,'B2');
xlswrite(file_name_exp, output.sh2, 2, 'B3');
xlswrite(file_name_exp, stat.front_inward, 2, 'H3');
xlswrite(file_name_exp, stat.back_inward, 2, 'I3');

% Third sheet for percentages
headers.sh3 = {'front_in', 'front_out', 'back_in', 'back_out', 'inside', 'outside'};
output.sh3 = [percent.front_inward, percent.back_inward, percent.inside_aster];
xlswrite(file_name_exp, headers.sh3, 3, 'B2');
xlswrite(file_name_exp, output.sh3, 3, 'B3');

% Forth sheet for force
headers.sh4 = {'weighted force [pN]'};
xlswrite(file_name_exp, headers.sh4, 4, 'B2');
xlswrite(file_name_exp, force_weighted, 4, 'B3');

%% Saving all Workspace
save('workspace_variables.mat');

%% Saving the Figures
if ~exist('plots', 'dir')
    mkdir('plots');
end
cd('plots');

savefig(f1,'trajectories_original.fig');
savefig(f2,'trajectories_rotated.fig');
savefig(f3,'histogram_average_speed.fig');
savefig(f4,'histogram_weighted_force.fig');
savefig(f5,'inwards_outwards_statistics.fig');


