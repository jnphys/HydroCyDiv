function [index, f1] = trackSelection(tracks)
% Function to select trajectories by clicking close to one end and pressing enter after finishing
% Input:
%   - tracks: Cell array containing trajectory data in the format of TrackMate
% Output:
%   - index: Indices of selected trajectories
%   - f1: Figure handle for the plot of selected tracks


%% Plot Tracks
f1 = figure;
distance = nan(size(tracks,1),2);                                                   % Initialize distance matrix
for ii = 1:size(tracks,1)
    plot(tracks{ii,1}(:,2), tracks{ii,1}(:,3), 'linewidth', 1);
    hold all;
    title('Select tracks by clicking at one end and pressing enter in the end ');
    xlabel('x');
    ylabel('y');
end
axis equal
set(gca, 'ydir', 'reverse');

% Get user input by clicking on the figure
[xs, ys] = getpts;

% Initialize index to store selected trajectory indices
index = nan(size(xs,1),1);

% Calculate distances from each end point of each track to the clicked points
for jj = 1:size(xs,1)
    for ii = 1:size(tracks,1)
        distance(ii,1) = sqrt((ys(jj) - tracks{ii}(end,3))^2 + ...
            (xs(jj) - tracks{ii}(end,2))^2);                                        % Distance from the end
        distance(ii,2) = sqrt((ys(jj) - tracks{ii}(1,3))^2 + ...
            (xs(jj) - tracks{ii}(1,2))^2);                                          % Distance from the beginning
    end
    % Find the index of the track with the minimum distance
    [~, index(jj,1)] = find(distance == min(distance(:)));
end

%% Replot Tracks to Confirm Selection
for ii = 1:length(index)
    plot(tracks{index(ii),1}(:,2), tracks{index(ii),1}(:,3), 'k', 'linewidth', 2);
    hold all;
end
pause(1.5)                                                                          % Pause to allow confirmation
end
