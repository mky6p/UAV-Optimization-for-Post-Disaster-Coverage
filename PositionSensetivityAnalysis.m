%% UAV Optimization with Greedy Algorithm (Improved with Gradient-Based Threshold)
clc;
clear;

% Load building and UAV data
buildings_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Buildings data\Building_joplin-tornado_00000139.xlsx";
uav_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Results\Result139\UAV_Inverse_FSPL_Random_Positions_and_Frequencies.xlsx";

buildings_data = readtable(buildings_file);
uav_data = readtable(uav_file);

% Constants
altitude = 150; % UAV altitude in meters
x_min = -300; x_max = 500; % X-axis range
y_min = -300; y_max = 500; % Y-axis range
grid_resolution = 300; % Number of grid points per axis
c = 3e8; % Speed of light in m/s

% Minimum distance between UAVs
min_distance = 30; % Adjust as needed

% Define grid and Gaussian kernel
x_range = linspace(x_min, x_max, grid_resolution);
y_range = linspace(y_min, y_max, grid_resolution);
[X, Y] = meshgrid(x_range, y_range);
Z = zeros(size(X)); % Initialize intensity grid

% Extract building coordinates and Gaussian kernel parameters
building_x = buildings_data.X_meters_; 
building_y = buildings_data.Y_meters_; 
gaussian_values = buildings_data.Intensity; 
sigma_values = buildings_data.Sigma_meters_; 

% Compute Gaussian kernel contributions
gaussian_kernel_3d = @(X, Y, x0, y0, sigma, intensity) ...
    intensity * exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2));

for i = 1:length(building_x)
    Z = Z + gaussian_kernel_3d(X, Y, building_x(i), building_y(i), sigma_values(i), gaussian_values(i));
end

% Normalize the Gaussian kernel
Z = Z / max(Z(:));

% --- Gradient-Based Threshold Optimization ---
% Step 1: Compute the gradient of the demand intensity (Z)
[grad_x, grad_y] = gradient(Z, x_range, y_range);
grad_magnitude = sqrt(grad_x.^2 + grad_y.^2);

% Step 2: Normalize the gradient magnitude
grad_magnitude = grad_magnitude / max(grad_magnitude(:));

% Set a gradient threshold to identify significant hotspots
grad_threshold = 0.3; % Adjust as needed
hotspot_mask = grad_magnitude > grad_threshold;

% Step 3: Dynamically adjust the intensity threshold
hotspot_intensities = Z(hotspot_mask);
if ~isempty(hotspot_intensities)
    high_intensity_threshold = prctile(hotspot_intensities, 75); % 75th percentile
else
    high_intensity_threshold = 0.5; % Fallback
end

disp(['Computed High Intensity Threshold: ', num2str(high_intensity_threshold)]);

% Step 4: Select candidate positions
high_intensity_positions = [X(Z > high_intensity_threshold), Y(Z > high_intensity_threshold)];
available_positions = high_intensity_positions;

% Extract UAV positions and frequencies from the file
uav_positions_initial = [uav_data.X_Position_m, uav_data.Y_Position_m]; % Initial positions from file
frequencies = uav_data.Frequency_Hz;
num_uavs = length(frequencies);

% Ensure the number of UAVs matches the file data
if size(uav_positions_initial, 1) ~= num_uavs
    error('Number of UAVs in position data does not match the number of frequencies.');
end

% Initialize UAV positions for greedy algorithm with file data
uav_positions = uav_positions_initial;

% Greedy placement of UAVs (refining initial positions)
for i = 1:num_uavs
    best_coverage_local = -Inf;
    best_position_local = [];
    
    temp_coverage = zeros(size(available_positions, 1), 1);
    temp_positions = zeros(size(available_positions, 1), 2);
    
    parfor j = 1:size(available_positions, 1)
        x_candidate = available_positions(j, 1);
        y_candidate = available_positions(j, 2);
        
        too_close = false;
        for k = 1:i-1
            if norm([x_candidate, y_candidate] - uav_positions(k, :)) < min_distance
                too_close = true;
                break;
            end
        end
        
        if ~too_close
            temp_pos = [uav_positions(1:i-1, :); x_candidate, y_candidate];
            temp_coverage(j) = compute_total_coverage(X, Y, Z, temp_pos, frequencies(1:i), altitude, c);
            temp_positions(j, :) = [x_candidate, y_candidate];
        end
    end
    
    [best_coverage_local, best_idx] = max(temp_coverage);
    if best_coverage_local > -Inf
        best_position_local = temp_positions(best_idx, :);
    end
    
    if ~isempty(best_position_local)
        uav_positions(i, :) = best_position_local;
    end
    
    available_positions(ismember(available_positions, best_position_local, 'rows'), :) = [];
end

% Display results
disp('Greedy Algorithm Optimal UAV Positions (starting from file):');
disp(uav_positions);

%% Visualization
% Fullscreen high-res figure
figure('Units','normalized','OuterPosition',[0 0 1 1]); 

markerSize = 80;
fontSize = 20;
xLimits = [min(x_range), max(x_range)];
yLimits = [min(y_range), max(y_range)];

% ---- Create subplots manually with tight control ----
t = tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact'); % No spacing

% 1. Gaussian Kernel Heatmap
nexttile;
imagesc(x_range, y_range, Z);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Gaussian Kernel Heatmap', 'FontSize', fontSize);

% 2. Optimal UAV Positions
nexttile;
imagesc(x_range, y_range, Z);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Optimal UAV Positions', 'FontSize', fontSize);
hold on;
scatter(uav_positions(:, 1), uav_positions(:, 2), markerSize, 'bx', 'LineWidth', 2);
hold off;

% 3. Signal Strength Distribution
[~, normalized_signal] = compute_total_coverage(X, Y, Z, uav_positions, frequencies, altitude, c);
nexttile;
imagesc(x_range, y_range, normalized_signal);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'parula');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Signal Strength Distribution', 'FontSize', fontSize);
hold on;
scatter(uav_positions(:, 1), uav_positions(:, 2), markerSize, 'rx', 'LineWidth', 2);
hold off;

% ---- Ensure no legend appears ----
legend off

% ---- Export high-resolution, clean figure ----
exportgraphics(gcf, 'Result_Final_NoLegend_TightLayout.png', 'Resolution', 300);

%% Sensitivity Analysis
% Define parameter ranges
min_distances = [20, 25, 30, 35, 40];
altitudes = [100, 125, 150, 175, 200];
percentiles = [50, 60, 70, 75, 80, 90];

coverages_min_distance = zeros(length(min_distances), 1);
coverages_altitude = zeros(length(altitudes), 1);
coverages_percentile = zeros(length(percentiles), 1);

% Store UAV positions for each parameter value
positions_min_distance = zeros(num_uavs, 2, length(min_distances));
positions_altitude = zeros(num_uavs, 2, length(altitudes));
positions_percentile = zeros(num_uavs, 2, length(percentiles));

% 1. Vary Minimum Distance
for m = 1:length(min_distances)
    min_distance = min_distances(m);
    altitude = 150;
    percentile = 75;
    [final_coverage, uav_pos] = run_optimization(X, Y, Z, hotspot_mask, min_distance, altitude, percentile, num_uavs, frequencies, uav_positions_initial, c);
    coverages_min_distance(m) = final_coverage;
    positions_min_distance(:,:,m) = uav_pos;
    fprintf('Min Distance = %d m, Coverage = %.4f\n', min_distance, final_coverage);
end

% 2. Vary Altitude
for a = 1:length(altitudes)
    altitude = altitudes(a);
    min_distance = 30;
    percentile = 75;
    [final_coverage, uav_pos] = run_optimization(X, Y, Z, hotspot_mask, min_distance, altitude, percentile, num_uavs, frequencies, uav_positions_initial, c);
    coverages_altitude(a) = final_coverage;
    positions_altitude(:,:,a) = uav_pos;
    fprintf('Altitude = %d m, Coverage = %.4f\n', altitude, final_coverage);
end

% 3. Vary Percentile
for p = 1:length(percentiles)
    percentile = percentiles(p);
    altitude = 150;
    min_distance = 30;
    [final_coverage, uav_pos] = run_optimization(X, Y, Z, hotspot_mask, min_distance, altitude, percentile, num_uavs, frequencies, uav_positions_initial, c);
    coverages_percentile(p) = final_coverage;
    positions_percentile(:,:,p) = uav_pos;
    fprintf('Percentile = %d, Coverage = %.4f\n', percentile, final_coverage);
end

% --- Save UAV Positions ---
% Create output directory if it doesn't exist
output_dir = 'C:\Users\mky6p\OneDrive - University of Missouri\Mill\Results\Result139\Sensitivity_Positions';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% 1. Save Minimum Distance Positions
min_distance_table = table();
for m = 1:length(min_distances)
    pos = positions_min_distance(:,:,m);
    temp_table = table((1:num_uavs)', repmat(min_distances(m), num_uavs, 1), pos(:,1), pos(:,2), ...
        'VariableNames', {'UAV_Index', 'Min_Distance', 'X_Position', 'Y_Position'});
    min_distance_table = [min_distance_table; temp_table];
end
writetable(min_distance_table, fullfile(output_dir, 'UAV_Positions_Min_Distance.xlsx'));

% 2. Save Altitude Positions
altitude_table = table();
for a = 1:length(altitudes)
    pos = positions_altitude(:,:,a);
    temp_table = table((1:num_uavs)', repmat(altitudes(a), num_uavs, 1), pos(:,1), pos(:,2), ...
        'VariableNames', {'UAV_Index', 'Altitude', 'X_Position', 'Y_Position'});
    altitude_table = [altitude_table; temp_table];
end
writetable(altitude_table, fullfile(output_dir, 'UAV_Positions_Altitude.xlsx'));

% 3. Save Percentile Positions
percentile_table = table();
for p = 1:length(percentiles)
    pos = positions_percentile(:,:,p);
    temp_table = table((1:num_uavs)', repmat(percentiles(p), num_uavs, 1), pos(:,1), pos(:,2), ...
        'VariableNames', {'UAV_Index', 'Percentile', 'X_Position', 'Y_Position'});
    percentile_table = [percentile_table; temp_table];
end
writetable(percentile_table, fullfile(output_dir, 'UAV_Positions_Percentile.xlsx'));

disp('UAV positions saved to:');
disp(fullfile(output_dir, 'UAV_Positions_Min_Distance.xlsx'));
disp(fullfile(output_dir, 'UAV_Positions_Altitude.xlsx'));
disp(fullfile(output_dir, 'UAV_Positions_Percentile.xlsx'));

% Plot sensitivity results (Coverage)
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
subplot(3,1,1);
plot(min_distances, coverages_min_distance, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Minimum Distance (m)', 'FontSize', 14);
ylabel('Coverage', 'FontSize', 14);
title('Sensitivity to Minimum Distance', 'FontSize', 16);
grid on;

subplot(3,1,2);
plot(altitudes, coverages_altitude, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Altitude (m)', 'FontSize', 14);
ylabel('Coverage', 'FontSize', 14);
title('Sensitivity to Altitude', 'FontSize', 16);
grid on;

subplot(3,1,3);
plot(percentiles, coverages_percentile, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Intensity Percentile', 'FontSize', 14);
ylabel('Coverage', 'FontSize', 14);
title('Sensitivity to Intensity Percentile', 'FontSize', 16);
grid on;

exportgraphics(gcf, 'Sensitivity_Analysis.png', 'Resolution', 300);

%% Visualize UAV 1 Position Trajectory for Altitude
% Initialize figure
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Plot demand heatmap
imagesc(x_range, y_range, Z);
set(gca, 'YDir', 'normal');
colormap('hot');
colorbar;
hold on;

% Extract UAV 1 positions for all altitudes
uav1_pos = squeeze(positions_altitude(1,:,:))'; % [5, 2] for x, y across 5 altitudes
altitudes = [100, 125, 150, 175, 200];
markers = {'*', 's', 'o', '^', 'd'}; % Star for 100 m, circle for 150 m, etc.

% Plot trajectory
for a = 1:length(altitudes)
    scatter(uav1_pos(a,1), uav1_pos(a,2), 50, 'k', 'Marker', markers{a}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Altitude = %d m', altitudes(a)));
    if altitudes(a) == 100
        text(uav1_pos(a,1) + 10, uav1_pos(a,2), '100 m', 'FontSize', 10, 'Color', 'k');
    end
end
plot(uav1_pos(:,1), uav1_pos(:,2), 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

% Customize plot
xlabel('X (meters)', 'FontSize', 14);
ylabel('Y (meters)', 'FontSize', 14);
title('UAV 1 Position Trajectory for Varying Altitude', 'FontSize', 16);
axis equal tight;
legend('Location', 'best');
%%
% Save figure
output_dir = 'C:\Users\mky6p\OneDrive - University of Missouri\Mill\Results\Result11\Sensitivity_Positions';
exportgraphics(gcf, fullfile(output_dir, 'UAV1_Position_Trajectory_Altitude.png'), 'Resolution', 300);
%% Visualize UAV Position Changes with Displacement Scatter Plot
% Calculate displacements relative to default positions
default_min_distance_idx = find(min_distances == 30);
default_altitude_idx = find(altitudes == 150);
default_percentile_idx = find(percentiles == 75);

% Extract default positions
ref_pos_min = positions_min_distance(:,:,default_min_distance_idx);
ref_pos_alt = positions_altitude(:,:,default_altitude_idx);
ref_pos_per = positions_percentile(:,:,default_percentile_idx);

% Verify reference positions are consistent
if ~isequal(ref_pos_min, ref_pos_alt, ref_pos_per)
    warning('Reference positions are not identical across parameters. Using min_distance reference.');
end
ref_pos = ref_pos_min; % Use min_distance as reference

% Initialize figure
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% 1. Min Distance
nexttile;
hold on;
colormap(parula(length(min_distances)));
for m = 1:length(min_distances)
    delta_pos = positions_min_distance(:,:,m) - ref_pos;
    scatter(delta_pos(:,1), delta_pos(:,2), 50, m*ones(num_uavs,1), 'filled', 'Marker', 'o', ...
        'DisplayName', sprintf('Min Distance = %d m', min_distances(m)));
end
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Reference (No Change)');
% Add reference circle (50-meter radius)
theta = linspace(0, 2*pi, 100);
plot(50*cos(theta), 50*sin(theta), 'k--', 'DisplayName', '50 m Radius');
xlabel('\DeltaX (meters)', 'FontSize', 14);
ylabel('\DeltaY (meters)', 'FontSize', 14);
title('Displacement (Min Distance)', 'FontSize', 16);
axis equal;
grid on;
colorbar('Ticks', 1:length(min_distances), 'TickLabels', min_distances);
xlim([-200, 200]);
ylim([-200, 200]);
legend('Location', 'best');

% 2. Altitude
nexttile;
hold on;
colormap(parula(length(altitudes)));
for a = 1:length(altitudes)
    delta_pos = positions_altitude(:,:,a) - ref_pos;
    scatter(delta_pos(:,1), delta_pos(:,2), 50, a*ones(num_uavs,1), 'filled', 'Marker', 'o', ...
        'DisplayName', sprintf('Altitude = %d m', altitudes(a)));
end
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Reference (No Change)');
plot(50*cos(theta), 50*sin(theta), 'k--', 'DisplayName', '50 m Radius');
xlabel('\DeltaX (meters)', 'FontSize', 14);
ylabel('\DeltaY (meters)', 'FontSize', 14);
title('Displacement (Altitude)', 'FontSize', 16);
axis equal;
grid on;
colorbar('Ticks', 1:length(altitudes), 'TickLabels', altitudes);
xlim([-200, 200]);
ylim([-200, 200]);
legend('Location', 'best');

% 3. Percentile
nexttile;
hold on;
colormap(parula(length(percentiles)));
for p = 1:length(percentiles)
    delta_pos = positions_percentile(:,:,p) - ref_pos;
    scatter(delta_pos(:,1), delta_pos(:,2), 50, p*ones(num_uavs,1), 'filled', 'Marker', 'o', ...
        'DisplayName', sprintf('Percentile = %d', percentiles(p)));
end
plot(0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Reference (No Change)');
plot(50*cos(theta), 50*sin(theta), 'k--', 'DisplayName', '50 m Radius');
xlabel('\DeltaX (meters)', 'FontSize', 14);
ylabel('\DeltaY (meters)', 'FontSize', 14);
title('Displacement (Percentile)', 'FontSize', 16);
axis equal;
grid on;
colorbar('Ticks', 1:length(percentiles), 'TickLabels', percentiles);
xlim([-200, 200]);
ylim([-200, 200]);
legend('Location', 'best');

% Save figure
exportgraphics(gcf, fullfile(output_dir, 'UAV_Position_Displacement_Scatter.png'), 'Resolution', 300);

%% Function Definitions
function [coverage, uav_positions] = run_optimization(X, Y, Z, hotspot_mask, min_distance, altitude, percentile, num_uavs, frequencies, uav_positions_initial, c)
    % Compute high-intensity threshold
    hotspot_intensities = Z(hotspot_mask);
    if ~isempty(hotspot_intensities)
        high_intensity_threshold = prctile(hotspot_intensities, percentile);
    else
        high_intensity_threshold = 0.5;
    end
    
    % Select candidate positions
    available_positions = [X(Z > high_intensity_threshold), Y(Z > high_intensity_threshold)];
    disp(['Number of available positions: ', num2str(size(available_positions, 1))]);
    
    % Start with initial positions
    uav_positions = uav_positions_initial;
    
    % Greedy placement to refine positions
    for i = 1:num_uavs
        best_coverage_local = -Inf;
        best_position_local = [];
        
        temp_coverage = zeros(size(available_positions, 1), 1);
        temp_positions = zeros(size(available_positions, 1), 2);
        
        parfor j = 1:size(available_positions, 1)
            x_candidate = available_positions(j, 1);
            y_candidate = available_positions(j, 2);
            
            too_close = false;
            for k = 1:i-1
                if norm([x_candidate, y_candidate] - uav_positions(k, :)) < min_distance
                    too_close = true;
                    break;
                end
            end
            
            if ~too_close
                temp_pos = [uav_positions(1:i-1, :); x_candidate, y_candidate];
                temp_coverage(j) = compute_total_coverage(X, Y, Z, temp_pos, frequencies(1:i), altitude, c);
                temp_positions(j, :) = [x_candidate, y_candidate];
            end
        end
        
        [best_coverage_local, best_idx] = max(temp_coverage);
        if best_coverage_local > -Inf
            best_position_local = temp_positions(best_idx, :);
        end
        
        if ~isempty(best_position_local)
            uav_positions(i, :) = best_position_local;
            available_positions(ismember(available_positions, best_position_local, 'rows'), :) = [];
        else
            disp(['No available position for UAV ', num2str(i), ', keeping initial position.']);
        end
    end
    
    % Compute coverage
    [coverage, ~] = compute_total_coverage(X, Y, Z, uav_positions, frequencies, altitude, c);
end

function [coverage, normalized_signal] = compute_total_coverage(X, Y, Z, uav_positions, frequencies, altitude, c)
    % Initialize signal strength grid (in dBm)
    num_uavs = size(uav_positions, 1);
    signal_strength = zeros(size(X, 1), size(X, 2), num_uavs);
    Pt = 20; % Transmit power in dBm (fixed for simplicity)

    % Calculate signal strength for each UAV
    for i = 1:num_uavs
        x_u = uav_positions(i, 1);
        y_u = uav_positions(i, 2);
        frequency = frequencies(i);
        
        % Compute distance from UAV to all grid points
        distance = sqrt((X - x_u).^2 + (Y - y_u).^2 + altitude^2);
        
        % Compute FSPL
        fspl_constant = 20 * log10(4 * pi / c) + 20 * log10(frequency);
        fspl = 20 * log10(distance) + fspl_constant;
        
        % Compute received signal strength in dBm
        signal_strength(:,:,i) = Pt - fspl;
    end
    
    % Find the strongest signal at each grid point
    [strongest_signal, ~] = max(signal_strength, [], 3);
    
    % Normalize the strongest signal (0 to 1)
    normalized_signal = (strongest_signal - min(strongest_signal(:))) / ...
                       (max(strongest_signal(:)) - min(strongest_signal(:)));
    
    % Compute coverage as the sum of demand intensity weighted by normalized signal
    coverage = sum(Z .* normalized_signal, 'all');
    
    % Add penalty for UAVs in low-intensity regions
    penalty_weight = 1;
    for i = 1:num_uavs
        [~, x_idx] = min(abs(X(1, :) - uav_positions(i, 1)));
        [~, y_idx] = min(abs(Y(:, 1) - uav_positions(i, 2)));
        intensity_at_uav = Z(y_idx, x_idx);
        
        if intensity_at_uav < 0.1 % Adjust threshold as needed
            coverage = coverage - penalty_weight * (0.1 - intensity_at_uav);
        end
    end
end