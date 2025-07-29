%% UAV Optimization with Greedy Algorithm (Comparing Thresholding Methods)
clc;
clear;

% Start timing for the entire script
tic;

% Load building and UAV data
buildings_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Buildings data\Building_joplin-tornado_00000011.xlsx";
uav_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Results\Result11\UAV_Inverse_FSPL_Random_Positions_and_Frequencies.xlsx";

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

% Extract UAV positions and frequencies from the file
uav_positions_initial = [uav_data.X_Position_m, uav_data.Y_Position_m]; % Initial positions from file
frequencies = uav_data.Frequency_Hz;
num_uavs = length(frequencies);

% Ensure the number of UAVs matches the file data
if size(uav_positions_initial, 1) ~= num_uavs
    error('Number of UAVs in position data does not match the number of frequencies.');
end

% --- Thresholding Methods ---

% 1. No Thresholding: Use all grid points as candidate positions
all_positions = [X(:), Y(:)];

% 2. Static Thresholding: Use a fixed intensity threshold
static_threshold = 0.5; % Fixed threshold for intensity
static_positions = [X(Z > static_threshold), Y(Z > static_threshold)];
disp(['Number of Candidate Positions (Static Thresholding): ', num2str(size(static_positions, 1))]);

% 3. Gradient-Based Thresholding
[grad_x, grad_y] = gradient(Z, x_range, y_range);
grad_magnitude = sqrt(grad_x.^2 + grad_y.^2);
grad_magnitude = grad_magnitude / max(grad_magnitude(:));
grad_threshold = 0.3; % Gradient threshold
hotspot_mask = grad_magnitude > grad_threshold;
hotspot_intensities = Z(hotspot_mask);
if ~isempty(hotspot_intensities)
    high_intensity_threshold = prctile(hotspot_intensities, 75); % 75th percentile
else
    high_intensity_threshold = 0.5; % Fallback
end
disp(['Computed High Intensity Threshold (Gradient-Based): ', num2str(high_intensity_threshold)]);
gradient_positions = [X(Z > high_intensity_threshold), Y(Z > high_intensity_threshold)];
disp(['Number of Candidate Positions (Gradient-Based Thresholding): ', num2str(size(gradient_positions, 1))]);

% --- Greedy Placement for Each Thresholding Method ---

% Initialize UAV positions for each method
uav_positions_no_threshold = uav_positions_initial;
uav_positions_static = uav_positions_initial;
uav_positions_gradient = uav_positions_initial;

% 1. No Thresholding
disp('Starting Greedy Placement with No Thresholding (Parallel)...');
available_positions = all_positions;
tic;
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
            if norm([x_candidate, y_candidate] - uav_positions_no_threshold(k, :)) < min_distance
                too_close = true;
                break;
            end
        end
        
        if ~too_close
            temp_pos = [uav_positions_no_threshold(1:i-1, :); x_candidate, y_candidate];
            temp_coverage(j) = compute_total_coverage(X, Y, Z, temp_pos, frequencies(1:i), altitude, c);
            temp_positions(j, :) = [x_candidate, y_candidate];
        end
    end
    
    [best_coverage_local, best_idx] = max(temp_coverage);
    if best_coverage_local > -Inf
        best_position_local = temp_positions(best_idx, :);
    end
    
    if ~isempty(best_position_local)
        uav_positions_no_threshold(i, :) = best_position_local;
    end
    
    available_positions(ismember(available_positions, best_position_local, 'rows'), :) = [];
end
time_no_threshold = toc;
disp(['Execution Time (No Thresholding, Parallel): ', num2str(time_no_threshold), ' seconds']);

% 2. Static Thresholding
disp('Starting Greedy Placement with Static Thresholding (Parallel)...');
available_positions = static_positions;
tic;
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
            if norm([x_candidate, y_candidate] - uav_positions_static(k, :)) < min_distance
                too_close = true;
                break;
            end
        end
        
        if ~too_close
            temp_pos = [uav_positions_static(1:i-1, :); x_candidate, y_candidate];
            temp_coverage(j) = compute_total_coverage(X, Y, Z, temp_pos, frequencies(1:i), altitude, c);
            temp_positions(j, :) = [x_candidate, y_candidate];
        end
    end
    
    [best_coverage_local, best_idx] = max(temp_coverage);
    if best_coverage_local > -Inf
        best_position_local = temp_positions(best_idx, :);
    end
    
    if ~isempty(best_position_local)
        uav_positions_static(i, :) = best_position_local;
    end
    
    available_positions(ismember(available_positions, best_position_local, 'rows'), :) = [];
end
time_static = toc;
disp(['Execution Time (Static Thresholding, Parallel): ', num2str(time_static), ' seconds']);

% 3. Gradient-Based Thresholding
disp('Starting Greedy Placement with Gradient-Based Thresholding (Parallel)...');
available_positions = gradient_positions;
tic;
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
            if norm([x_candidate, y_candidate] - uav_positions_gradient(k, :)) < min_distance
                too_close = true;
                break;
            end
        end
        
        if ~too_close
            temp_pos = [uav_positions_gradient(1:i-1, :); x_candidate, y_candidate];
            temp_coverage(j) = compute_total_coverage(X, Y, Z, temp_pos, frequencies(1:i), altitude, c);
            temp_positions(j, :) = [x_candidate, y_candidate];
        end
    end
    
    [best_coverage_local, best_idx] = max(temp_coverage);
    if best_coverage_local > -Inf
        best_position_local = temp_positions(best_idx, :);
    end
    
    if ~isempty(best_position_local)
        uav_positions_gradient(i, :) = best_position_local;
    end
    
    available_positions(ismember(available_positions, best_position_local, 'rows'), :) = [];
end
time_gradient = toc;
disp(['Execution Time (Gradient-Based Thresholding, Parallel): ', num2str(time_gradient), ' seconds']);

% Display UAV positions for each method
disp('Greedy Algorithm Optimal UAV Positions (No Thresholding):');
disp(uav_positions_no_threshold);
disp('Greedy Algorithm Optimal UAV Positions (Static Thresholding):');
disp(uav_positions_static);
disp('Greedy Algorithm Optimal UAV Positions (Gradient-Based Thresholding):');
disp(uav_positions_gradient);

% --- Visualization of UAV Placements ---
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
markerSize = 80;
fontSize = 20;
xLimits = [min(x_range), max(x_range)];
yLimits = [min(y_range), max(y_range)];

t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% 1. No Thresholding
nexttile;
imagesc(x_range, y_range, Z);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('No Thresholding', 'FontSize', fontSize);
hold on;
scatter(uav_positions_no_threshold(:, 1), uav_positions_no_threshold(:, 2), markerSize, 'bx', 'LineWidth', 2);
hold off;

% 2. Static Thresholding
nexttile;
imagesc(x_range, y_range, Z);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Static Thresholding', 'FontSize', fontSize);
hold on;
scatter(uav_positions_static(:, 1), uav_positions_static(:, 2), markerSize, 'bx', 'LineWidth', 2);
hold off;

% 3. Gradient-Based Thresholding
nexttile;
imagesc(x_range, y_range, Z);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Gradient-Based Thresholding', 'FontSize', fontSize);
hold on;
scatter(uav_positions_gradient(:, 1), uav_positions_gradient(:, 2), markerSize, 'bx', 'LineWidth', 2);
hold off;

% Export the figure
exportgraphics(gcf, 'uav_placements_thresholding_comparison.png', 'Resolution', 300);

% --- Compute Performance Metrics for Each Method ---
% No Thresholding
[coverage_no_threshold, normalized_signal_no_threshold] = compute_total_coverage(X, Y, Z, uav_positions_no_threshold, frequencies, altitude, c);
aae_no_threshold = (1/numel(Z)) * sum(abs(Z(:) - normalized_signal_no_threshold(:)));
frobenius_no_threshold = sqrt(sum((Z(:) - normalized_signal_no_threshold(:)).^2));
peak_error_no_threshold = abs(max(Z(:)) - max(normalized_signal_no_threshold(:)));
coverage_ratio_no_threshold = sum(normalized_signal_no_threshold(:) > 0.5) / numel(Z);

% Static Thresholding
[coverage_static, normalized_signal_static] = compute_total_coverage(X, Y, Z, uav_positions_static, frequencies, altitude, c);
aae_static = (1/numel(Z)) * sum(abs(Z(:) - normalized_signal_static(:)));
frobenius_static = sqrt(sum((Z(:) - normalized_signal_static(:)).^2));
peak_error_static = abs(max(Z(:)) - max(normalized_signal_static(:)));
coverage_ratio_static = sum(normalized_signal_static(:) > 0.5) / numel(Z);

% Gradient-Based Thresholding
[coverage_gradient, normalized_signal_gradient] = compute_total_coverage(X, Y, Z, uav_positions_gradient, frequencies, altitude, c);
aae_gradient = (1/numel(Z)) * sum(abs(Z(:) - normalized_signal_gradient(:)));
frobenius_gradient = sqrt(sum((Z(:) - normalized_signal_gradient(:)).^2));
peak_error_gradient = abs(max(Z(:)) - max(normalized_signal_gradient(:)));
coverage_ratio_gradient = sum(normalized_signal_gradient(:) > 0.5) / numel(Z);

% Display Performance Metrics
disp('--- Performance Metrics ---');
disp('No Thresholding:');
disp(['Weighted Signal Strength: ', num2str(-coverage_no_threshold)]);
disp(['AAE: ', num2str(aae_no_threshold)]);
disp(['Frobenius Norm Error: ', num2str(frobenius_no_threshold)]);
disp(['Peak Error: ', num2str(peak_error_no_threshold)]);
disp(['Signal Coverage Ratio: ', num2str(coverage_ratio_no_threshold)]);

disp('Static Thresholding:');
disp(['Weighted Signal Strength: ', num2str(-coverage_static)]);
disp(['AAE: ', num2str(aae_static)]);
disp(['Frobenius Norm Error: ', num2str(frobenius_static)]);
disp(['Peak Error: ', num2str(peak_error_static)]);
disp(['Signal Coverage Ratio: ', num2str(coverage_ratio_static)]);

disp('Gradient-Based Thresholding:');
disp(['Weighted Signal Strength: ', num2str(-coverage_gradient)]);
disp(['AAE: ', num2str(aae_gradient)]);
disp(['Frobenius Norm Error: ', num2str(frobenius_gradient)]);
disp(['Peak Error: ', num2str(peak_error_gradient)]);
disp(['Signal Coverage Ratio: ', num2str(coverage_ratio_gradient)]);

%% Report Total Execution Time and Computer Configuration
total_time = toc;
disp('--- Execution Time Summary ---');
disp(['Total Script Execution Time: ', num2str(total_time), ' seconds']);
disp(['Greedy Placement (No Thresholding, Parallel): ', num2str(time_no_threshold), ' seconds']);
disp(['Greedy Placement (Static Thresholding, Parallel): ', num2str(time_static), ' seconds']);
disp(['Greedy Placement (Gradient-Based Thresholding, Parallel): ', num2str(time_gradient), ' seconds']);

disp('--- Computer Configuration ---');
[~, sys_info] = system('wmic computersystem get manufacturer, model, totalphysicalmemory /format:list');
disp('System Information:');
disp(sys_info);
disp(['MATLAB Version: ', version]);
disp(['Number of Workers (Parallel Pool): ', num2str(parallel.defaultClusterProfile)]);


%% Function Definitions
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