%% UAV Optimization with Greedy Algorithm (Improved with Gradient-Based Threshold)
clc;
clear;

% Load building and UAV data
buildings_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Buildings data\Building_joplin-tornado_00000008.xlsx";
uav_file = "C:\Users\mky6p\OneDrive - University of Missouri\Mill\Results\Result8\UAV_Inverse_FSPL_Random_Positions_and_Frequencies.xlsx";

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
figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); 

markerSize = 80;
fontSize = 20;
xLimits = [min(x_range), max(x_range)];
yLimits = [min(y_range), max(y_range)];

% ---- Create subplots with a 2x3 layout to accommodate new plots ----
t = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

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
legend('UAV Positions', 'FontSize', fontSize); % Add legend for UAV positions
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
legend('UAV Positions', 'FontSize', fontSize); % Add legend for UAV positions
hold off;

% 4. Gradient Magnitude Heatmap with Threshold Overlay
nexttile;
imagesc(x_range, y_range, grad_magnitude);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'jet');
colorbar;
hold on;
contour(x_range, y_range, grad_magnitude, [grad_threshold grad_threshold], 'k-', 'LineWidth', 2);
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title(['Gradient Magnitude (Threshold = ', num2str(grad_threshold), ')'], 'FontSize', fontSize);

% 5. Hotspot Mask
nexttile;
imagesc(x_range, y_range, hotspot_mask);
axis equal tight;
set(gca, 'YDir', 'normal');
colormap(gca, 'gray');
colorbar;
xlabel('X (meters)', 'FontSize', fontSize);
ylabel('Y (meters)', 'FontSize', fontSize);
title('Hotspot Mask', 'FontSize', fontSize);

% 6. Histogram of Hotspot Intensities with Intensity Threshold
nexttile;
histogram(hotspot_intensities, 50, 'FaceColor', 'b', 'EdgeColor', 'k');
hold on;
xline(high_intensity_threshold, 'r-', 'LineWidth', 2, 'Label', ['Threshold = ', num2str(high_intensity_threshold)]);
xlabel('Intensity', 'FontSize', fontSize);
ylabel('Frequency', 'FontSize', fontSize);
title('Hotspot Intensities Distribution', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);

% ---- Ensure no legend appears in heatmaps ----
legend off

% ---- Export high-resolution, clean figure ----
exportgraphics(gcf, 'Result_Final_With_Threshold_Visualizations.png', 'Resolution', 300);

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