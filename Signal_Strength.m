clc
clear all
close all

% Constants
c = 3e8;  % Speed of light in m/s
altitude = 150;  % UAV altitude in meters

% Define grid dimensions in meters
x_range = [-300, 500];  % X-axis range in meters
y_range = [-300, 500];  % Y-axis range in meters

% Define grid
grid_points = 500;
x = linspace(x_range(1), x_range(2), grid_points);
y = linspace(y_range(1), y_range(2), grid_points);
[X, Y] = meshgrid(x, y);

% User input for number of UAVs
num_uavs = input('Enter the number of UAVs: ');
while num_uavs < 1
    num_uavs = input('Please enter a positive number of UAVs: ');
end

% Assign frequencies linearly from 2 GHz to 3 GHz
freq_start = 2e9; % 2 GHz
freq_end = 3e9;   % 3 GHz
uav_frequencies = linspace(freq_start, freq_end, num_uavs); % Unique frequencies

% Generate random positions for UAVs within the specified range
uav_positions = [ 
    x_range(1) + (x_range(2) - x_range(1)) * rand(1, num_uavs); 
    y_range(1) + (y_range(2) - y_range(1)) * rand(1, num_uavs);
    altitude * ones(1, num_uavs); % Fixed altitude
    20 * ones(1, num_uavs)        % Transmit power (20 dBm)
];

% Initialize signal strength matrix (in dBm) for each UAV
signal_strength = zeros(grid_points, grid_points, num_uavs);

% Calculate FSPL and signal strength for each UAV
for k = 1:num_uavs
    % UAV position and properties
    uav_x = uav_positions(1, k);
    uav_y = uav_positions(2, k);
    uav_z = uav_positions(3, k);
    Pt = uav_positions(4, k); % Transmit power in dBm
    freq = uav_frequencies(k); % Unique frequency

    % Distance from UAV to each grid point (3D Euclidean distance)
    d = sqrt((X - uav_x).^2 + (Y - uav_y).^2 + uav_z.^2);

    % Free Space Path Loss (FSPL) in dB
    FSPL = 20 * log10(d) + 20 * log10(freq) + 20 * log10(4 * pi / c);

    % Received power in dBm
    signal_strength(:,:,k) = Pt - FSPL;
end

% Determine which UAV has the highest signal strength at each grid point
[strongest_signal, strongest_uav] = max(signal_strength, [], 3);

% Normalize the strongest signal strength (0 to 1)
normalized_signal = (strongest_signal - min(strongest_signal(:))) / ...
                   (max(strongest_signal(:)) - min(strongest_signal(:)));

% === Plot normalized combined heat map and dominant UAV map next to each other ===
figure('Position', [100, 100, 1200, 500]);

% --- Subplot 1: Normalized Strongest Signal Strength ---
subplot(1,2,1); % 1 row, 2 columns, 1st plot
imagesc(x_range, y_range, normalized_signal);
colormap('parula');
colorbar;
title('Normalized Strongest Signal Strength');
xlabel('X (m)');
ylabel('Y (m)');
axis xy equal tight;
caxis([0, 1]); % Normalized scale
hold on;
scatter(uav_positions(1,:), uav_positions(2,:), 36, 'ko', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
for k = 1:num_uavs
    text(uav_positions(1,k) + 10, uav_positions(2,k), ['UAV ', num2str(k)], 'Color', 'white');
end
hold off;

% --- Subplot 2: Dominant UAV per Grid Point ---
subplot(1,2,2); % 1 row, 2 columns, 2nd plot
imagesc(x_range, y_range, strongest_uav);
colormap('parula');
colorbar('Ticks', 1:num_uavs, 'TickLabels', arrayfun(@(x) ['UAV ', num2str(x)], 1:num_uavs, 'UniformOutput', false));
title('Dominant UAV per Grid Point');
xlabel('X (m)');
ylabel('Y (m)');
axis xy equal tight;
hold on;
scatter(uav_positions(1,:), uav_positions(2,:), 36, 'ko', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
hold off;
