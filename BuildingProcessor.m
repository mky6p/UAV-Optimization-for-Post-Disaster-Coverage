% MATLAB code equivalent to the provided Python script with updated headers

% Load JSON file (assumes file is in current directory or specify full path)
json_path = "C:\Users\mky6p\OneDrive - University of Missouri\AIAA Aviation 2025\Labels\joplin-tornado_00000001_pre_disaster.json"; % Update path as needed
json_data = jsondecode(fileread(json_path));

% Step 1: Extract building data from JSON
building_data = struct();

% Process each feature
features = json_data.features.lng_lat;
for k = 1:length(features)
    properties = features(k).properties;
    uid = properties.uid;
    
    % Check if 'subtype' exists; default for pre-disaster
    if isfield(properties, 'subtype')
        damage_level = properties.subtype;
    else
        damage_level = 'no-damage'; % Or 'un-classified' if preferred for pre-disaster
    end
    if isempty(damage_level)
        damage_level = 'no-damage';
    end
    
    % Parse WKT (simple polygon parsing, assumes POLYGON format)
    wkt_str = features(k).wkt;
    coords_str = extractBetween(wkt_str, '((', '))');
    coords_cell = split(coords_str{1}, ',');
    coords = zeros(length(coords_cell), 2);
    for i = 1:length(coords_cell)
        xy = split(strtrim(coords_cell{i}));
        coords(i, :) = [str2double(xy{2}), str2double(xy{1})]; % [lat, lon]
    end
    coords = coords(1:end-1, :); % Remove duplicate last point
    
    % Calculate area
    real_area = calculate_real_area(coords);
    
    % Calculate centroid (simple average for this example)
    centroid_lat = mean(coords(:, 1));
    centroid_lon = mean(coords(:, 2));
    
    % Store data
    building_data(k).UID = uid;
    building_data(k).DamageLevel = damage_level;
    building_data(k).RealArea_sqm = real_area;
    building_data(k).Latitude = centroid_lat;
    building_data(k).Longitude = centroid_lon;
end

% Convert struct to table (MATLAB equivalent of DataFrame)
df = struct2table(building_data);

% Calculate size categories
q1 = quantile(df.RealArea_sqm, 0.33);
q2 = quantile(df.RealArea_sqm, 0.66);

df.SizeCategory = arrayfun(@(x) classify_relative_size(x, q1, q2), df.RealArea_sqm, 'UniformOutput', false);

% Map damage levels to intensity
subtype_to_intensity = containers.Map(...
    {'destroyed', 'major-damage', 'minor-damage', 'no-damage', 'un-classified'}, ...
    [1.0, 0.8, 0.5, 0.05, 0.01]);
df.Intensity = cellfun(@(x) subtype_to_intensity(x), df.DamageLevel);

% Calculate sigma and divide by 100,000 (as per previous request)
df.Sigma = sqrt(df.RealArea_sqm) / 100000;

% Save initial table to Excel
initial_file_path = 'updated_building_data_new_json.xlsx';
writetable(df, initial_file_path);
fprintf('Initial Excel file saved at: %s\n', initial_file_path);

% Step 2: Extract latitude and longitude bounds
min_lat = min(df.Latitude);
max_lat = max(df.Latitude);
min_lon = min(df.Longitude);
max_lon = max(df.Longitude);

% Calculate reference point
lat_ref = (min_lat + max_lat) / 2;
lon_ref = (min_lon + max_lon) / 2;

% Conversion factors
m_per_deg_lat = 111320;
m_per_deg_lon = 111320 * cosd(lat_ref);

% Define grid for the heatmap
x = linspace(min_lon, max_lon, 300);
y = linspace(min_lat, max_lat, 300);
x_m = (x - lon_ref) * m_per_deg_lon;
y_m = (y - lat_ref) * m_per_deg_lat;
[X, Y] = meshgrid(x_m, y_m);
Z = zeros(size(X));

% Define target dimensions for pixel range
target_width = 600;
target_height = 600;
x_min = min(x_m);
x_max = max(x_m);
y_min = min(y_m);
y_max = max(y_m);
x_scale = target_width / (x_max - x_min);
y_scale = target_height / (y_max - y_min);

% Calculate Gaussian and update Z for each building
output_data = struct();
for k = 1:height(df)
    x0_m = (df.Longitude(k) - lon_ref) * m_per_deg_lon;
    y0_m = (df.Latitude(k) - lat_ref) * m_per_deg_lat;
    intensity = df.Intensity(k);
    sigma_m = df.Sigma(k) * m_per_deg_lat; % Sigma already divided by 100,000
    damage_level = df.DamageLevel{k};
    
    x0_pixel = (x0_m - x_min) * x_scale;
    y0_pixel = (y0_m - y_min) * y_scale;
    
    Z_building = gaussian_kernel(X, Y, x0_m, y0_m, sigma_m, intensity);
    Z = Z + Z_building;
    
    formula = sprintf('%f * exp(-((x - %f)^2 + (y - %f)^2) / (2 * %f^2))', ...
        intensity, x0_m, y0_m, sigma_m);
    total_fx_y = sum(Z_building, 'all');
    
    output_data(k).X_meters = x0_m;            % Matches "X (meters)"
    output_data(k).Y_meters = y0_m;            % Matches "Y (meters)"
    output_data(k).X_pixels = x0_pixel;        % Matches "X (pixels)"
    output_data(k).Y_pixels = y0_pixel;        % Matches "Y (pixels)"
    output_data(k).Intensity = intensity;      % Matches "Intensity"
    output_data(k).Sigma_meters = sigma_m;        % Matches "ma (meters)" (intended as Sigma)
    output_data(k).DamageLevel = damage_level; % Matches "Damage Level"
    output_data(k).Formula_f_x_y = formula;    % Matches "Formula f(x,y)"
    output_data(k).Total_f_x_y = total_fx_y;   % Matches "Total f(x,y)"
end

% Convert output struct to table and save to Excel with custom headers
output_df = struct2table(output_data);
% Rename the columns to match the desired headers exactly
output_df.Properties.VariableNames = {
    'X (meters)', ...
    'Y (meters)', ...
    'X (pixels)', ...
    'Y (pixels)', ...
    'Intensity', ...
    'Sigma (meters)', ...
    'Damage Level', ...
    'Formula f(x,y)', ...
    'Total f(x,y)'};

output_file = "C:\Users\mky6p\OneDrive - University of Missouri\AIAA Aviation 2025\Labels\Building_joplin-tornado_00000001.xlsx";
writetable(output_df, output_file);
fprintf('Output saved to: %s\n', output_file);

% Optional: Display first few rows (equivalent to head())
disp(head(output_df, 5));

% --- Function Definitions ---

function area = calculate_real_area(coords)
    if size(coords, 1) < 3
        area = 0.0;
        return;
    end
    lat_ref = mean(coords(:, 1)); % Approximate latitude for conversion
    meters_per_lat = 111320;
    meters_per_lon = 111320 * cosd(lat_ref); % cosd uses degrees
    area = 0.0;
    for i = 1:size(coords, 1)
        j = mod(i, size(coords, 1)) + 1;
        x1 = coords(i, 2) * meters_per_lon;
        y1 = coords(i, 1) * meters_per_lat;
        x2 = coords(j, 2) * meters_per_lon;
        y2 = coords(j, 1) * meters_per_lat;
        area = area + (x1 * y2) - (x2 * y1);
    end
    area = abs(area) / 2; % Area in square meters
end

function size_cat = classify_relative_size(area, q1, q2)
    if area <= q1
        size_cat = 'small';
    elseif area <= q2
        size_cat = 'medium';
    else
        size_cat = 'large';
    end
end

function z = gaussian_kernel(x, y, x0, y0, sigma, intensity)
    z = intensity * exp(-((x - x0).^2 + (y - y0).^2) / (2 * sigma^2));
end