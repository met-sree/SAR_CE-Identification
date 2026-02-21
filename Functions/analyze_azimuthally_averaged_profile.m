function [condition_result, eyewall_data, moat_data] = analyze_azimuthally_averaged_profile(wsp_sar_polar, radial_dist_km)
% Analyze azimuthally averaged wind profile
% Input:
%   wsp_sar_polar: SAR polar wind speed data in KNOTS (angle × radial distance)
%   radial_dist_km: radial distance array (km)
% Output:
%   condition_result: 1 if concentric eyewall detected, 0 otherwise
%   eyewall_data: [primary_dist, primary_wind, secondary_dist, secondary_wind] in knots
%   moat_data: [moat_dist, moat_wind] in knots (NaN if no moat)

% Minimum wind speed threshold in knots
MIN_WIND_KNOTS = 17;

% Ensure radial distance length matches data
if size(wsp_sar_polar, 2) ~= length(radial_dist_km)
    error('Radial distance dimension mismatch! wsp_sar_polar size: %d×%d, radial_dist_km length: %d', ...
        size(wsp_sar_polar, 1), size(wsp_sar_polar, 2), length(radial_dist_km));
end

% ========== 1. Calculate azimuthally averaged profile ==========
fprintf('Calculating azimuthally averaged profile...\n');

% Calculate azimuthally averaged profile (average across all angles)
azm_averaged = nanmean(wsp_sar_polar, 1)';  % Mean across all angles
azm_averaged_std = std(wsp_sar_polar, 0, 1, 'omitnan')';  % Standard deviation

fprintf('Profile length: %d points, Radial range: %.1f to %.1f km\n', ...
    length(radial_dist_km), min(radial_dist_km), max(radial_dist_km));

% ========== 2. Detect features ==========
fprintf('\nDetecting eyewall and moat features in azimuthally averaged profile...\n');
fprintf('Filtering criterion: Eyewalls AND moats must have wind speed >= %.1f knots\n', MIN_WIND_KNOTS);

% Create figure window
%figure('WindowState', 'maximized', 'Name', 'Azimuthally Averaged Wind Profile Analysis');
fig = figure('Position', [100, 100, 1400, 800]);
% Detection parameters
min_prominence = 2;  % Minimum prominence (in knots)
min_separation = 2;   % Minimum separation (km)

% Initialize outputs
condition_result = 0;
eyewall_data = [NaN, NaN, NaN, NaN];  % [primary_dist, primary_wind, secondary_dist, secondary_wind]
moat_data = [NaN, NaN];  % [moat_dist, moat_wind]

% Select data
aa = azm_averaged;      % Mean wind speed in knots
bb = azm_averaged_std;  % Standard deviation in knots

% Find maximum wind speed location
[max_wind, max_idx] = max(aa);
max_dist = radial_dist_km(max_idx);
fprintf('Maximum wind: %.1f knots at %.1f km\n', max_wind, max_dist);

% Calculate eyewall positions (local maxima)
isEyewall = islocalmax(aa, ...
    'MinProminence', min_prominence, ...
    'MinSeparation', min_separation, ...
    'SamplePoints', radial_dist_km);
eyewall_indices = find(isEyewall);
eyewall_values = aa(isEyewall);
eyewall_distances = radial_dist_km(eyewall_indices);

% Calculate moat positions (local minima)
isMoat = islocalmin(aa, ...
    'MinProminence', min_prominence, ...
    'MinSeparation', min_separation, ...
    'SamplePoints', radial_dist_km);
moat_indices = find(isMoat);
moat_values = aa(isMoat);
moat_distances = radial_dist_km(moat_indices);

% ========== FILTER EYEWALLS: Remove those with wind speed < 17 knots ==========
valid_eyewall_mask = eyewall_values >= MIN_WIND_KNOTS;
eyewall_indices = eyewall_indices(valid_eyewall_mask);
eyewall_values = eyewall_values(valid_eyewall_mask);
eyewall_distances = eyewall_distances(valid_eyewall_mask);

% ========== FILTER MOATS: Remove those with wind speed < 17 knots ==========
valid_moat_mask = moat_values >= MIN_WIND_KNOTS;
moat_indices = moat_indices(valid_moat_mask);
moat_values = moat_values(valid_moat_mask);
moat_distances = moat_distances(valid_moat_mask);

% Display detected features
fprintf('\nDetected %d eyewall(s) and %d moat(s) (after filtering)\n', ...
    length(eyewall_values), length(moat_values));
if ~isempty(eyewall_values)
    fprintf('Eyewalls at: ');
    for i = 1:length(eyewall_values)
        fprintf('%.1f km (%.1f knots) ', eyewall_distances(i), eyewall_values(i));
    end
    fprintf('\n');
end
if ~isempty(moat_values)
    fprintf('Moats at: ');
    for i = 1:length(moat_values)
        fprintf('%.1f km (%.1f knots) ', moat_distances(i), moat_values(i));
    end
    fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% for finding is it concentric eyewall or not %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables
primary_info = []; secondary_info = []; % Store [distance, windspeed] in knots
moat_info = []; % Store [distance, windspeed] in knots
primary_std = NaN; secondary_std = NaN; moat_wind = NaN;

% Check if we have at least 2 eyewalls AND at least 1 moat (after filtering)
if length(eyewall_distances) >= 2 && ~isempty(moat_distances)
    
    % Sort eyewalls by distance (ascending order)
    [sorted_eyewalls, sort_idx] = sort(eyewall_distances);
    sorted_eyewall_winds = eyewall_values(sort_idx);
    
    % Check all consecutive eyewall pairs
    for i = 1:(length(sorted_eyewalls) - 1)
        inner_eyewall_dist = sorted_eyewalls(i);
        outer_eyewall_dist = sorted_eyewalls(i + 1);
        inner_eyewall_wind = sorted_eyewall_winds(i);
        outer_eyewall_wind = sorted_eyewall_winds(i + 1);
        
        % Find moats between this specific pair
        moat_idx = find(moat_distances > inner_eyewall_dist & ...
            moat_distances < outer_eyewall_dist);
        
        if ~isempty(moat_idx)
            % Found moat between this eyewall pair
            
            % Store eyewall information for this pair
            primary_info = [inner_eyewall_dist, inner_eyewall_wind];
            secondary_info = [outer_eyewall_dist, outer_eyewall_wind];
            
            % Store all moats between this pair (take the first moat)
            moat_info = [moat_distances(moat_idx(1)), moat_values(moat_idx(1))];
            
            % Get the standard deviation at eyewall positions
            if length(eyewall_indices) >= 2
                primary_std = bb(eyewall_indices(sort_idx(i)));
                secondary_std = bb(eyewall_indices(sort_idx(i+1)));
            end
            
            % Get moat wind speed
            moat_wind = moat_info(1, 2);
            
            % Check the condition: (eyewall wind - eyewall std) > moat wind
            condition1 = (inner_eyewall_wind - primary_std) > 0.5*moat_wind;
            condition2 = (outer_eyewall_wind - secondary_std) > 0.5*moat_wind;
            
            if condition1 && condition2
                condition_result = 1;
                
                % Store data
                eyewall_data = [primary_info, secondary_info];
                moat_data = moat_info;
            else
                condition_result = 0;
                % Set to NaN since conditions not met
                eyewall_data = [NaN, NaN, NaN, NaN];
                moat_data = [NaN, NaN];
            end
            
            break; % Stop at first pair with moat
        end
    end
end

% If no moat found between any eyewall pair, set all to NaN
if condition_result == 0 && isempty(moat_info)
    eyewall_data = [NaN, NaN, NaN, NaN];
    moat_data = [NaN, NaN];
end

% Display results
if condition_result == 1
    fprintf('\n=== CONCENTRIC EYEWALL DETECTED in azimuthally averaged profile ===\n');
    fprintf('Primary Eyewall: Distance = %.2f km, Wind Speed = %.2f knots, Std = %.2f knots\n', ...
        primary_info(1), primary_info(2), primary_std);
    fprintf('Secondary Eyewall: Distance = %.2f km, Wind Speed = %.2f knots, Std = %.2f knots\n', ...
        secondary_info(1), secondary_info(2), secondary_std);
    fprintf('Moat: Distance = %.2f km, Wind Speed = %.2f knots\n', ...
        moat_info(1), moat_info(2));
    fprintf('Condition check: (%.1f-%.1f=%.1f) > %.1f AND (%.1f-%.1f=%.1f) > %.1f = TRUE\n', ...
        inner_eyewall_wind, primary_std, inner_eyewall_wind-primary_std, moat_wind, ...
        outer_eyewall_wind, secondary_std, outer_eyewall_wind-secondary_std, moat_wind);
else
    if length(eyewall_distances) >= 2 && ~isempty(moat_distances)
        fprintf('\n=== No concentric eyewall in azimuthally averaged profile (condition not met) ===\n');
    else
        fprintf('\n=== No concentric eyewall in azimuthally averaged profile: ');
        if length(eyewall_distances) < 2
            fprintf('Not enough eyewalls (need ≥2, found %d)', length(eyewall_distances));
        else
            fprintf('No valid moats (all filtered out or none found)');
        end
        fprintf(' ===\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ========== 3. Create plot ==========
% Plot base curve (with error bars) - all in knots
plot(radial_dist_km, aa, 'b-', 'LineWidth', 2);
hold on;

% Add standard deviation shading
upper_bound = aa + bb;
lower_bound = aa - bb;
fill([radial_dist_km, fliplr(radial_dist_km)], ...
    [upper_bound', fliplr(lower_bound')], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'DisplayName', '±1 Std Dev');

% ========== SYMBOL PLOTTING ==========
alpha = 0.5;

if condition_result == 1
    % Primary eyewall: MAGENTA asterisk
    scatter(primary_info(1), primary_info(2), 100, 'm*', ...
        'MarkerEdgeAlpha', alpha, ...
        'MarkerFaceAlpha', alpha, ...
        'LineWidth', 2, ...
        'DisplayName', 'Primary Eyewall');
    
    % Secondary eyewall: RED asterisk
    scatter(secondary_info(1), secondary_info(2), 100, 'r*', ...
        'MarkerEdgeAlpha', alpha, ...
        'MarkerFaceAlpha', alpha, ...
        'LineWidth', 2, ...
        'DisplayName', 'Secondary Eyewall');
    
    % Moat: BLUE asterisk
    scatter(moat_info(1), moat_info(2), 100, 'b*', ...
        'MarkerEdgeAlpha', alpha, ...
        'MarkerFaceAlpha', alpha, ...
        'LineWidth', 2, ...
        'DisplayName', 'Moat');
else
    % No concentric eyewall detected - Show maximum wind speed
    scatter(max_dist, max_wind, 100, 'm*', ...
        'MarkerEdgeAlpha', alpha, ...
        'MarkerFaceAlpha', alpha, ...
        'LineWidth', 2, ...
        'DisplayName', 'Maximum Wind');
end

% Plot all detected eyewalls and moats (lightly, for reference)
if ~isempty(eyewall_distances)
    scatter(eyewall_distances, eyewall_values, 50, 'k^', ...
        'MarkerFaceColor', [0.5 0.5 0.5], ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.3, ...
        'LineWidth', 1, ...
        'DisplayName', 'All Eyewalls');
end

if ~isempty(moat_distances)
    scatter(moat_distances, moat_values, 50, 'kv', ...
        'MarkerFaceColor', [0.5 0.5 0.5], ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', 0.3, ...
        'LineWidth', 1, ...
        'DisplayName', 'All Moats');
end

% Set title and labels
title_str = 'Azimuthally Averaged Wind Profile';
if condition_result == 1
    title_str = sprintf('%s - CONCENTRIC EYEWALL DETECTED (1)', title_str);
else
    title_str = sprintf('%s - NO CONCENTRIC EYEWALL (0)', title_str);
end

title(title_str, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Radial Distance from Center (km)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Wind Speed (knots)', 'FontSize', 12, 'FontWeight', 'bold');

% Add grid
grid on;

% Add legend
legend('Location', 'best', 'FontSize', 10);

% Set axis limits
xlim([0, max(radial_dist_km)]);
y_limits = [0, max(aa)*1.2];
ylim(y_limits);

% Add statistics text
stats_text = sprintf('Mean: %.1f±%.1f knots\nMax: %.1f knots at %.0fkm\nCondition: %d', ...
    mean(aa, 'omitnan'), std(aa, 'omitnan'), ...
    max_wind, max_dist, condition_result);

if condition_result == 1
    stats_text = sprintf('%s\nPrimary Eyewall: %.1fkm (%.1f knots)\nSecondary Eyewall: %.1fkm (%.1f knots)\nMoat: %.1fkm (%.1f knots)', ...
        stats_text, primary_info(1), primary_info(2), ...
        secondary_info(1), secondary_info(2), moat_info(1), moat_info(2));
end

text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 10, ...
    'BackgroundColor', [1 1 1 0.7], ...
    'FontWeight', 'bold');

% Set axis properties
ax = gca;
set(ax, ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'TickLength', [0.02 0.02], ...
    'LineWidth', 1.5, ...
    'Box', 'on');

hold off;

% ========== 4. Display summary ==========
fprintf('\n========== SUMMARY ==========\n');
fprintf('Azimuthally Averaged Profile Analysis:\n');
fprintf('Condition: %d\n', condition_result);
if condition_result == 1
    fprintf('Primary Eyewall: %.1f km (%.1f knots)\n', eyewall_data(1), eyewall_data(2));
    fprintf('Secondary Eyewall: %.1f km (%.1f knots)\n', eyewall_data(3), eyewall_data(4));
    fprintf('Moat: %.1f km (%.1f knots)\n', moat_data(1), moat_data(2));
else
    fprintf('No concentric eyewall structure detected.\n');
end
fprintf('=============================\n');

% Set figure background to white
set(gcf, 'Color', 'w');
end
