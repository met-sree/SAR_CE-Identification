function [fig, condition_results, eyewall_data, moat_data, sector_coverage_percentage] = analyze_azimuthal_profiles(wsp_sar_polar, radial_dist_km, varargin)
% Analyze azimuthal profiles and plot features for 8 quadrants
% 
% Input:
%   wsp_sar_polar: SAR polar wind speed data in KNOTS (angle × radial distance)
%   radial_dist_km: radial distance array (km)
%   Optional Name-Value pairs:
%     'PlotResults' - Boolean, whether to create plots (default: true)
%     'MinWindKnots' - Minimum wind speed threshold (default: 17)
%     'MinProminence' - Minimum prominence for detection (default: 4)
%     'MinSeparation' - Minimum separation for detection (default: 5)
%     'OutlierThreshold' - Threshold for secondary eyewall outlier detection (default: 50)
%     'MaxWindTolerance' - Tolerance for comparing wind speeds (default: 0.1)
% 
% Output:
%   fig: Figure handle (empty if PlotResults=false)
%   condition_results: 8×1 array with 1 if quadrant meets all conditions, 0 otherwise
%   eyewall_data: 8×4 matrix with [primary_dist, primary_wind; secondary_dist, secondary_wind; NaN, NaN] for each quadrant
%   moat_data: 8×2 matrix with [moat_dist, moat_wind] for each quadrant (NaN if no moat)
%   sector_coverage_percentage: 8×n matrix with coverage percentage for each radial distance

% ========== 1. Parse Input Parameters ==========
p = inputParser;
addRequired(p, 'wsp_sar_polar', @isnumeric);
addRequired(p, 'radial_dist_km', @isnumeric);
addParameter(p, 'PlotResults', true, @islogical);
addParameter(p, 'MinWindKnots', 17, @isnumeric);
addParameter(p, 'MinProminence', 4, @isnumeric);
addParameter(p, 'MinSeparation', 5, @isnumeric);
addParameter(p, 'OutlierThreshold', 50, @isnumeric);
addParameter(p, 'MaxWindTolerance', 0.1, @isnumeric);
parse(p, wsp_sar_polar, radial_dist_km, varargin{:});

% Extract parameters
PLOT_RESULTS = p.Results.PlotResults;
MIN_WIND_KNOTS = p.Results.MinWindKnots;
MIN_PROMINENCE = p.Results.MinProminence;
MIN_SEPARATION = p.Results.MinSeparation;
OUTLIER_THRESHOLD = p.Results.OutlierThreshold;
MAX_WIND_TOLERANCE = p.Results.MaxWindTolerance;

% ========== 2. Data Validation ==========
fprintf('Validating input data...\n');
if size(wsp_sar_polar, 2) ~= length(radial_dist_km)
    error('Radial distance dimension mismatch! wsp_sar_polar size: %d×%d, radial_dist_km length: %d', ...
        size(wsp_sar_polar, 1), size(wsp_sar_polar, 2), length(radial_dist_km));
end

if size(wsp_sar_polar, 1) ~= 360
    warning('Expected 360 azimuthal angles, found %d. Adjusting quadrant calculations.', size(wsp_sar_polar, 1));
end

% ========== 3. Data Processing Functions ==========

    % ========== Calculate azimuthal profiles ==========
    function [azm_profile, azm_profile_std, sector_coverage_percentage] = calculate_azimuthal_profiles(wsp_sar_polar)
        fprintf('Calculating azimuthal profiles...\n');
        
        num_angles = size(wsp_sar_polar, 1);
        num_radial = size(wsp_sar_polar, 2);
        quadrant_size = floor(num_angles / 8);
        
        azm_profile = zeros(num_radial, 8);
        azm_profile_std = zeros(num_radial, 8);
        sector_coverage_percentage = zeros(num_radial, 8);
        
        for quad = 1:8
            % Calculate start and end indices for this quadrant
            start_idx = (quad - 1) * quadrant_size + 1;
            end_idx = min(quad * quadrant_size, num_angles);
            
            % Extract data for current quadrant
            quadrant_data = wsp_sar_polar(start_idx:end_idx, :);
            
           % Calculate coverage percentage
           coverage_percent = (sum(~isnan(quadrant_data), 1) / size(quadrant_data, 1))' * 100;
           sector_coverage_percentage(:, quad) = coverage_percent;
    
           % Calculate mean and std for all data first
           temp_mean = nanmean(quadrant_data, 1)';
           temp_std = std(quadrant_data, 0, 1, 'omitnan')';
    
          % Set to NaN where coverage < 60%
          insufficient_coverage = coverage_percent < 30;
          temp_mean(insufficient_coverage) = NaN;
          temp_std(insufficient_coverage) = NaN;
    
          % Store results
          azm_profile(:, quad) = temp_mean;
          azm_profile_std(:, quad) = temp_std;	    
        end
    end

    % ========== Detect features for a single quadrant ==========
    function [condition_result, eyewall_info, moat_info, features] = detect_quadrant_features(aa, bb, radial_dist_km, quad_num)
        % Initialize outputs
        condition_result = 0;
        eyewall_info = struct('primary', struct('dist', NaN, 'wind', NaN, 'std', NaN, 'is_max_wind', false), ...
                             'secondary', struct('dist', NaN, 'wind', NaN, 'std', NaN, 'is_max_wind', false));
        moat_info = struct('dist', NaN, 'wind', NaN);
        features = struct();
        
        % Find maximum wind speed location
        [max_wind, max_idx] = max(aa);
        max_dist = radial_dist_km(max_idx);
        
        % Calculate eyewall positions (local maxima)
        isEyewall = islocalmax(aa, ...
            'MinProminence', MIN_PROMINENCE, ...
            'MinSeparation', MIN_SEPARATION, ...
            'SamplePoints', radial_dist_km);
        eyewall_indices = find(isEyewall);
        eyewall_values = aa(isEyewall);
        
        % Calculate moat positions (local minima)
        isMoat = islocalmin(aa, ...
            'MinProminence', MIN_PROMINENCE, ...
            'MinSeparation', MIN_SEPARATION, ...
            'SamplePoints', radial_dist_km);
        moat_indices = find(isMoat);
        moat_values = aa(isMoat);
        
        % Store raw features
        features.eyewall_positions = radial_dist_km(eyewall_indices);
        features.moat_positions = radial_dist_km(moat_indices);
        features.eyewall_winds = eyewall_values;
        features.moat_winds = moat_values;
        features.max_wind = max_wind;
        features.max_dist = max_dist;
        features.mean_wind = mean(aa, 'omitnan');
        features.std_wind = std(aa, 'omitnan');
        
        % Filter eyewalls: Remove those with wind speed < MIN_WIND_KNOTS
        valid_eyewall_mask = eyewall_values >= MIN_WIND_KNOTS;
        eyewall_indices = eyewall_indices(valid_eyewall_mask);
        eyewall_values = eyewall_values(valid_eyewall_mask);
        features.eyewall_positions = features.eyewall_positions(valid_eyewall_mask);
        features.eyewall_winds = features.eyewall_winds(valid_eyewall_mask);
        
        % Check if we have at least 2 eyewalls AND at least 1 moat
        if length(features.eyewall_positions) >= 2 && ~isempty(features.moat_positions)
            
            % Sort eyewalls by distance (ascending = inner to outer)
            [sorted_eyewalls, sort_idx] = sort(features.eyewall_positions);
            sorted_eyewall_winds = features.eyewall_winds(sort_idx);
            
            % Check all consecutive eyewall pairs
            for i = 1:(length(sorted_eyewalls) - 1)
                inner_eyewall_dist = sorted_eyewalls(i);
                outer_eyewall_dist = sorted_eyewalls(i + 1);
                inner_eyewall_wind = sorted_eyewall_winds(i);
                outer_eyewall_wind = sorted_eyewall_winds(i + 1);
                
                % Find moats between this specific pair
                moat_idx = find(features.moat_positions > inner_eyewall_dist & ...
                              features.moat_positions < outer_eyewall_dist);
                
                if ~isempty(moat_idx)
                    % Get moat info (take the first moat)
                    moat_dist = features.moat_positions(moat_idx(1));
                    moat_wind = features.moat_winds(moat_idx(1));
                    
                    % Get standard deviations at eyewall positions
                    if length(eyewall_indices) >= 2
                        primary_std = bb(eyewall_indices(sort_idx(i)));
                        secondary_std = bb(eyewall_indices(sort_idx(i+1)));
                    else
                        primary_std = NaN;
                        secondary_std = NaN;
                    end
                    
                    % Check the condition: (eyewall wind - eyewall std) > moat wind
                    condition1 = (inner_eyewall_wind - primary_std) > moat_wind;
                    condition2 = (outer_eyewall_wind - secondary_std) > moat_wind;
                    
                    if condition1 && condition2
                        % Check if either eyewall is at maximum wind location
                        inner_is_max = abs(inner_eyewall_wind - max_wind) <= MAX_WIND_TOLERANCE;
                        outer_is_max = abs(outer_eyewall_wind - max_wind) <= MAX_WIND_TOLERANCE;
                        
                        if inner_is_max || outer_is_max
                            condition_result = 1;
                            
                            % ALWAYS: inner eyewall = primary, outer eyewall = secondary
                            % Regardless of which one has max wind
                            eyewall_info.primary.dist = inner_eyewall_dist;
                            eyewall_info.primary.wind = inner_eyewall_wind;
                            eyewall_info.primary.std = primary_std;
                            eyewall_info.primary.is_max_wind = inner_is_max;
                            
                            eyewall_info.secondary.dist = outer_eyewall_dist;
                            eyewall_info.secondary.wind = outer_eyewall_wind;
                            eyewall_info.secondary.std = secondary_std;
                            eyewall_info.secondary.is_max_wind = outer_is_max;
                            
                            % Store moat information
                            moat_info.dist = moat_dist;
                            moat_info.wind = moat_wind;
                            
                            features.eyewall_max_wind_check = struct('inner_is_max', inner_is_max, ...
                                                                    'outer_is_max', outer_is_max, ...
                                                                    'max_wind', max_wind);
                        else
                            % Neither eyewall is at maximum wind location
                            condition_result = 0;
                            features.eyewall_max_wind_check = struct('inner_is_max', inner_is_max, ...
                                                                    'outer_is_max', outer_is_max, ...
                                                                    'max_wind', max_wind, ...
                                                                    'reason', 'neither_eyewall_at_max');
                        end
                        break;
                    end
                end
            end
        end
        
        % If no concentric eyewall found or condition not met, use maximum wind as primary
        if condition_result == 0
            eyewall_info.primary.dist = max_dist;
            eyewall_info.primary.wind = max_wind;
            eyewall_info.primary.std = bb(max_idx);
            eyewall_info.primary.is_max_wind = true;
            eyewall_info.secondary.is_max_wind = false;
        end
        
        % Store additional features
        features.condition_result = condition_result;
        features.eyewall_info = eyewall_info;
        features.moat_info = moat_info;
    end

    % ========== Apply secondary eyewall outlier filter ==========
    function [condition_results, outlier_info] = apply_outlier_filter(condition_results, eyewall_data, outlier_threshold)
        % Apply outlier filter based on secondary eyewall wind speeds
        outlier_info = struct();
        
        % Extract secondary eyewall wind speeds (column 4)
        SE_data = eyewall_data(:, 4);
        
        % Find valid secondary eyewall wind speeds (non-NaN and condition=1)
        valid_mask = ~isnan(SE_data) & condition_results == 1;
        valid_SE_data = SE_data(valid_mask);
        
        outlier_info.original_condition_results = condition_results;
        outlier_info.original_SE_data = SE_data;
        outlier_info.valid_SE_data = valid_SE_data;
        
        if isempty(valid_SE_data) || length(valid_SE_data) < 2
            % Not enough data for outlier detection
            outlier_info.median_val = NaN;
            outlier_info.idx_outliers = [];
            outlier_info.num_outliers = 0;
            outlier_info.outlier_quadrants = [];
            return;
        end
        
        % Calculate median of valid secondary eyewall wind speeds
        median_val = median(valid_SE_data, 'omitnan');
        outlier_info.median_val = median_val;
        
        % Find indices of outliers (outside median ± threshold range)
        outlier_mask = valid_SE_data < (median_val - outlier_threshold) | ...
                      valid_SE_data > (median_val + outlier_threshold);
        
        % Get the actual quadrant indices of outliers
        valid_quadrants = find(valid_mask);
        outlier_quadrants = valid_quadrants(outlier_mask);
        
        outlier_info.idx_outliers = outlier_quadrants;
        outlier_info.num_outliers = length(outlier_quadrants);
        outlier_info.outlier_quadrants = outlier_quadrants;
        
        % Update condition results for outliers
        if outlier_info.num_outliers > 0
            condition_results(outlier_quadrants) = 0;
            outlier_info.updated_condition_results = condition_results;
        else
            outlier_info.updated_condition_results = condition_results;
        end
    end

    % ========== Plot quadrant results with coverage shading ==========
    function plot_quadrant_results(ax, quad, radial_dist_km, aa, azm_profile_std, coverage_percent, features, is_outlier)
        % Plot base curve
        h_wind = plot(ax, radial_dist_km, aa, 'b-', 'LineWidth', 1.5);
        hold(ax, 'on');
        
        % Add standard deviation shading
        upper_bound = aa + azm_profile_std(:, quad);
        lower_bound = aa - azm_profile_std(:, quad);
        fill(ax, [radial_dist_km, fliplr(radial_dist_km)], ...
            [upper_bound', fliplr(lower_bound')], ...
            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Create second y-axis for coverage percentage
        yyaxis(ax, 'right');
        
        % Plot coverage percentage as gray line with shading
        h_coverage = plot(ax, radial_dist_km, coverage_percent(:, quad), ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'LineStyle', '-', ...
            'DisplayName', 'Coverage');
        
        % Add gray shading under coverage line
        fill(ax, [radial_dist_km, fliplr(radial_dist_km)], ...
            [coverage_percent(:, quad)', zeros(1, length(coverage_percent(:, quad)))], ...
            [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
            'DisplayName', 'Coverage Area');
        
        % Configure right y-axis
        ylim(ax, [0 100]);
        ax.YAxis(2).Color = [0.3 0.3 0.3];
        
        % Determine column position (1 or 2)
        col = mod(quad-1, 2) + 1;
        
        % Configure axis labels based on column
        if col == 1
            % First column: Show left y-axis labels, hide right y-axis labels
            ylabel(ax, 'Wind Speed (knots)', 'FontSize', 10, 'FontWeight', 'bold');
            ax.YAxis(2).Label.String = '';
            ax.YAxis(2).TickLabels = {};
        else
            % Second column: Show right y-axis labels, hide left y-axis labels
            ylabel(ax, 'Coverage (%)', 'FontSize', 10, 'FontWeight', 'bold');
            ax.YAxis(1).Label.String = '';
            ax.YAxis(1).TickLabels = {};
        end
        
        grid(ax, 'on');
        
        % Switch back to left y-axis
        yyaxis(ax, 'left');
        
        % Plot symbols based on condition
        alpha = 0.5;
        
        if features.condition_result == 1 && ~is_outlier
            % Plot concentric eyewall features (if not an outlier)
            % PRIMARY EYEWALL: Magenta asterisk (always, regardless of max wind)
            scatter(ax, features.eyewall_info.primary.dist, features.eyewall_info.primary.wind, ...
                100, 'm*', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, ...
                'LineWidth', 2, 'DisplayName', 'PE');
            
            % SECONDARY EYEWALL: Red asterisk (always, regardless of max wind)
            scatter(ax, features.eyewall_info.secondary.dist, features.eyewall_info.secondary.wind, ...
                100, 'r*', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, ...
                'LineWidth', 2, 'DisplayName', 'SE');
            
            % MOAT: Green asterisk
            scatter(ax, features.moat_info.dist, features.moat_info.wind, ...
                100, 'g*', 'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, ...
                'LineWidth', 2, 'DisplayName', 'Moat');
        else
            % Plot only primary eyewall (maximum wind location) - Magenta asterisk
            scatter(ax, features.max_dist, features.max_wind, 50, 'm*', ...
                'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha, ...
                'LineWidth', 2, 'DisplayName', 'PE');
        end
        
        % Set title with additional info
        title_str = sprintf('Sector %d', quad);
        if features.condition_result == 1
            title_str = sprintf('%s - (CE Flag = 1)', title_str);
        else
            title_str = sprintf('%s - (CE Flag = 0)', title_str);
        end
        
        text(ax, 0.02, 0.98, title_str, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 8, ...
            'BackgroundColor', [1 1 1 0.7], 'FontWeight', 'bold');
        
        % Configure left y-axis
        xlim(ax, [0, max(radial_dist_km)]);
        y_limits = [0, max(aa)*1.2];
        %ylim(ax, 300);
        ax.YAxis(1).Color = 'b';
        
        % Set x-axis labels based on row position
        row = ceil(quad/2);
        
        if row == 4
            xlabel(ax, 'Radial Distance (km)', 'FontSize', 10, 'FontWeight', 'bold');
        else
            xlabel(ax, '');
            set(ax, 'XTickLabel', []);
        end

         if col == 1
            % First column: Show left y-axis labels, hide right y-axis labels
            ylabel(ax, 'Wind Speed (knots)', 'FontSize', 10, 'FontWeight', 'bold');
            % ax.YAxis(2).Label.String = '';
            % ax.YAxis(2).TickLabels = {};
         end
        
        % Update wind speed line display name
        set(h_wind, 'DisplayName', 'Vmax');
        
        % Create legend (only left y-axis items)
        if features.condition_result == 1 && ~is_outlier
            % For concentric eyewall: Vmax, PE, SE, Moat
            legend_items = [h_wind, ...
                           findobj(ax.Children, 'DisplayName', 'PE'), ...
                           findobj(ax.Children, 'DisplayName', 'SE'), ...
                           findobj(ax.Children, 'DisplayName', 'Moat')];
            legend_names = {'Vmax', 'PE', 'SE', 'Moat'};
            
            legend(ax, legend_items, legend_names, ...
                   'Location', 'southwest', 'FontSize', 8, 'NumColumns', 2, ...
                   'Box', 'off', 'FontWeight', 'bold');
        else
            % For single eyewall: Vmax, PE
            legend(ax, [h_wind, findobj(ax.Children, 'DisplayName', 'PE')], ...
                   {'Vmax', 'PE'}, ...
                   'Location', 'southwest', 'FontSize', 8, 'NumColumns', 2, ...
                   'Box', 'off', 'FontWeight', 'bold');
        end
        
        % Add statistics text
        stats_text = sprintf('Mean: %.1f±%.1f knots\nMax: %.1f knots at %.0fkm\nCondition: %d', ...
            features.mean_wind, features.std_wind, ...
            features.max_wind, features.max_dist, features.condition_result);
        
        if features.condition_result == 1
            if is_outlier
                stats_text = sprintf('%s\n[OUTLIER FILTERED]', stats_text);
            else
                max_location = '';
                if features.eyewall_info.primary.is_max_wind
                    max_location = 'Primary';
                elseif features.eyewall_info.secondary.is_max_wind
                    max_location = 'Secondary';
                end
                
                stats_text = sprintf('%s\nMax wind location: %s eyewall\nPrimary: %.1fkm (%.1f knots)\nSecondary: %.1fkm (%.1f knots)\nMoat: %.1fkm (%.1f knots)', ...
                    stats_text, max_location, ...
                    features.eyewall_info.primary.dist, features.eyewall_info.primary.wind, ...
                    features.eyewall_info.secondary.dist, features.eyewall_info.secondary.wind, ...
                    features.moat_info.dist, features.moat_info.wind);
            end
        elseif isfield(features, 'eyewall_max_wind_check') && ...
               isfield(features.eyewall_max_wind_check, 'reason') && ...
               strcmp(features.eyewall_max_wind_check.reason, 'neither_eyewall_at_max')
            stats_text = sprintf('%s\n[Neither eyewall at max wind]', stats_text);
        end
        
        text(ax, 0.58, 0.98, stats_text, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 8, ...
            'FontWeight', 'bold');
        
        hold(ax, 'off');
    end

    % ========== Create summary figure ==========
    function fig = create_summary_figure(radial_dist_km, azm_profile, azm_profile_std, coverage_percent, all_features, outlier_quadrants)
        % Create figure window
        fig = figure('WindowState', 'maximized');
        
        % Store subplot handles
        ax_handles = gobjects(8, 1);
        
        % Create subplots for each quadrant
        for quad = 1:8
            ax_handles(quad) = subplot(4, 2, quad);
            
            % Check if this quadrant is an outlier
            is_outlier = ismember(quad, outlier_quadrants);
            
            plot_quadrant_results(ax_handles(quad), quad, radial_dist_km, ...
                azm_profile(:, quad), azm_profile_std, coverage_percent, ...
                all_features{quad}, is_outlier);
        end
        
        % Define custom subplot positions
        subplot_positions = {
            [0.15, 0.72, 0.35, 0.2],  % Q1
            [0.53, 0.72, 0.35, 0.2],  % Q2
            [0.15, 0.5,  0.35, 0.2],  % Q3
            [0.53, 0.5,  0.35, 0.2],  % Q4
            [0.15, 0.28, 0.35, 0.2],  % Q5
            [0.53, 0.28, 0.35, 0.2],  % Q6
            [0.15, 0.06, 0.35, 0.2],  % Q7
            [0.53, 0.06, 0.35, 0.2]   % Q8
        };
        
        % Apply custom positions
        for quad = 1:8
            set(ax_handles(quad), 'Position', subplot_positions{quad}, 'FontWeight', 'bold');
            set(ax_handles(quad), 'LooseInset', get(ax_handles(quad), 'TightInset'));
        end
        
        % Add main title
        sgtitle('Average azimuthal wind profiles in each 45 degree sector', ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        % Set figure properties
        set(fig, 'Color', 'w');
        ha = findobj(fig, 'type', 'axes');
        set(ha, 'FontSize', 9);
    end

    % ========== Display analysis results ==========
    function display_analysis_results(condition_results, eyewall_data, moat_data, all_features, outlier_info)
        % Display quadrant-by-quadrant results
        fprintf('\n\n========== SUMMARY OF CONCENTRIC EYEWALL DETECTION ==========\n');
        for quad = 1:8
            if condition_results(quad) == 1
                if ismember(quad, outlier_info.outlier_quadrants)
                    fprintf('Quadrant %d: Condition = 1 [FILTERED - OUTLIER]\n', quad);
                else
                    fprintf('Quadrant %d: Condition = 1', quad);
                    if all_features{quad}.eyewall_info.primary.is_max_wind
                        fprintf(' (Primary at max wind)\n');
                    elseif all_features{quad}.eyewall_info.secondary.is_max_wind
                        fprintf(' (Secondary at max wind)\n');
                    else
                        fprintf(' (Neither eyewall at max wind)\n');
                    end
                end
                fprintf('  Primary: %.1f km (%.1f knots)\n', eyewall_data(quad,1), eyewall_data(quad,2));
                fprintf('  Secondary: %.1f km (%.1f knots)\n', eyewall_data(quad,3), eyewall_data(quad,4));
                fprintf('  Moat: %.1f km (%.1f knots)\n', moat_data(quad,1), moat_data(quad,2));
            else
                fprintf('Quadrant %d: Condition = 0', quad);
                if ismember(quad, outlier_info.outlier_quadrants)
                    fprintf(' [OUTLIER FILTERED]\n');
                elseif isfield(all_features{quad}, 'eyewall_max_wind_check') && ...
                       isfield(all_features{quad}.eyewall_max_wind_check, 'reason') && ...
                       strcmp(all_features{quad}.eyewall_max_wind_check.reason, 'neither_eyewall_at_max')
                    fprintf(' [Neither eyewall at max wind]\n');
                else
                    fprintf(' (No valid data)\n');
                end
            end
        end
        
        % Display outlier detection information
        if ~isempty(outlier_info.outlier_quadrants)
            fprintf('\n========== OUTLIER DETECTION RESULTS ==========\n');
            fprintf('Median of secondary eyewall wind speeds: %.1f knots\n', outlier_info.median_val);
            fprintf('Outlier threshold: ±%.1f knots from median\n', OUTLIER_THRESHOLD);
            fprintf('Valid range: %.1f to %.1f knots\n', ...
                outlier_info.median_val - OUTLIER_THRESHOLD, ...
                outlier_info.median_val + OUTLIER_THRESHOLD);
            fprintf('Number of outliers detected: %d\n', outlier_info.num_outliers);
            fprintf('Outlier quadrants: %s\n', mat2str(outlier_info.outlier_quadrants'));
            fprintf('================================================\n');
        end
        
        % Summary table
        fprintf('\n========== DATA SUMMARY TABLE (All values in knots) ==========\n');
        fprintf('Quadrant | Condition | Primary (km, knots) | Secondary (km, knots) | Moat (km, knots) | Max Wind Location\n');
        fprintf('---------|-----------|---------------------|-----------------------|------------------|------------------\n');
        for quad = 1:8
            if condition_results(quad) == 1
                max_location = '';
                if all_features{quad}.eyewall_info.primary.is_max_wind
                    max_location = 'Primary';
                elseif all_features{quad}.eyewall_info.secondary.is_max_wind
                    max_location = 'Secondary';
                end
                
                if ismember(quad, outlier_info.outlier_quadrants)
                    fprintf('   %2d    |    1*     | %5.1f, %5.1f     | %6.1f, %5.1f       | %5.1f, %5.1f | %s*\n', ...
                        quad, eyewall_data(quad,1), eyewall_data(quad,2), ...
                        eyewall_data(quad,3), eyewall_data(quad,4), ...
                        moat_data(quad,1), moat_data(quad,2), max_location);
                else
                    fprintf('   %2d    |     1     | %5.1f, %5.1f     | %6.1f, %5.1f       | %5.1f, %5.1f | %s\n', ...
                        quad, eyewall_data(quad,1), eyewall_data(quad,2), ...
                        eyewall_data(quad,3), eyewall_data(quad,4), ...
                        moat_data(quad,1), moat_data(quad,2), max_location);
                end
            else
                fprintf('   %2d    |     0     |       NaN, NaN        |       NaN, NaN          |       NaN, NaN | -\n', quad);
            end
        end
        fprintf('========================================================================\n');
        
        % Final classification
        num_valid_quadrants = sum(condition_results);
        fprintf('\n========== FINAL CLASSIFICATION ==========\n');
        if num_valid_quadrants >= 5
            fprintf('✓ SYMMETRIC CONCENTRIC EYEWALL STRUCTURE DETECTED\n');
            fprintf('  %d out of 8 quadrants meet all conditions (≥5 required)\n', num_valid_quadrants);
            
            % Check radial distance between outer eyewalls
            valid_quadrants = find(condition_results == 1);
            outer_eyewall_distances = eyewall_data(valid_quadrants, 3);
            
            if length(outer_eyewall_distances) >= 2
                max_distance_diff = max(outer_eyewall_distances) - min(outer_eyewall_distances);
                fprintf('  Maximum radial distance between outer eyewalls: %.1f km\n', max_distance_diff);
                
                if max_distance_diff < 50
                    fprintf('  ✓ Outer eyewalls are coherent (distance < 50 km)\n');
                else
                    fprintf('  ✗ Outer eyewalls may not be coherent (distance ≥ 50 km)\n');
                end
            end
        else
            fprintf('✗ NO SYMMETRIC CONCENTRIC EYEWALL STRUCTURE\n');
            fprintf('  Only %d out of 8 quadrants meet all conditions (<5 required)\n', num_valid_quadrants);
        end
        fprintf('==========================================\n');
    end

% ========== 4. Main Data Processing ==========

% Calculate azimuthal profiles
[azm_profile, azm_profile_std, sector_coverage_percentage] = calculate_azimuthal_profiles(wsp_sar_polar);

% Initialize storage
all_features = cell(8, 1);
condition_results = zeros(8, 1);
eyewall_data = nan(8, 4);
moat_data = nan(8, 2);

% Detect features for each quadrant
fprintf('\nDetecting eyewall and moat features...\n');
fprintf('Filtering criterion: Eyewalls AND moats must have wind speed >= %.1f knots\n', MIN_WIND_KNOTS);
fprintf('Max wind condition: At least one eyewall must be at maximum wind location (tolerance: ±%.1f knots)\n', MAX_WIND_TOLERANCE);

for quad = 1:8
    fprintf('\n=== Processing Quadrant %d ===\n', quad);
    
    % Get data for current quadrant
    aa = azm_profile(:, quad);
    bb = azm_profile_std(:, quad);
    
    % Detect features
    [condition_result, eyewall_info, moat_info, features] = ...
        detect_quadrant_features(aa, bb, radial_dist_km, quad);
    
    % Store results
    all_features{quad} = features;
    condition_results(quad) = condition_result;
    
    % Store in output matrices
    eyewall_data(quad, :) = [eyewall_info.primary.dist, eyewall_info.primary.wind, ...
                            eyewall_info.secondary.dist, eyewall_info.secondary.wind];
    moat_data(quad, :) = [moat_info.dist, moat_info.wind];
    
    % Display quadrant results
    if condition_result == 1
        fprintf('✓ CONCENTRIC EYEWALL DETECTED\n');
        fprintf('  Primary (inner): %.1f km, %.1f knots (std: %.1f)\n', ...
            eyewall_info.primary.dist, eyewall_info.primary.wind, eyewall_info.primary.std);
        fprintf('  Secondary (outer): %.1f km, %.1f knots (std: %.1f)\n', ...
            eyewall_info.secondary.dist, eyewall_info.secondary.wind, eyewall_info.secondary.std);
        fprintf('  Moat: %.1f km, %.1f knots\n', moat_info.dist, moat_info.wind);
        
        if eyewall_info.primary.is_max_wind
            fprintf('  ✓ Primary eyewall is at maximum wind location (%.1f knots)\n', features.max_wind);
        elseif eyewall_info.secondary.is_max_wind
            fprintf('  ✓ Secondary eyewall is at maximum wind location (%.1f knots)\n', features.max_wind);
        end
    else
        fprintf('✗ No concentric eyewall detected\n');
        if isfield(features, 'eyewall_max_wind_check') && ...
           isfield(features.eyewall_max_wind_check, 'reason') && ...
           strcmp(features.eyewall_max_wind_check.reason, 'neither_eyewall_at_max')
            fprintf('  Reason: Neither eyewall is at maximum wind location (Max wind: %.1f knots)\n', features.max_wind);
        elseif ~isempty(features.eyewall_positions) && length(features.eyewall_positions) < 2
            fprintf('  Reason: Not enough eyewalls (found %d, need ≥2)\n', length(features.eyewall_positions));
        elseif ~isempty(features.eyewall_positions) && isempty(features.moat_positions)
            fprintf('  Reason: No moats found between eyewalls\n');
        elseif ~isempty(features.eyewall_positions) && ~isempty(features.moat_positions)
            fprintf('  Reason: Conditions not met (wind difference too small)\n');
        else
            fprintf('  Reason: Insufficient features detected\n');
        end
    end
end

% ========== 5. Apply Outlier Filter ==========
fprintf('\nApplying secondary eyewall outlier filter...\n');
fprintf('Outlier threshold: ±%.1f knots from median\n', OUTLIER_THRESHOLD);

[condition_results, outlier_info] = apply_outlier_filter(condition_results, eyewall_data, OUTLIER_THRESHOLD);

% Display outlier detection results
if outlier_info.num_outliers > 0
    fprintf('Detected %d outlier(s): Quadrants %s\n', ...
        outlier_info.num_outliers, mat2str(outlier_info.outlier_quadrants'));
    fprintf('Median secondary eyewall wind speed: %.1f knots\n', outlier_info.median_val);
    fprintf('Valid range: %.1f to %.1f knots\n', ...
        outlier_info.median_val - OUTLIER_THRESHOLD, ...
        outlier_info.median_val + OUTLIER_THRESHOLD);
    
    % Display outlier details
    for i = 1:length(outlier_info.outlier_quadrants)
        quad = outlier_info.outlier_quadrants(i);
        SE_wind = eyewall_data(quad, 4);
        fprintf('  Quadrant %d: Secondary eyewall wind = %.1f knots (outside valid range)\n', ...
            quad, SE_wind);
    end
else
    fprintf('No outliers detected.\n');
end

% ========== 6. Plot Results (if requested) ==========
if PLOT_RESULTS
    fprintf('\nCreating summary figure...\n');
    fig = create_summary_figure(radial_dist_km, azm_profile, azm_profile_std, ...
        sector_coverage_percentage, all_features, outlier_info.outlier_quadrants);
else
    fig = [];
end

% ========== 7. Display Summary ==========
display_analysis_results(condition_results, eyewall_data, moat_data, all_features, outlier_info);

fprintf('\nAnalysis complete!\n');

end
