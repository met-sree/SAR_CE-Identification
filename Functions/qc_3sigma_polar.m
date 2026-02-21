function [wsp_wind_qc, qc_flag, stats] = qc_3sigma_polar(wsp_wind, varargin)
% QC_3SIGMA_POLAR - Perform 3-sigma quality control on polar wind data
%
% Syntax:
%   [wsp_wind_qc, qc_flag, stats] = qc_3sigma_polar(wsp_wind)
%   [wsp_wind_qc, qc_flag, stats] = qc_3sigma_polar(wsp_wind, 'Name', Value, ...)
%
% Input:
%   wsp_wind - 2D matrix of wind speeds (azimuth × radius)
%              Typically 361×300 or similar (angles × distances)
%
% Optional Name-Value pairs:
%   'SigmaThreshold' - Number of standard deviations (default: 3)
%   'MinDataPoints'  - Minimum valid points needed per radius (default: 5)
%   'Verbose'        - Display progress and statistics (default: true)
%   'Plot'           - Generate diagnostic plots (default: false)
%   'KeepOriginal'   - Return original values for non-outliers (default: true)
%                      If false, sets non-outliers to mean value
%
% Output:
%   wsp_wind_qc - Quality-controlled wind matrix (NaNs for outliers)
%   qc_flag     - Binary flag matrix (1=good, 0=outlier)
%   stats       - Structure containing QC statistics
%
% Example:
%   % Basic usage
%   [wsp_qc, flag, stats] = qc_3sigma_polar(wsp_wind);
%
%   % With custom threshold and plots
%   [wsp_qc, flag, stats] = qc_3sigma_polar(wsp_wind, ...
%       'SigmaThreshold', 2.5, 'Plot', true);
%
% See also: nanmean, nanstd, isoutlier

% Parse input parameters
p = inputParser;
addRequired(p, 'wsp_wind', @isnumeric);
addParameter(p, 'SigmaThreshold', 3, @(x) isnumeric(x) && x>0);
addParameter(p, 'MinDataPoints', 5, @(x) isnumeric(x) && x>1);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'Plot', false, @islogical);
addParameter(p, 'KeepOriginal', true, @islogical);
parse(p, wsp_wind, varargin{:});

% Extract parameters
sigma_thresh = p.Results.SigmaThreshold;
min_points = p.Results.MinDataPoints;
verbose = p.Results.Verbose;
do_plot = p.Results.Plot;
keep_original = p.Results.KeepOriginal;

% Initialize variables
[n_azimuth, n_radius] = size(wsp_wind);
qc_flag = ones(size(wsp_wind));
wsp_wind_qc = wsp_wind;  % Start with original data

% Pre-allocate statistics arrays
stats.mu = NaN(n_radius, 1);
stats.sigma = NaN(n_radius, 1);
stats.lower_thresh = NaN(n_radius, 1);
stats.upper_thresh = NaN(n_radius, 1);
stats.n_outliers = zeros(n_radius, 1);
stats.n_valid_original = zeros(n_radius, 1);

if verbose
    fprintf('Performing 3-sigma QC on %d×%d wind matrix...\n', n_azimuth, n_radius);
    fprintf('Threshold: %.1fσ, Min points per radius: %d\n', sigma_thresh, min_points);
end

% Loop through each radius
for r = 1:n_radius
    % Get wind speeds at current radius
    winds = wsp_wind(:, r);
    
    % Find valid (non-NaN) points
    valid_idx = ~isnan(winds);
    winds_valid = winds(valid_idx);
    stats.n_valid_original(r) = sum(valid_idx);
    
    % Skip if insufficient data
    if length(winds_valid) < min_points
        if verbose && length(winds_valid) > 0
            fprintf('  Radius %d: Insufficient data (%d points, need %d)\n', ...
                r, length(winds_valid), min_points);
        end
        continue;
    end
    
    % Calculate statistics
    mu = mean(winds_valid);
    sigma = std(winds_valid);
    
    % Store statistics
    stats.mu(r) = mu;
    stats.sigma(r) = sigma;
    stats.lower_thresh(r) = mu - sigma_thresh * sigma;
    stats.upper_thresh(r) = mu + sigma_thresh * sigma;
    
    % Check for valid sigma
    if sigma <= 0 || isnan(sigma)
        if verbose
            fprintf('  Radius %d: Invalid sigma (%.2f)\n', r, sigma);
        end
        continue;
    end
    
    % Identify outliers (only among originally valid points)
    is_outlier = false(n_azimuth, 1);
    is_outlier(valid_idx) = winds_valid < stats.lower_thresh(r) | ...
                           winds_valid > stats.upper_thresh(r);
    
    % Count outliers
    stats.n_outliers(r) = sum(is_outlier);
    
    % Update QC flag
    qc_flag(is_outlier, r) = 0;
    
    % Apply QC to data
    if keep_original
        % Replace outliers with NaN (keep good data as-is)
        wsp_wind_qc(is_outlier, r) = NaN;
    else
        % Replace outliers with mean value
        wsp_wind_qc(is_outlier, r) = mu;
    end
end

% Calculate summary statistics
total_points = n_azimuth * n_radius;
total_outliers = sum(stats.n_outliers);
total_valid_original = sum(stats.n_valid_original);
outlier_percentage = 100 * total_outliers / total_valid_original;

% Store summary in stats structure
stats.total_points = total_points;
stats.total_outliers = total_outliers;
stats.total_valid_original = total_valid_original;
stats.outlier_percentage = outlier_percentage;
stats.sigma_threshold = sigma_thresh;

if verbose
    fprintf('\nQC Summary:\n');
    fprintf('  Total data points: %d\n', total_points);
    fprintf('  Originally valid points: %d\n', total_valid_original);
    fprintf('  Outliers removed: %d (%.2f%% of valid data)\n', ...
        total_outliers, outlier_percentage);
    fprintf('  Final valid points: %d\n', sum(~isnan(wsp_wind_qc(:))));
    
    % Display radius with most outliers
    [max_outliers, max_radius] = max(stats.n_outliers);
    if max_outliers > 0
        fprintf('  Radius with most outliers: %d (%d outliers, %.1f%%)\n', ...
            max_radius, max_outliers, 100*max_outliers/n_azimuth);
    end
end

% Generate diagnostic plots if requested
if do_plot
    plot_qc_results(wsp_wind, wsp_wind_qc, qc_flag, stats);
end

end

% -------------------------------------------------------------------------
% Helper function for plotting
% -------------------------------------------------------------------------
function plot_qc_results(original, qced, qc_flag, stats)
% PLOT_QC_RESULTS - Generate diagnostic plots for QC results

figure('Position', [100, 100, 1400, 800], 'Name', '3-sigma QC Diagnostics');

% 1. Original vs QC'ed data
subplot(2, 3, 1);
imagesc(original);
colorbar;
title('Original Wind Speed');
xlabel('Radius index');
ylabel('Azimuth index');
clim([0 70]);  % Adjust based on typical wind speeds

subplot(2, 3, 2);
imagesc(qced);
colorbar;
title('After 3-\sigma QC');
xlabel('Radius index');
ylabel('Azimuth index');
clim([0 70]);

% 2. QC flags
subplot(2, 3, 3);
imagesc(qc_flag);
colorbar;
title(sprintf('QC Flags (1=good, 0=outlier)\nTotal outliers: %d (%.1f%%)', ...
    stats.total_outliers, stats.outlier_percentage));
xlabel('Radius index');
ylabel('Azimuth index');
colormap(gca, [1 0.5 0.5; 0.5 1 0.5]);  % Red for bad, green for good

% 3. Statistics per radius
subplot(2, 3, 4);
hold on;
plot(stats.mu, 'b-', 'LineWidth', 2);
plot(stats.upper_thresh, 'r--', 'LineWidth', 1.5);
plot(stats.lower_thresh, 'r--', 'LineWidth', 1.5);
xlabel('Radius index');
ylabel('Wind speed (m/s)');
title(sprintf('Mean ± %.1fσ at each radius', stats.sigma_threshold));
legend('Mean', sprintf('Mean+%.1fσ', stats.sigma_threshold), ...
    sprintf('Mean-%.1fσ', stats.sigma_threshold), 'Location', 'best');
grid on;

% 4. Outliers per radius
subplot(2, 3, 5);
bar(stats.n_outliers, 'EdgeColor', 'none');
xlabel('Radius index');
ylabel('Number of outliers');
title('Outliers removed per radius');
grid on;

% 5. Percentage of outliers per radius
subplot(2, 3, 6);
percentage = 100 * stats.n_outliers / size(original, 1);
bar(percentage, 'EdgeColor', 'none');
xlabel('Radius index');
ylabel('Outliers (%)');
title('Percentage of outliers per radius');
grid on;
ylim([0 min(100, max(percentage)*1.1)]);

% Add overall title
sgtitle(sprintf('3-Sigma Quality Control (Threshold = %.1fσ)', stats.sigma_threshold), ...
    'FontSize', 12, 'FontWeight', 'bold');
end