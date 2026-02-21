function varargout = plot_wind_profile_with_coverage(radial_dist_km, azi_wsp, non_nan_percentage, tit)
% PLOT_WIND_PROFILE_WITH_COVERAGE Create a dual-axis plot of wind speed and data coverage
%
% Syntax:
%   plot_wind_profile_with_coverage(radial_dist_km, azi_wsp, non_nan_percentage)
%
% Inputs:
%   radial_dist_km      - Radial distance from center in km (vector)
%   azi_wsp             - Azimuthal wind speed in m/s (vector, same size as radial_dist_km)
%   non_nan_percentage  - Percentage of data coverage (vector, same size as radial_dist_km)
%
% Outputs:
%   Creates a figure with dual y-axes showing wind speed (left) and data coverage (right)
%
% Example:
%   radial_dist_km = 0:0.1:100;
%   azi_wsp = 50 * exp(-radial_dist_km/30) + 20;
%   non_nan_percentage = 100 * (1 - exp(-radial_dist_km/50));
%   plot_wind_profile_with_coverage(radial_dist_km, azi_wsp, non_nan_percentage)

    % Input validation
    if nargin < 3
        error('Function requires three input arguments: radial_dist_km, azi_wsp, non_nan_percentage');
    end
    
    if ~isequal(size(radial_dist_km), size(azi_wsp), size(non_nan_percentage))
        error('All input vectors must have the same size');
    end
    
    % Ensure column vectors
    radial_dist_km = radial_dist_km(:);
    azi_wsp = azi_wsp(:);
    non_nan_percentage = non_nan_percentage(:);
    
    % Create maximized figure
    %fig = figure('WindowState', 'maximized');
     fig = figure('Position', [100, 100, 1400, 800]);
    
    % Left y-axis: Wind speed
    yyaxis left;
    h1 = plot(radial_dist_km, azi_wsp, 'b-', 'LineWidth', 3);
    ylabel('Azimuthal Wind Speed (m s^{-1})', 'FontSize', 14, 'Color', 'b', 'FontWeight', 'bold');
    
    % Auto-adjust ylim with 5% margin
    maxw = max(azi_wsp(:));
    %ylim([0 maxw * 1.05]);  % 5% margin above max value
    grid on;
    
    % Set left axis color to blue
    ax = gca;
    ax.YAxis(1).Color = 'b';           % Left y-axis color
    ax.YAxis(1).Label.Color = 'b';     % Left y-axis label color
    
    % Right y-axis: Data coverage percentage
    yyaxis right;
    h2 = plot(radial_dist_km, non_nan_percentage, 'k--', 'LineWidth', 3);
    ylabel('Percentage of data coverage (%)', 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    ylim([0 100]);
    
    % Set right axis color to black
    ax.YAxis(2).Color = 'k';           % Right y-axis color
    ax.YAxis(2).Label.Color = 'k';     % Right y-axis label color
    
    % Simple shaded area (no hatching) for data coverage
    hold on;
    x_fill = [radial_dist_km; flipud(radial_dist_km)];
    y_fill = [non_nan_percentage; zeros(size(non_nan_percentage))];
    h_fill = fill(x_fill, y_fill, [0.6 0.6 0.6], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
        'DisplayName', 'Data Completeness Area');
    
    % Labels and title
    xlabel('Radial Distance from Center (km)', 'FontSize', 14, 'FontWeight', 'bold');
    title(strcat('Azimuthal Wind Profile and Data Coverage (',string(tit),')'), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Legend (excluding the shaded area for cleaner look)
    legend([h1, h2], {'Azimuthally averaged wind profile', 'Percentage of data coverage'}, ...
        'Location', 'best', 'FontSize', 12, 'FontWeight', 'bold');
    
    % Additional option: Also change tick label colors for consistency
   % ax.YAxis(1).TickLabelColor = 'b';  % Left tick labels in blue
    %ax.YAxis(2).TickLabelColor = 'k';  % Right tick labels in black
    
    % Set all axis properties at once
    set(ax, ...
        'FontSize', 14, ...            % Tick label size for all axes
        'FontWeight', 'bold', ...      % Tick label bold for all axes
        'TickLength', [0.01 0.01], ... % Tick length for all axes
        'LineWidth', 1.5, ...          % Axis line width for all axes
        'Box', 'on', ...
        'TickDir', 'out', ...
        'Position', [0.15, 0.15, 0.6, 0.75]); % Plot area position
    
    % Optional: Save figure (uncomment and modify path as needed)
    % print(fig, 'D:\SRI_NEW\SRI\Plots\Coastal_buffer.tiff', '-dtiff', '-r300');
    
    % Ensure all changes are applied
    drawnow;
    
    % Optional: Return handles if needed
    if nargout >= 1
        varargout{1} = fig;
    end
    if nargout >= 2
        varargout{2} = ax;
    end
    if nargout >= 3
        varargout{3} = [h1, h2, h_fill];
    end
end
