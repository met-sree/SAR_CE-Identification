function [wsp_sar_polar, xx_geo, yy_geo] = cart_to_polar(sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat)
% CREATE_POLAR_INTERPOLATION Convert scattered SAR data to polar grid
%
% Syntax:
%   [wsp_sar_polar, xx_geo, yy_geo] = create_polar_interpolation(sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat)
%
% Inputs:
%   sar_lon    - Longitude coordinates of SAR data (scattered)
%   sar_lat    - Latitude coordinates of SAR data (scattered)
%   wsp_sar    - Wind speed values at scattered points
%   maxrds     - Maximum radial distance in kilometers
%   res        - Radial resolution in kilometers
%   cen_lon    - Center longitude of polar grid
%   cen_lat    - Center latitude of polar grid
%
% Outputs:
%   wsp_sar_polar - Interpolated wind speed on polar grid
%   xx_geo        - Longitude coordinates of polar grid (geographic)
%   yy_geo        - Latitude coordinates of polar grid (geographic)
%
% Example:
%   [polar_data, lon_grid, lat_grid] = create_polar_interpolation(sar_lon, sar_lat, wsp_sar, 300, 0.5, 120.5, 23.7);

    % Constants
    KM_PER_DEGREE = 111.1;  % Approximate km per degree at equator
    
    % Validate inputs
    if nargin < 7
        error('All 7 input arguments are required: sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat');
    end
    
    % Calculate number of radial steps
    num_radial_steps = round(maxrds/res) + 1;
    radial_distances_km = 0:res:maxrds;
    
    % Calculate number of azimuthal steps (360 degrees with 1-degree resolution)
    azimuth_steps = 360;
    
    % Initialize polar coordinate arrays
    xx = zeros(azimuth_steps, num_radial_steps);
    yy = zeros(azimuth_steps, num_radial_steps);
    
    % Create polar grid in local coordinates (degrees offset from center)
    rd = 1;
    for rad_km = radial_distances_km
        % Convert km to degrees
        rad_deg = rad_km / KM_PER_DEGREE;
        
        for the = 0:(azimuth_steps)
            xx(the+1, rd) = rad_deg * cosd(the);
            yy(the+1, rd) = rad_deg * sind(the);
        end
        rd = rd + 1;
    end
    
    % Create geographic coordinates by adding center position
    xx_geo = xx + cen_lon;  % xx_geo is what you called xx+cen_lon
    yy_geo = yy + cen_lat;  % yy_geo is what you called yy+cen_lat
    
    % Create scattered interpolant
    F = scatteredInterpolant(double(sar_lon(:)), double(sar_lat(:)), double(wsp_sar(:)));
    F.Method = 'nearest';
    F.ExtrapolationMethod = 'none';
    
    % Interpolate to polar grid
    wsp_sar_polar = F(xx_geo, yy_geo);
    
    % Display summary
    fprintf('Polar interpolation complete:\n');
    fprintf('  Grid size: %d azimuths × %d radii\n', azimuth_steps, num_radial_steps);
    fprintf('  Radial range: 0 to %.1f km (%.1f to %.1f degrees)\n', ...
        maxrds, 0, maxrds/KM_PER_DEGREE);
    fprintf('  Center: (%.4f°E, %.4f°N)\n', cen_lon, cen_lat);
end