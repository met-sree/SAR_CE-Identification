function [wsp_sar_polar, xx_geo, yy_geo] = cart_to_polar(sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat)
% CART_TO_POLAR Convert scattered SAR data to polar grid (GLOBAL VERSION)
%
% This version uses geodesic calculations to work correctly anywhere on Earth.
%
% Syntax:
%   [wsp_sar_polar, xx_geo, yy_geo] = cart_to_polar(sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat)
%
% Inputs:
%   sar_lon    - Longitude coordinates of SAR data (scattered)
%   sar_lat    - Latitude coordinates of SAR data (scattered)
%   sar_lon and sar_lat should be in degrees (standard -180 to 180 or 0 to 360)
%   wsp_sar    - Wind speed values at scattered points
%   maxrds     - Maximum radial distance in kilometers
%   res        - Radial resolution in kilometers
%   cen_lon    - Center longitude of polar grid (degrees)
%   cen_lat    - Center latitude of polar grid (degrees)
%
% Outputs:
%   wsp_sar_polar - Interpolated wind speed on polar grid
%   xx_geo        - Longitude coordinates of polar grid (geographic)
%   yy_geo        - Latitude coordinates of polar grid (geographic)
%
% Example:
%   [polar_data, lon_grid, lat_grid] = cart_to_polar(sar_lon, sar_lat, wsp_sar, 300, 0.5, 120.5, 23.7);

    % Validate inputs
    if nargin < 7
        error('All 7 input arguments are required: sar_lon, sar_lat, wsp_sar, maxrds, res, cen_lon, cen_lat');
    end
    
    % Normalize longitudes to 0-360 range for consistency
    sar_lon = mod(sar_lon, 360);
    cen_lon = mod(cen_lon, 360);
    
    % Calculate number of radial steps
    num_radial_steps = round(maxrds/res) + 1;
    radial_distances_km = 0:res:maxrds;
    
    % Calculate number of azimuthal steps (360 degrees with 1-degree resolution)
    azimuth_steps = 361;
    
    % Initialize output grids
    xx_geo = zeros(azimuth_steps, num_radial_steps);
    yy_geo = zeros(azimuth_steps, num_radial_steps);
    
    % Create polar grid using geodesic calculations
    % Reference ellipsoid (WGS84)
    wgs84 = wgs84Ellipsoid('kilometer');
    
    for rd = 1:num_radial_steps
        distance_km = radial_distances_km(rd);
        
        for az = 1:azimuth_steps
            % Convert azimuth from grid index to bearing (degrees clockwise from north)
            bearing = az - 1;  % 0 to 359 degrees
            
            if distance_km == 0
                % At center point
                xx_geo(az, rd) = cen_lon;
                yy_geo(az, rd) = cen_lat;
            else
                % Calculate destination point using geodesic calculations
                [lat2, lon2] = reckon(cen_lat, cen_lon, distance_km, bearing, wgs84);
                
                % Normalize longitude
                lon2 = mod(lon2, 360);
                
                xx_geo(az, rd) = lon2;
                yy_geo(az, rd) = lat2;
            end
        end
    end
    
    % Handle data near the longitude discontinuity (around 0/360°)
    % For interpolation, we might need to adjust points to be consistent
    
    % Prepare data for interpolation - handle longitude wrapping
    interp_lon = mod(sar_lon, 360);
    
    % For points that might cross the 0/360 boundary near the center,
    % we create a duplicated dataset shifted by 360 degrees
    center_lon = mod(cen_lon, 360);
    
    % If center is near 0°, we need to handle wrap-around
    if center_lon < 10 || center_lon > 350
        % Create duplicate points shifted by 360 degrees
        lon_shifted = interp_lon;
        lon_shifted(interp_lon > 180) = interp_lon(interp_lon > 180) - 360;
        lon_shifted(interp_lon <= 180) = interp_lon(interp_lon <= 180) + 360;
        
        % Combine original and shifted points
        all_lon = [interp_lon(:); lon_shifted(:)];
        all_lat = [sar_lat(:); sar_lat(:)];
        all_wsp = [wsp_sar(:); wsp_sar(:)];
    else
        all_lon = interp_lon(:);
        all_lat = sar_lat(:);
        all_wsp = wsp_sar(:);
    end
    
    % Remove any NaN values
    valid_idx = ~isnan(all_lon) & ~isnan(all_lat) & ~isnan(all_wsp);
    all_lon = all_lon(valid_idx);
    all_lat = all_lat(valid_idx);
    all_wsp = all_wsp(valid_idx);
    
    % Check if we have enough data points
    if length(all_lon) < 3
        error('Insufficient valid data points for interpolation (need at least 3)');
    end
    
    % Create scattered interpolant
    try
        F = scatteredInterpolant(double(all_lon), double(all_lat), double(all_wsp));
        F.Method = 'linear';  % 'linear' is usually better than 'nearest' for sparse data
        F.ExtrapolationMethod = 'none';
        
        % Interpolate to polar grid
        % For interpolation grid, we might need to adjust longitudes
        interp_grid_lon = xx_geo;
        
        % If center is near boundary, adjust interpolation grid too
        if center_lon < 10 || center_lon > 350
            % Create two versions of the grid
            grid_lon_1 = mod(xx_geo, 360);
            grid_lon_2 = grid_lon_1;
            grid_lon_2(grid_lon_1 > 180) = grid_lon_1(grid_lon_1 > 180) - 360;
            grid_lon_2(grid_lon_1 <= 180) = grid_lon_1(grid_lon_1 <= 180) + 360;
            
            % Try interpolation on both grids
            try
                wsp1 = F(grid_lon_1, yy_geo);
                valid1 = ~isnan(wsp1);
            catch
                wsp1 = NaN(size(grid_lon_1));
                valid1 = false(size(grid_lon_1));
            end
            
            try
                wsp2 = F(grid_lon_2, yy_geo);
                valid2 = ~isnan(wsp2);
            catch
                wsp2 = NaN(size(grid_lon_2));
                valid2 = false(size(grid_lon_2));
            end
            
            % Combine results
            wsp_sar_polar = NaN(size(xx_geo));
            wsp_sar_polar(valid1) = wsp1(valid1);
            wsp_sar_polar(valid2 & isnan(wsp_sar_polar)) = wsp2(valid2 & isnan(wsp_sar_polar));
        else
            % Simple case - center not near boundary
            wsp_sar_polar = F(xx_geo, yy_geo);
        end
        
    catch ME
        warning('Interpolation failed: %s', ME.message);
        wsp_sar_polar = NaN(size(xx_geo));
    end
    
  
end
