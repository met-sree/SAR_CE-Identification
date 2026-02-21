function process_tc_case2(case_name)
% PROCESS_TC_CASE Process TC case data and create MAT file
%   process_tc_case('WP272024_USAGI') processes the specified TC case

clc

% Display the case name being processed
fprintf('Processing case: %s\n', case_name);

addpath('/scratch/ynekkali/Data_Processesing/Paper_data/Functions/');
maxrds=300; %%% interpolation of radius upto (km)
res=1; %%%% resolution of interpolation (degree)
radial_dist_km= 0:res:maxrds;
%%%%read the file %%%%%%%%%%%%%%%%%%%%%%%%
% Set folder path (using the input case_name)
folder_path = strcat('/scratch/ynekkali/NOAA_STAR_RAW/Named_TC/',case_name,'/');
tc_name=case_name(5:end);
ibfilePath="/scratch/ynekkali/IBTrACS.ALL.v04r01.nc";
mat_folderPath="/scratch/ynekkali/Data_Processesing/Paper_data/Processed_Data/";
coverge_plot_location="/scratch/ynekkali/Data_Processesing/Paper_data/Plots/Plots_azim/";
quad_plot_location= "/scratch/ynekkali/Data_Processesing/Paper_data/Plots/Plots_quad/";

% Create directories if they don't exist
if ~exist(strcat(coverge_plot_location,tc_name), 'dir')
    mkdir(strcat(coverge_plot_location,tc_name))
end
if ~exist(strcat(quad_plot_location,tc_name), 'dir')
    mkdir(strcat(quad_plot_location,tc_name))
end

%%%%%%%%%%%%% extrct IBtracks Data %%%%%%%%%%%%%%%%%%%
name_ibtracs = ncread(ibfilePath, 'name');
sid_ibtracs = ncread(ibfilePath, 'sid');
sid_str = cellstr(sid_ibtracs');
name_str = cellstr(name_ibtracs');
ib_year = cellfun(@(x) x(1:min(4, length(x))), sid_str, 'UniformOutput', false);
tc_year_name = cellfun(@(x,y) [x '_' y], ib_year, name_str, 'UniformOutput', false);
[is_member, member_locations] = ismember(tc_year_name, tc_name);
matching_tc_indices = find(is_member);
IBTrACS=create_IBTrACKS(ibfilePath, matching_tc_indices);
valid_idx = ~isnan(IBTrACS.usa_lon) & ~isnan(IBTrACS.usa_lat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if folder exists
if ~exist(folder_path, 'dir')
    error('Folder not found: %s', folder_path);
end

% Get list of .txt files
txt_files = dir(fullfile(folder_path, '*.txt'));
txt_names = {txt_files.name};

% Initialize arrays (preallocate if possible)
if ~isempty(txt_names)
    % Initialize cell arrays and matrices
    num_files = length(txt_names);
    sar_lon_polar = [];
    sar_lat_polar = [];
    wsp_raw_polar=[];
    wsp_bufr_polar=[];
    wsp_qc_polar=[];
    sar_snap_cen_lat = [];
    sar_snap_cen_lon = [];
    sar_snap_vmax = [];
    sar_snap_rmax = [];
    sar_snap_time = NaT(num_files, 1);
    sar_data_coverage=[];
    sector_eyewall_info=[];
    sector_moat_info=[];
    sector_CE_flag =[];
    CE_flag = [];
    sar_vmax = [];
    sar_rmax = [];
    ib_lat = [];
    ib_lon = [];
    ib_rmax = [];
    ib_vmax = [];
    ib_sshs =[];
else
    warning('No text files found in: %s', folder_path);
    return;
end

% Loop through each .txt file
for i = 1:length(txt_names)
    fprintf('Processing file %d of %d: %s\n', i, length(txt_names), txt_names{i});
    
    % Get the base name without extension
    [~, base_name, ~] = fileparts(txt_names{i});
    dt = datetime(regexp(base_name, '\d+_\d+_\d+_\d+_\d+_\d+', 'match'), ...
        'InputFormat', 'yyyy_MM_dd_HH_mm_ss');
    %%%% extract lat and lon from SAR text file %%%%%%%%%%%%%%%
    [snap_cen_lat, snap_cen_lon, snap_date] = ExtractStromInfoFromTXTfile(folder_path, base_name);
    % Remove '_NOAAinfo_summary' from the base name
    base_name = strrep(base_name, '_NOAAinfo_summary', '');
    %%%%%% read the netcdf file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nc_file = strcat(folder_path,base_name,'_wind_level2.nc');
    
    % Check if netcdf file exists
    if ~exist(nc_file, 'file')
        warning('NetCDF file not found: %s', nc_file);
        continue;
    end
    
    sar_lat = ncread(nc_file, 'latitude');
    sar_lon=ncread(nc_file, 'longitude');
    wsp_sar = ncread(nc_file, 'sar_wind')*1.94384449;
    %%%% mask the data %%%%%%%%%%%%%%%%
    mask = ncread(nc_file, 'mask');
    wsp_sar_masked=wsp_sar;
    wsp_sar_masked(mask == 1) = NaN;
    [wsp_raw_polar1, xx, yy] = cart_to_polar(sar_lon, sar_lat, wsp_sar_masked, maxrds, res, snap_cen_lon, snap_cen_lat); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wsp_sar1=wsp_sar;
    buffer_distance_km1=1; %%% 1 km bufr line
    mask_with_buffer1=create_buffer_mask(mask, buffer_distance_km1, 0.5);
    wsp_sar1(mask_with_buffer1 >=0)=NaN;
    wsp_sar3=wsp_sar;
    buffer_distance_km3=3; %%% 1 km bufr line
    mask_with_buffer3=create_buffer_mask(mask, buffer_distance_km3, 0.5);
    wsp_sar3(mask_with_buffer3 >=0)=NaN;
    %%%%%%%%%%% convert cartsian to polar %%%%%%%%%%%%%
    [wsp_bufr1km_polar1, xx, yy] = cart_to_polar(sar_lon, sar_lat, wsp_sar1, maxrds, res, snap_cen_lon, snap_cen_lat);
    [wsp_bufr3km_polar1, xx, yy] = cart_to_polar(sar_lon, sar_lat, wsp_sar3, maxrds, res, snap_cen_lon, snap_cen_lat);
    %%%%%%%%%%%%% quality check %%%%%%%%%%%%%%%%%%%%%
    [wsp_qc_polar1, flag, stats] = qc_3sigma_polar(wsp_bufr3km_polar1);
    %%%%%%%%%%%%%%% calculate the coverge percentage %%%%%%%%%%%%%%
    [n_azimuth, n_radius] = size(wsp_qc_polar1);
    non_nan_percentage = 100 * (sum(~isnan(wsp_qc_polar1), 1) / n_azimuth);
    coverge_in_150km=nanmean(non_nan_percentage(1:150));
    %%%%%%%%%%%%%%%%% Averge in each 45 degree sector %%%%%%%%%%%%
    [fig2, sector_ce_flag1, eyewall_data, moat_data,sector_coverage_percentage1]=analyze_azimuthal_profiles(wsp_qc_polar1,radial_dist_km);
    print(fig2,strcat(quad_plot_location,tc_name,'/',tc_name,'_',string(i),'.jpeg'),'-djpeg')
    close(fig2)
    
    %%%%%%%%%%%%% Create Mat file %%%%%%%%%%%%%%%%%%%%%%
    azi_wsp=nanmean(wsp_qc_polar1, 1);
    vmax=max(azi_wsp(:));
    [id1,id2] = find(azi_wsp==vmax);

    % Store data
    if i == 1
        sar_wsp_polar = zeros(size(wsp_qc_polar1, 1), size(wsp_qc_polar1, 2), num_files);
        sar_lon_polar = zeros(size(xx, 1), size(xx, 2), num_files);
        sar_lat_polar = zeros(size(yy, 1), size(yy, 2), num_files);
        qud_eyewalls_info = zeros(size(eyewall_data, 1), size(eyewall_data, 2), num_files);
        qud_moat_info = zeros(size(moat_data, 1), size(moat_data, 2), num_files);
    end
    
    sar_lon_polar(:,:,i)=xx;
    sar_lat_polar(:,:,i)=yy;
    wsp_raw_polar(:,:,i)=wsp_raw_polar1;
    wsp_bufr1km_polar(:,:,i)=wsp_bufr1km_polar1;
    wsp_bufr3km_polar(:,:,i)=wsp_bufr3km_polar1;
    wsp_qc_polar(:,:,i)=wsp_qc_polar1;
    sar_snap_cen_lat(i,1)=snap_cen_lat;
    sar_snap_cen_lon(i,1)=snap_cen_lon;
    sar_snap_time(i,1)=snap_date;
    sector_eyewall_info(:,:,i)=eyewall_data;
    sector_moat_info(:,:,i)=moat_data;
    sector_coverage_percentage(:,:,i)=sector_coverage_percentage1;
    %%%%%%%%%%%%%%% for removing data with rain bands %%%%%%%%%%%%%%%%%%%%%%%%%
    %idx_rainbnd=find(abs(diff(eyewall_data(:,4)))>50);
    %condition_results(idx_rainbnd)=0;
    %qud_cond_info(:,i)=condition_results;
    sector_ce_flag(:,i)=sector_ce_flag1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% for checking at least 5 azimutal sections having CE %%%%%%%%%%%%%%%%%
    if sum(sector_ce_flag1) >= 5
        ce_flag(i,1) = 1;
    else
        ce_flag(i,1) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sar_data_coverage(:,i)=non_nan_percentage;

    %%%%%%%%%%%%%%%%% extract IBTrACS data for snap time %%%%%%%
    if any(valid_idx)
        ib_lat(i,1)= interp1(datenum(IBTrACS.datetimes(valid_idx)), IBTrACS.usa_lat(valid_idx), datenum(dt), 'linear');
        ib_lon(i,1)= interp1(datenum(IBTrACS.datetimes(valid_idx)), IBTrACS.usa_lon(valid_idx), datenum(dt), 'linear');
        ib_rmax(i,1)= interp1(datenum(IBTrACS.datetimes(valid_idx)), IBTrACS.usa_rmw(valid_idx), datenum(dt), 'linear');
        ib_vmax(i,1)= interp1(datenum(IBTrACS.datetimes(valid_idx)), IBTrACS.usa_wind(valid_idx), datenum(dt), 'linear');
	ib_sshs(i,1)= interp1(datenum(IBTrACS.datetimes(valid_idx)), IBTrACS.sshs(valid_idx), datenum(dt), 'nearest');
        ib_sshs_lmi=max(IBTrACS.sshs(valid_idx));
    else
        ib_lat(i,1)= NaN;
        ib_lon(i,1)= NaN;
        ib_rmax(i,1)= NaN;
        ib_vmax(i,1)= NaN;
	ib_sshs(i,1)=NaN;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%% plotting radial proilfe and data coverge %%%%%
    [fig1, ax]=plot_wind_profile_with_coverage(radial_dist_km, azi_wsp, non_nan_percentage, ce_flag(i,1));
    print(fig1,strcat(coverge_plot_location,tc_name,'/',tc_name,'_',string(i),'.jpeg'),'-djpeg')
    close(fig1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

end

% Define variable names with descriptions
var_data = {
    'sar_snap_time', 'SAR Time';
    'sar_lon_polar', 'Longitude coordinates in polar co-ordinate';
    'sar_lat_polar', 'Latitude coordinates in polar co-ordinate';
    'wsp_raw_polar', 'Interpolated SAR wind speed in polar co-ordinate';
    'wsp_bufr1km_polar','Interpolated SAR wind after 1km bufr along coastline'
    'wsp_bufr3km_polar', 'Interpolated SAR wind after 3 km bufr along coastline'
    'wsp_qc_polar',' Interpolated SAR wind after bufr and 3 sigm quality check'
    'sar_snap_cen_lat', 'Center latitudes from text file';
    'sar_snap_cen_lon', 'Center longitudes from text file';
    'sector_eyewall_info', 'Each azimuthal 45 degree two eyewalls wind speed and its radius';
    'sector_moat_info', 'Each azimuthal 45 degree moat wind speed and its radius';
    'sector_coverage_percentage','Sectorwise data coverge percentage'
    'sector_ce_flag', 'Each azimuthal 45 degree CE information (0 is no CE and 1 for CE)';
    'ce_flag', 'Overall CE for snap (0 is no CE and 1 for CE)';
    'sar_data_coverage', 'Coverage percentage';
    'ib_lat', 'IBTrACS latitude at snap time';
    'ib_lon', 'IBTrACS longitude at snap time';
    'ib_rmax', 'IBTrACS rmax at snap time';
    'ib_vmax', 'IBTrACS vmax at snap time';
    'ib_sshs', 'IBTrACS intensity at snap time';
    'ib_sshs_lmi', 'IBTrACS maximum intensity of TC';
    };

% Create STAR structure
STAR = struct();

% Assign each variable
for i = 1:size(var_data, 1)
    var_name = var_data{i, 1};
    if exist(var_name, 'var')
        STAR.(var_name) = eval(var_name);
    else
        warning('Variable %s not found', var_name);
    end
end

%%%%%%%%%%%merge two mat file and create Mat file %%%%%%%%%%%%%%%
mat_file_path = fullfile(mat_folderPath, [tc_name '.mat']);
save(mat_file_path, 'IBTrACS', 'STAR', '-v7.3');
fprintf('MAT file saved: %s\n', mat_file_path);

end
