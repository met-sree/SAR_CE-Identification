function mask_with_buffer = create_buffer_mask(mask, buffer_distance_km, resolution_km)
% CREATE_BUFFER_MASK Create a buffer zone around land mask
%
% Syntax:
%   mask_with_buffer = create_buffer_mask(mask, buffer_distance_km, resolution_km)
%
% Inputs:
%   mask - Original mask matrix (typically with values: -1=ocean, 0=buffer, 1=land)
%   buffer_distance_km - Desired buffer distance in kilometers
%   resolution_km - Grid resolution in kilometers
%
% Output:
%   mask_with_buffer - Updated mask with buffer zone (values: -1=ocean, 0=buffer, 1=land)
%
% Example:
%   mask_buf = create_buffer_mask(mask, 1, 0.5);

    % Calculate buffer size in pixels
    buffer_pixels = round(buffer_distance_km / resolution_km);
    
    fprintf('Creating %.1f km buffer zone (%.0f pixels)\n', buffer_distance_km, buffer_pixels);
    
    % Create land mask (where mask == 1)
    land_mask = mask == 1;
    
    % Fast buffer creation using convolution
    kernel_size = 2 * buffer_pixels + 1;
    kernel = ones(kernel_size);
    
    % Create buffer zone using convolution
    convolved = conv2(double(land_mask), kernel, 'same');
    buffer_zone = convolved > 0 & ~land_mask;
    
    % Create new mask with 3 classes
    mask_with_buffer = mask;  % Start with original
    mask_with_buffer(buffer_zone) = 0;  % Buffer zone = 0
    
    % Display summary
    fprintf('Buffer creation complete:\n');
    fprintf('  Ocean pixels: %d\n', sum(mask_with_buffer(:) == -1));
    fprintf('  Buffer pixels: %d\n', sum(mask_with_buffer(:) == 0));
    fprintf('  Land pixels: %d\n', sum(mask_with_buffer(:) == 1));
end