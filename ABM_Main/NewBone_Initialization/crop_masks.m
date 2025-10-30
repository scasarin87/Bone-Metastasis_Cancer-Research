%% Crop Masks
%  This function is used to crop all the masks previously defined
%  respecting the cortical bone mask size and taking a certain delta value
%
% Input  -> cortical_bone_mask : = 1 if cortical bone site
%           outer_mask         : = 1 if outer site
%
% Output -> cortical_bone_mask : cropped mask
%           outer_mask         : cropped mask

function [cortical_bone_mask, bone_marrow_mask, outer_mask, rows, columns] = crop_masks(cortical_bone_mask, bone_marrow_mask, outer_mask)
    
    % Sites to keep as margin
    delta_sites = 5; 
    % Find minimum and maximum values for cortical bone mask row and col
    [row, col] = find(cortical_bone_mask == 1);
    min_row = max(1, min(row) - delta_sites);
    max_row = min(size(cortical_bone_mask, 1), max(row) + delta_sites);
    min_col = max(1, min(col) - delta_sites);
    max_col = min(size(cortical_bone_mask, 2), max(col) + delta_sites);
    % Crop masks accordingly
    cortical_bone_mask = cortical_bone_mask(min_row : max_row, min_col : max_col);
    bone_marrow_mask = bone_marrow_mask(min_row : max_row, min_col : max_col);
    outer_mask = outer_mask(min_row : max_row, min_col : max_col);
    % Define new row, col size
    [rows, columns] = size(cortical_bone_mask);
    
end