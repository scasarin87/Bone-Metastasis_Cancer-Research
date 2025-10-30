%% Mask Definition
%  This function creates all the indipendent masks that individually
%  compose the final ABM structure.

% Input  -> processed_image: takes the image previously binarized, reshape,
%                            padded and returns
%           X, Y           : hexagonal grid sites
%           T_cell_division: time-span between consecutive cellular events
%
% Output -> cortical_bone_mask: 1 -> cortical bone site
%                               0 -> bone marrow and outer sites
%                                                                 
%           outer_mask        : 1 -> sites surrounding the cortical bone
%                               0 -> cortical bone and bone marrow
%
%           bone_marrow_mask  : 1 -> site belonging to bone marrow
%                               0 -> cortical bone and outer

function [cortical_bone_mask, outer_mask, bone_marrow_mask] = MaskDefinition(processed_image) 
   
    
    % 1) CORTICAL BONE MASK
    cortical_bone_mask = not(processed_image);
    conn = conndef(2,'maximal'); % pixels connectivity definition
    cortical_bone_mask = imfill(cortical_bone_mask, conn, 'holes'); % fill holes within the cortical bone
    
    % 2) OUTER MASK
    outer_mask = not(cortical_bone_mask);
    
    % 3) BONE MARROW MASK
    cortical_bone_mask = not(processed_image);
    bone_marrow_mask = not(cortical_bone_mask | outer_mask);
    
end



% OLD CODE: it doesn't work properly
%
% function [cortical_bone_mask, outer_mask, bone_marrow_mask] = MaskDefinition(processed_image, X, Y, site_dim, T_cell_division)
%
%     % 1) CORTICAL BONE MASK
%     cortical_bone_mask = not(processed_image);
%     conn = conndef(2,'maximal'); % pixels connectivity definition
%     % cortical_bone_mask = imfill(cortical_bone_mask, conn, 'holes'); % to remove holes within the cortical bone
% 
%     % 2) OUTER MASK (a smarter algorithm must be found)
%     rows       = size(cortical_bone_mask, 1);
%     columns    = size(cortical_bone_mask, 1);
%     outer_mask = zeros(rows, columns);
% 
%     % Check from left to right until I find the cortical edge
%     for i_row = 1:rows
%         for i_col = 1:columns
%             if cortical_bone_mask(i_row,i_col) == 1
%                 break
%             end
%             if cortical_bone_mask(i_row, i_col) == 0
%                 outer_mask(i_row, i_col) = 1;
%             end            
%         end 
%     end 
% 
%     % Check from right to left until I find the cortical edge
%     for i_row = 1:rows
%         for i_col = 1:columns
%             if cortical_bone_mask(rows + 1 - i_row, columns + 1 - i_col) == 1
%                 break
%             end
%             if cortical_bone_mask(rows + 1 - i_row, columns + 1 - i_col) == 0
%                 outer_mask(rows + 1 - i_row, columns + 1 - i_col) = 1;
%             end            
%         end 
%     end 
% 
%     % Check from up to down until I find the cortical edge
%     for i_col = 1:columns
%         for i_row = 1:rows
%             if cortical_bone_mask(i_row,i_col) == 1
%                 break
%             end
%             if cortical_bone_mask(i_row, i_col) == 0
%                 outer_mask(i_row, i_col) = 1;
%             end            
%         end 
%     end 
% 
%     % Check from down to up until I find the cortical edge
%     for i_col = 1:columns
%         for i_row = 1:rows
%             if cortical_bone_mask(rows + 1 - i_row, columns + 1 - i_col) == 1
%                 break
%             end
%             if cortical_bone_mask(rows + 1 - i_row, columns + 1 - i_col) == 0
%                 outer_mask(rows + 1 - i_row, columns + 1 - i_col) = 1;
%             end            
%         end 
%     end 
% 
%     clear i_col i_row
% 
%     % BONE MARROW MASK
%     bone_marrow_mask = zeros(rows, columns);
%     % find the coordinates of the combined ones between the two masks
%     marrow = union(find(cortical_bone_mask == 1), find(outer_mask == 1));
%     % dedine bone marrow mask
%     bone_marrow_mask(marrow) = 1; 
%     bone_marrow_mask = not(bone_marrow_mask); % need to invert 1 and 0 to keep consistency with other masks
% end


