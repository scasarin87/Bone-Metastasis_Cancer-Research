%% Mask Merging
%  This function is used to merge all the masks previously defined and thus
%  initilizing the ABM
%
% Input  -> cortical_bone_mask : = 1 if cortical bone site
%           outer_mask         : = 1 if outer site
%
% Output -> bone matrix        : -2 -> outer site
%                                -1 -> cortical bone
%                                 0 -> bone marrow
%                                 2 -> tumor site
%                                 3 -> vessel site
%                                 4 -> near bone tumor site
%                                 5 -> Osteoblast site
%                                 6 -> Osteoclast site

function [bone] = MaskMerging(cortical_bone_mask, outer_mask, rows, columns)
    
    bone = zeros(rows, columns); % Initialize bone matrix to 0
    
    OUTER_SITE          = -2;
    CORTICAL_SITE       = -1;
    BONEMARROW_SITE     =  0;
    TUMOR_SITE          =  2;
    VESSEL_SITE         =  3;
    TUMOR_NEARBONE_SITE =  4;
    OSTEOBLAST_SITE     =  5;
    OSTEOCLAST_SITE     =  6;
    
    for i = 1:rows
        for j = 1:columns
            
            if outer_mask(i, j) == 1  
                bone(i, j) = OUTER_SITE;
            end
            
            if cortical_bone_mask(i, j) == 1
                bone(i, j) = CORTICAL_SITE;
            end  
                        
        end 
    end 
    
end 