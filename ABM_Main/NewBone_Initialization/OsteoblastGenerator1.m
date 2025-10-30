%% Osteoblast Generator %%

%  This function computes the internal edge of the cortical bone and then
%  it generates an OBs matrix according to the threshold (rand_value)

% Input  -> rows, columns  : x and y size of the image
%           OBs_percentage : % of cortical bone internal edge to be covered
%                            by OBs
%           cort_bone_mask : = 1 if site belongs to cortical bone
%           outer_mask     : = 1 if site belongs to background
%           tumor          : = 1 if site belongs to tumor cell
%
% Output -> OsteoblastMask : = 1 if site belongs to OBs


function [OsteoblastMask] = OsteoblastGenerator1(rows, columns, OBs_percentage, cortical_bone_mask, bone_marrow_mask, tumor)
   
    OsteoblastMask = zeros(rows, columns); % OBs mask initialization
    
    for i_mask = 2:rows-1
        for j_mask = 2:columns-1

            if cortical_bone_mask(i_mask,j_mask) == 1

                % Ugly way to define the edge of the cortical bone
                if (bone_marrow_mask(i_mask+1,j_mask) == 1 || bone_marrow_mask(i_mask-1,j_mask) == 1 ||  bone_marrow_mask(i_mask,j_mask-1) == 1 ||  bone_marrow_mask(i_mask,j_mask + 1) == 1)

                    if  rand(1) <= OBs_percentage

                        OsteoblastMask(i_mask,j_mask) = 1;

                    end 
                end 
            end 
        end 
    end 
    
end