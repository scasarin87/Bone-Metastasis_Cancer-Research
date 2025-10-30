%% Set ABM %%

%  This function is used to merge all the masks previously defined and thus
%  initilizing the ABM
%
% Input  -> cortical_bone : = 1 if cortical bone site
%           bone_marrow   : = 1 if bone marrow site
%           osteoblasts   : = 1 if osteoblast site
%           tumor         : = 1 if tumor site
%           vessels       : = 1 if vessel site
%           rows, columns : size of ABM grid
%
% Output -> bone          : -2 -> outer site
%                           -1 -> cortical bone
%                            0 -> bone marrow
%                            2 -> tumor site
%                            3 -> vessel site
%                            4 -> near bone tumor site
%                            5 -> Osteoblast site
%                            6 -> Osteoclast site
%           BONE          : records bone matrix over follow-up time
%           site          : struct containing ABM sites names
%           tumor_area_st : starting tumor area 

function [bone, site, tumor_area_start] = set_abm(cortical_bone, bone_marrow, ...
                                                         tumor, vessels, rows, columns)
  
    % ABM Sites Struct Definition
    DefineABM_sites;
    
    % Initialize bone matrix, each grid value = -2 (outer)
    bone = site.outer * ones(rows, columns);
    
    % Filling the grid
    for i = 1:rows
        for j = 1:columns
            
            if (cortical_bone(i, j) == 1 && bone(i, j) == site.outer)
                bone(i, j) = site.cortical_bone;
            end 
            
            if (bone_marrow(i, j) == 1 && bone(i, j) == site.outer)
                bone(i, j) = site.bone_marrow;
            end
            
            if (tumor(i, j) == 1 && bone(i, j) == site.bone_marrow)
                bone(i, j) = site.tumor;
            end
            
            if (vessels(i, j) == 1 && bone(i, j) == site.tumor)
                bone(i, j) = site.vessel;
            end       
            
        end 
    end   
    
    clear i j
    
    % Define Starting Tumor Area
    tumor_area_start = sum(sum(bone == site.tumor | bone == site.vessel));
                      
end
