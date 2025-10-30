%% Display ABM %%

%  This function displays the ABM by identifying each model site with a color
%  according to value assumed by the "bone" matrix

% Input  -> rows, columns  : x and y size of the image
%           bone           : main ABM matrix, each site is labeled according 
%                            to the previous mask merging function
%           show_abm       : if = 1 the abm is shown, otherwise not
%
% Output -> ABM_matrix     : We now assign a color to each bone matrix site
%                            and store it in this matrix 

function [ABM_matrix] = display_abm(rows, columns, bone, show_abm)

    global current_geometry
    
    if show_abm
    
        % Color Struct Definition
        DefineColors;
    
        % ABM Sites Struct Definition
        DefineABM_sites;
    
        channels   = 3; % RGB
        ABM_matrix = zeros(rows, columns, channels);
    
        for i = 1:rows
            for j = 1:columns
    
                switch bone(i, j)
                    case site.bone_marrow
                        ABM_matrix(i, j, :) = color.light_grey;
    
                    case site.cortical_bone
                        ABM_matrix(i, j, :) = color.dark_grey;

                    case site.cortical_bone_induced
                        ABM_matrix(i, j, :) = color.dark_grey;

                    case site.outer
                        ABM_matrix(i, j, :) = color.white;
    
                    case {site.tumor_edge}
                        ABM_matrix(i, j, :) = color.light_pink;
    
                    case {site.tumor}
                        ABM_matrix(i, j, :) = color.mill_pink;
    
                    case site.vessel
                        ABM_matrix(i, j, :) = color.red;
    
                    case site.vessel_cabo
                        ABM_matrix(i, j, :) = color.purple;
    
                    case site.osteoblast
                        ABM_matrix(i, j, :) = color.yellow;
    
                    case site.osteoclast
                        ABM_matrix(i, j, :) = color.ocs_blue;
    
                    case site.mitotic_cells
                        ABM_matrix(i, j, :) = color.green;
    
                    case site.apoptotic_cells
                        ABM_matrix(i, j, :) = color.orange;
    
                    otherwise('Unexpected ABM site type! Check bone matrix');
    
                end
            end
        end
    
        fig = gcf;
        if (strcmp(current_geometry, 'femur_ls_1') == 1) || (strcmp(current_geometry, 'calvaria_ls_1') == 1)
            set(fig, 'Units', 'pixels', 'Position', [450, 50, 400, 710]);
        elseif (strcmp(current_geometry, 'tibia_ls_1') == 1) || (strcmp(current_geometry, 'tibia_ls_1_lr') == 1)
            set(fig, 'Units', 'pixels', 'Position', [450, 50, 400, 710]);
        elseif (strcmp(current_geometry, 'tibia_cs_1') == 1) || (strcmp(current_geometry, 'tibia_cs_2') == 1) || (strcmp(current_geometry, 'tibia_cs_3') == 1) || (strcmp(current_geometry, 'vertebra_cs_1' == 1))
            set(fig, 'Units', 'pixels', 'Position', [450, 80, 600, 500]);
        end
    
        pause(0.01)
    
        imagesc(ABM_matrix)
        axis off
    
    else
        ABM_matrix = 0; % Returns void matrix
    end
    
    end