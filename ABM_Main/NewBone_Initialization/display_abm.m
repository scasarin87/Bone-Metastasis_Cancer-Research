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

                    case site.outer 
                        ABM_matrix(i, j, :) = color.white;
                    
                    case {site.tumor_edge}
                        ABM_matrix(i, j, :) = color.light_blue;

                    case {site.tumor}
                        ABM_matrix(i, j, :) = color.blue;                        

                    case site.vessel
                        ABM_matrix(i, j, :) = color.red;
                    
                    case site.vessel_cabo
                        ABM_matrix(i, j, :) = color.purple;

                    case site.osteoblast
                        ABM_matrix(i, j, :) = color.yellow;

                    case site.osteoclast
                        ABM_matrix(i, j, :) = color.black;

                    otherwise('Unexpected ABM site type! Check bone matrix');

                end  
            end  
        end  

        image(ABM_matrix);
        axis off
        pause(0.01)
        
    else
        ABM_matrix = 0; % Returns void matrix 
    end 

end