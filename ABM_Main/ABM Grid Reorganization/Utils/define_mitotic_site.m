%% Define Mitotic Site %%

%  This function is used to define the bone marrow site adjacent to the
%  tumor that is either becoming a PCa cell during mitosis or that defines
%  the closest PCa cell that will undergo apoptosis. We subdivide tumor
%  into 4 quadrants, if a PCa undergoes mitosis, the new cell will be
%  placed in the bone marrow site adjacent to the tumor in the
%  corresponding quadrant.

%  Input  -> bone           : current bone matrix
%            X, Y           : hexagonal grid
%            site           : struct defining the ABM sites
%            row/col_agent  : current PCa cell row and col 

%  Output -> row/col_mitosis: mitosis site

function [row_mitosis, col_mitosis] = define_mitotic_site(bone, site, X, Y, ...
                                row_agent, col_agent, Rad, hour, start_rad, end_time, flag_cabo)    
    
        % Define Whole Tumor and Complementary Tumor Mask Matrix
        [curr_tumor, no_tumor] = def_tumor_masks(bone, site);
        
        % Define Properties
        % def_current_region(curr_tumor, row_agent, col_agent);

        % Define Current Tumor External Boundary 
        [marrow_edge] = def_tumor_external_boundary(no_tumor, bone, site);

        % Bone Marrow Site Selection
        [row_mitosis, col_mitosis] = bm_site_selection(marrow_edge, row_agent, col_agent, X, Y, Rad, hour, start_rad, end_time, flag_cabo);      

end