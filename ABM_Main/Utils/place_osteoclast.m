%% Place Osteoclast %%

%  This function generate a single OCs mask and place it within the bone
%  matrix

function [bone] = place_osteoclast(bone, row_cortbone_edge, col_cortbone_edge, ...
                                            curr_cb_site, method, cortical_bone_edge, site)
    
    % Here we increase the size of the single OC, so that it will occupy 3 ABM sites at least (60 um)
    osteoclast_matrix = ones(size(bone, 1), size(bone, 2));
    osteoclast_matrix(row_cortbone_edge(curr_cb_site), col_cortbone_edge(curr_cb_site)) = 0;
    osteoclast_matrix = imerode(osteoclast_matrix, method);
    osteoclast_matrix = not(osteoclast_matrix);
    osteoclast_matrix = (osteoclast_matrix & cortical_bone_edge);
    [row_ocs, col_ocs] = find(osteoclast_matrix);
    % Here I Place the OC cell 
    for oc_site = 1 : size(row_ocs, 1)
        bone(row_ocs(oc_site), col_ocs(oc_site)) = site.osteoclast;
    end       


end