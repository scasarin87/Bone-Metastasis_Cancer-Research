%% ABM Initialization %%

%  This function initializes the ABM, defines its most important variables
%  and loads the correct bone geometry according to the user selction

% Fixed Parameters Definition + RAD dynamics parameters + Cabo Flag
FixedParameters; 

% Load the Desired Bone Geometry
[cortical_bone, bone_marrow, rows, columns] = load_geometry();

% Compute initial bone marrow area
bone_marrow_start_area = sum(bone_marrow, 'all');

% ABM Hexagonal Grid Building
[X, Y, ax, ay, bx, by, ...
    central_col, central_row, directionx, directiony] = hexagonal_grid(rows, columns);

% Tumor Mask Building
[tumor, vess_prob, internal_clock, a_tum, b_tum] = set_tumor(rows, columns, X, Y, T_cell_division);

% Vessels Mask Building
% C42B vessels
if strcmp(cell_line, 'c42b')
    [vessels, center_vessels, vessels_number, single_vessel_masks] = vessels_mask_c42b(rows, columns, ...
                                                                                  tumor, vess_prob, ...
                                                                                  site_dim, X, Y);
% PC3 vessels
elseif strcmp(cell_line, 'pc3')
    [vessels, single_vessel_masks, center_vessels, vessels_number] = vessels_mask_pc3(rows, columns, ...
                                                                                      tumor, site_dim, ...
                                                                                      X, Y);

% PC3 vessels
elseif strcmp(cell_line, 'renca')
    [vessels, single_vessel_masks, center_vessels, vessels_number, vess_degree_range] = vessels_mask_renca(rows, columns, ...
                                                                                                        tumor, site_dim, ...
                                                                                                        X, Y);
    vess_degree_range = [121 + 15, 181 - 15];
end

% ABM building ('bone' matrix set)        
[bone, site, tumor_area_start] = set_abm(cortical_bone, bone_marrow, ...
                                             tumor, vessels, rows, columns);

% Add Osteoblast and Osteoclast
[bone, BONE, MITOTIC, APOPTOTIC, ...
 papo_cell_near_obs, obs_influenced_cells, ...
 obs_number_start, cortical_bone_edge, ...
 no_bone_marrow, method, internal_clock_ocs] = osteoblast_osteoclast_generator(bone, site, ...
                                                                               X, Y, rows, ...
                                                                               columns, site_dim, ...
                                                                               T_ocs_resorption);                                                                            

% Show Initialized ABM
[abm_matrix] = display_abm(rows, columns, bone, show_abm);

% Fixed Matrices Definition (-> updated in ABM_Tumor_Dynamics.m)
FixedMatrix;

% If cabo is administrated we track just 15 days follow up with pc3
if strcmp(cell_line, 'pc3') && flag_cabozantinib == 1
    follow_up = 360;
    enter_condition = 1; % Needed in pc3_boundary_vessel_dynamics to set cabo vessel ratio
end



