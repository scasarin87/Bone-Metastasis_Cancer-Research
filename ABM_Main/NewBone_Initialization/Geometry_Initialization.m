%% ABM Geometry Initialization

% Loading and precessing of bone images + cortical bone, bone marrow tumor and 
% vessels mask creation: starting from the bone images returns all the masks 
%(ex: cortical_bone_mask = 1 if site belonging to cortical bone, 0 otherwise). 
% Finally, it merges the masks and intializes the ABM.

clear all

bone_image_filename = 'tibia_ls_1.jpg'; % Path to be specified in SetDirectories function
nameNewFolder       = "tibia_ls_1";% Folder where I'll save all the output masks

flag_save   = 0; % = 1 if you want to save the created masks in nameNewFolder
x_real_size = 2400; % in questa variabile vado a mettere la size (in um) dell'asse x dell'immagine che mi da eleonora, quindi quanti um rappresenta sull'asse x l'immagine che mi da

parameters; % Fixed parameters

% Load image
bone_image = load_BoneImage(bone_image_filename);

% Processing Image: binarization, reshape, padding
[processed_image, rows, columns, X, Y] = processing_BoneImage(bone_image, x_real_size, site_dim);

% Grid Features Definition
[ax, ay, bx, by, Cx, Cy, directionx, directiony] = grid_FeaturesDefinition(X, Y);

% Mask Creation
[cortical_bone_mask, outer_mask, ...
    bone_marrow_mask] = MaskDefinition(processed_image, X, Y, site_dim, T_cell_division);

% [cortical_bone_mask, outer_mask, ...
%  bone_marrow_mask] = MaskDefinition_Gerardo(processed_image);

% Crop matrices
[cortical_bone_mask, bone_marrow_mask, outer_mask, rows, columns] = crop_masks(cortical_bone_mask, bone_marrow_mask, outer_mask);                              

% Bone matrix initialization
[bone] = MaskMerging(cortical_bone_mask, outer_mask, rows, columns);

% Show ABM
[ABM_matrix] = display_abm(rows, columns, bone, 1);
DisplayMasks;

% Manual Refinement Step
% AfterManualRefinement;

% Save Masks Matrices in nameNewFolder
if flag_save == 1
   SaveVariables;
end


