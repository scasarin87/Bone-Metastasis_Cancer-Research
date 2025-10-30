%% Bone Coordinates Generator %%

%  This function is used to extract a txt file with all the cortical bone
%  coordinates needed in Abaqus

% Extract cortical bone coordinates
[cb_rows, cb_cols] = find(bone == site.cortical_bone);
% Define nodes                                
nodes = [cb_cols cb_rows];                                 
% Set the output file name
filename_format = 'cb_coordinates_%.2d_cx_%d_cy_%d_d_%d.txt';
filename = sprintf(filename_format, epoch, CxT, CyT, tumor_dim);
% Define folder name
folder_name = 'CoordinatesExtraction';
% Txt file generation function
Matlab2Abaqus(nodes, filename);

% save('Results/final_tumor_area.mat', 'final_tumor_area');
% save('Results/resorption_area.mat', 'resorption_area');


