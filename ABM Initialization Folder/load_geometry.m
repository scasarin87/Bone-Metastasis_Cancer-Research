%% Load Geometry %%

%  This function scans the geometries folder in my local environment, and
%  loads the cortical bone, bone marrow and osteoblast masks.

%  Input  -> current_geometry: string which contains the name of the
%                              desired geometry to be implemented
%
%  Output -> cortical_bone   : cortical bone mask -> = 1 if site belongs to CB
%            bone_marrow     : bone marrow mask   -> = 1 if site belongs to BM
%            bone_marrow_edge: = 1 if site is within the OBs influence to
%                              PCa cells when cabo therapy starts
%            osteoblasts     : osteoblasts mask   -> = 1 if site is OBs   
%            rows, columns   : size of the ABM grid

function [cortical_bone, bone_marrow, rows, columns] = load_geometry()
    
    global current_geometry;
    
    SetDirectories;
    
    cd(GeometryDirectory); % Defined in SetDirectories.mat
    cd(current_geometry);  % Directory that stores masks matrices
    
    load('cortical_bone')
    cortical_bone = cortical_bone_mask;
    
    load('bone_marrow.mat')
    bone_marrow = bone_marrow_mask;
    
    bound = 0;
    cortical_bone = padarray(cortical_bone, [bound, bound], 0, 'pre');
    cortical_bone = padarray(cortical_bone, [bound, bound], 0, 'post');
    
    bone_marrow = padarray(bone_marrow, [bound, bound], 0, 'pre');
    bone_marrow = padarray(bone_marrow, [bound, bound], 0, 'post');
    
    clear cortical_bone_mask bone_marrow_mask 
    
    [rows, columns] = size(cortical_bone);
    
    cd(codeDirectory)

end