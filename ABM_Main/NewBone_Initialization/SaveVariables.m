%% Save current geometry variables

SetDirectoriesGeometry;

cd(outputDirectory)      % Move to geometry directory
mkdir(nameNewFolder)    % Create the directory according to the given name
cd(nameNewFolder)       % Enter the new folder

% Save geometries
save("cortical_bone", "cortical_bone_mask") 
save("bone_marrow"  , "bone_marrow_mask")
save("outer"        , "outer_mask")

cd(codeDirectory)          % Back to the code directory
