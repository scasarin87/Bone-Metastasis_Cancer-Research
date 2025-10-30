%% Load Bone Image from Directory

% Input : filename -> name of the image to be loaded
% Output: image    -> the variable containing the bone image

function [image] = load_BoneImage(filename)
     
     SetDirectoriesGeometry;
     BoneImageDirectory = strcat(imageDirectory, filename);
     image = imread(filename);

end

