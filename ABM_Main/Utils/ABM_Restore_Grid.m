%% ABM Restore Grid

%  This function is used to restore the cortical bone and osteoblast ABM 
%  sites if some computation lead to a previous change

% Filling the grid
for i = 1:rows
    for j = 1:columns

        if cortical_bone(i, j) == 1 
            bone(i, j) = site.cortical_bone;
        end 

        if (osteoblasts(i, j) == 1 && bone(i, j) == site.cortical_bone)
            bone(i, j) = site.osteoblast;
        end 

    end 
end   

clear i j