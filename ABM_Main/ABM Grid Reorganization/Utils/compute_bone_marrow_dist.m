%% Compute Bone Marrow Distance %%

%  This function is used to compute the distance between 

%  Input  -> bone             : current bone matrix
%            bone_marrow_dist : matrix containing the distance between
%                               current tumor site and the bone marrow
%                               edges
%            row/col_agent    : row/col of the current PCa site considered
%            site             : struct defining the ABM sites
%            X, Y             : hexagonal grid 

%  Output -> bone_marrow_dist   : matrix containing the distance between
%                                current tumor site and the bone marrow
%                                edges

function [bone_marrow_dist] = compute_bone_marrow_dist(bone, bone_marrow_dist, ...
                                 row_agent, col_agent, rows, columns, X, Y, site)
                             
    % Scan all the grid to fill the 'bone_marrow_dist' matrix
    for row_marrow = 2 : rows - 1
        for col_marrow = 2 : columns - 1                             

            % We are interested only in the bone marrow site adjacent to the tumor (da riscrivere questo if, che fa schifo) 
            if bone(row_marrow, col_marrow) == site.bone_marrow && ((bone(row_marrow + 1, col_marrow) == site.tumor || ...
                        bone(row_marrow - 1, col_marrow) == site.tumor || ...
                            bone(row_marrow, col_marrow + 1) == site.tumor || ...
                                bone(row_marrow, col_marrow - 1) == site.tumor || ...
                                    bone(row_marrow + 1, col_marrow + 1) == site.tumor || ...
                                         bone(row_marrow + 1, col_marrow - 1) == site.tumor || ...
                                            bone(row_marrow - 1, col_marrow + 1) == site.tumor || ...
                                                bone(row_marrow - 1, col_marrow - 1) == site.tumor) || (bone(row_marrow + 1, col_marrow) == site.tumor_edge || bone(row_marrow - 1, col_marrow) == site.tumor_edge || bone(row_marrow, col_marrow + 1) == site.tumor_edge || bone(row_marrow, col_marrow - 1) == site.tumor_edge || bone(row_marrow + 1, col_marrow + 1) == site.tumor_edge || bone(row_marrow + 1, col_marrow - 1) == site.tumor_edge || bone(row_marrow - 1, col_marrow + 1) == site.tumor_edge || bone(row_marrow - 1, col_marrow - 1) == site.tumor_edge))

               % Get the tumor-bone marrow distance
               bone_marrow_dist(row_marrow, col_marrow) = compute_distance(X, Y, row_agent, col_agent, row_marrow, col_marrow);
           
            end 
        end
    end

    clear row_marrow col_marrow
    
end