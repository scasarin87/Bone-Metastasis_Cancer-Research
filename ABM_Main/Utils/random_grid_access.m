%% Random Grid Access

%  This function is used to define randomly, at every iteration, the order
%  of the ABM site to be investigated adding stochasticity.

%  Input  -> rows, columns     : Y and X size of the grid
%
%  Output -> random_grid_order : rows*columns vector containing the site 
%                                indexes in order of investigation

function [random_grid_order] = random_grid_access(rows, columns)

    rand_order = rand(rows * columns, 1);
    [~, random_grid_order] = sort(rand_order);

end