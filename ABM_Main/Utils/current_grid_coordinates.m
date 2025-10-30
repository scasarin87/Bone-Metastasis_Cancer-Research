%% Current Grid Coordinates %%

%  This function is used to define the current grid coordinates; remember
%  that we are accessing randomly to the ABM grid to improve model
%  stochasticity throughout simulations. The order is defined by the
%  'random_grid_order' vector.

%  Input  -> rows, columns        : Y and X size of the grid
%            agent                : current agent index
%            random_grid_order    : rows*columns vector containing the site 
%                                   indexes in order of investigation
%
%  Output -> row_agent, col_agent : Y and X coordinates of the current
%  sited

function [row_agent, col_agent] = current_grid_coordinates(rows, columns, ...
                                                              agent, random_grid_order)
                                                          
    % row_agent = ceil(random_grid_order(agent) / rows);
    % col_agent = random_grid_order(agent) - (row_agent - 1) * rows;
    row_agent = ceil(random_grid_order(agent) / columns);
    col_agent = random_grid_order(agent) - (row_agent - 1) * columns;
    
end