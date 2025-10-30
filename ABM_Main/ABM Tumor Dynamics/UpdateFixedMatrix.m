%% Update Fixed Matrix %%

%  This function updates the matrices declared in ABM_Initialization.m
%  regarding the mitotic/apoptotic activity of the PCa cells

% Update the matrix that keep track of the mitotic/apoptotic events
CHANGE_CELL(:, : ,hour) = change_cell;
    
% Also Update Mitosis and Apoptosis Events
cells_matrix(2, hour) = length(find(change_cell ==  1)); % Mitotic events at this iteration
cells_matrix(3, hour) = length(find(change_cell == -1)); % Apoptotic events at this iteration
cells_matrix(4, hour) = sum(sum(internal_clock == T_cell_division)) - ...
                            (cells_matrix(2, hour) + cells_matrix(3, hour)); % Quiescent cells at this iteration