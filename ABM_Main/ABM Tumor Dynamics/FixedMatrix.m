%% Fixed Matrix %%

%  This function initializes some fixed matrices that I'll fill during the
%  ABM simulations

%% 1) Follow-Up Data Matrices

% Initialize at every cycle the matrix that keeps track of the cellular activity

% This matrix can assume the following values:
% 0: if the tumor cell in rows, columns is quiescent at hour == follow_up
% 1: if the tumor cell in rows, columns undergo mitosis at hour == follow_up
%-1: if the tumor cell in rows, columns undergo apoptosis at hour == follow_up

CHANGE_CELL = zeros(rows, columns, follow_up);                                               

% Initialize the temporal dynamics output
% 1st row --> Number of tumoral cells
% 2nd row --> Number of potential mitotic cells
% 3rd row --> Number of potential apoptotic cells
% 4th row --> Number of quiescent cells
cells_matrix = zeros(4, follow_up); 
