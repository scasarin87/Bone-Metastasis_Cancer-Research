%% Apoptosis Event %%

%  This function is used to remove a cell (apoptosis) at the end of the 
%  path previously defined

%  Input  -> old_bone           : current bone matrix
%            old_change_cell    : currently contains sites which undergo
%                                 mitosis or apoptosis
%            old_internal_clock : current matrix with PCa cell life time
%            site               : struct defining the ABM sites
%            row/col_agent      : current ABM row and col 

%  Output -> new_bone           : updated bone matrix after mitosis event
%            new_change_cell    : updated matrix
%            new_internal_clock : updated internal clock (once mitosis
%                                 happens, it is set to 0)


function [new_bone, new_change_cell, new_internal_clock] = apoptosis_event(old_bone, old_change_cell, ...
            old_internal_clock, site, row_agent, col_agent, row_apoptosis, col_apoptosis)
        
    % First, the site (row_agent, col_agent) vacates the cell and it is 
    % immediately replaced by another one
    
    old_internal_clock(row_agent, col_agent) = 0; % reset the clock ...
    old_change_cell(row_agent, col_agent)    = 0; % ... and the dynamics matrix
    
    old_bone(row_apoptosis, col_apoptosis)   = site.bone_marrow;
    old_internal_clock(row_apoptosis, col_apoptosis) = 0;
    old_change_cell(row_apoptosis, col_apoptosis)    = 0;
        
    new_bone           = old_bone(:, :);
    new_change_cell    = old_change_cell(:, :);
    new_internal_clock = old_internal_clock(:, :);       
        
end        