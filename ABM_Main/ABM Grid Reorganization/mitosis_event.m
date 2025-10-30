%% Mitosis Event %%

%  This function is used to add a cell (mitosis) at the end of the path
%  previously defined

%  Input  -> old_bone           : current bone matrix
%            old_change_cell    : currently contains sites which undergo
%                                 mitosis or apoptosis
%            old_internal_clock : current matrix with PCa cell life time
%            site               : struct defining the ABM sites
%            row/col_agent      : current ABM row and col 
%            row/col_mitosis    : mitotic site

%  Output -> new_bone           : updated bone matrix after mitosis event
%            new_change_cell    : updated matrix
%            new_internal_clock : updated internal clock (once mitosis
%                                 happens, it is set to 0)

function [new_bone, new_change_cell, new_internal_clock] = mitosis_event(old_bone, old_change_cell, ...
            old_internal_clock, site, row_agent, col_agent, row_mitosis, col_mitosis)
        
        %Mitosis -> We add a cell at the found bone marrow site
        old_bone(row_mitosis, col_mitosis)           = site.tumor;
        old_change_cell(row_mitosis, col_mitosis)    = 0;
        old_internal_clock(row_mitosis, col_mitosis) = 0;
        
        % Set to 0 the clock of the dividing cell
        old_internal_clock(row_agent, col_agent) = 0;
        old_change_cell(row_agent, col_agent)    = 0;  
                       
        new_bone           = old_bone(:, :);
        new_change_cell    = old_change_cell(:, :);
        new_internal_clock = old_internal_clock(:, :);        

end        