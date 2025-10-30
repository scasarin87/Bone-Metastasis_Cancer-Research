%% Update Internal Clock %%

%  This function is used to update PCa cells internal clock at the start of
%  every iteration (a new hour in tumor development). When clock reaches
%  the 24h the cells will undergo a mitotic or apoptotic event depending on
%  several condition related to the mitosis and apoptosis probabilities.

%  Input  -> internal_clock: matrix that stores each PCa cell "age" [h]
%            rows, columns : Y and X size of the grid
%            bone          : main ABM matrix, each agent has its code
%                            identifier
%            site          : struct containing the site names
%            T_cell_divis  : After * hours a PCa will undergo mit/apo event
%
%  Output -> updated_internal_clock : updated_internal_clock = internal_clock + 1h

function [updated_internal_clock, updated_internal_clock_ocs] = update_internal_clock(internal_clock, internal_clock_ocs, rows, columns, ...
                                                    bone, site, T_cell_division, T_ocs_resorption)

    updated_internal_clock = internal_clock;
    updated_internal_clock_ocs = internal_clock_ocs;
    
    for jj = 1:rows
        for kk = 1 : columns
            
            % Check if ABM site is tumor
            if bone(jj, kk) == site.tumor || bone(jj, kk) == site.tumor_edge
                if internal_clock(jj, kk) < T_cell_division % if cell is still within its cycle,...
                    updated_internal_clock(jj,kk) = internal_clock(jj, kk) + 1; % ... then just update the clock
                else
                    updated_internal_clock(jj, kk) = 0; % otherwise, reset the clock
                end
            end
            
            % Check if ABM site is ocs
            if bone(jj, kk) == site.osteoclast
                if internal_clock_ocs(jj, kk) < T_ocs_resorption % if cell is still within its cycle,...
                    updated_internal_clock_ocs(jj,kk) = internal_clock_ocs(jj, kk) + 1; % ... then just update the clock
                else
                    updated_internal_clock_ocs(jj, kk) = 0; % otherwise, reset the clock
                end
            end
        end  
    end 
    
    clear jj kk
    
end
