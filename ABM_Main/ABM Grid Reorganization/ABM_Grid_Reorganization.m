%% ABM Grid Re-Organization %%

%  This function performs all the cellular events determined by the Monte
%  Carlo method. The response of the tumor mass will change according to 
%  the position of the mitotic/apoptotic cell. 

% Compute the total number of mitotic/apoptotic events
number_of_events = sum(sum(abs(change_cell)));

% Keep track of the coordinates of each site undergoing mitosis/apoptosis
[row_event, col_event] = find(change_cell ~= 0);

% Define random access order to the potential mitotic/apoptotic sites
[random_grid_order] = random_grid_access(number_of_events, 1);

% For loop across every cellular event at the current simulation hour
for event = 1 : number_of_events
    
    % Control: check if N.events > free bone marrow sites. If so... break
    if event > sum(sum(bone == site.bone_marrow))
        break
    end
    
    % Define coordinates (row and column) of the current active site 
    row_agent = row_event(random_grid_order(event));
    col_agent = col_event(random_grid_order(event));
    
    % Double control if agent is tumor site and perform cellular dynamics
    if (bone(row_agent, col_agent) == site.tumor || bone(row_agent, col_agent) == site.tumor_edge)
        
        % Define Mitotic Site
        [row_mitosis, col_mitosis] = define_mitotic_site(bone, site, X, Y, row_agent, col_agent, Rad, hour, start_therapy_rad, follow_up, flag_cabo);
        
        if size(row_mitosis, 1) > 0
            
            % Define Apoptotic Site    
            [row_apoptosis, col_apoptosis] = define_apoptotic_site(bone, site, X, Y, row_mitosis, col_mitosis, ...
                                                               directionx, directiony, flag_cabo, alpha,...
                                                               center_vessels, vessels_time, vessel_retard, site_dim);
        
        
            % Mitosis Event (change_cell(x, y) = 1 -> (x,y) tumor mitosis)
            if change_cell(row_agent, col_agent) == 1 

                % Perform Mitotic Event
                [bone, change_cell, internal_clock] = mitosis_event(bone, change_cell, ...
                    internal_clock, site, row_agent, col_agent, row_mitosis, col_mitosis);

            end 

            % Apoptosis Event (change_cell(x, y) = - 1 -> (x,y) tumor apoptosis)
            if change_cell(row_agent, col_agent) == -1 && row_apoptosis ~= 0

                % Perform Apoptotic Event
                [bone, change_cell, internal_clock] = apoptosis_event(bone, change_cell, ...
                    internal_clock, site, row_agent, col_agent, row_apoptosis, col_apoptosis);

            end 
        end
    end 

    if flag_cabozantinib == 1
        [obs_influenced_cells, papo_cell_near_obs] = change_influence_matrix(bone, site, site_dim, X, Y, rows, columns);
    end
    
end