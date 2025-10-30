%% ABM Tumor Dynamics %%

%  This function is used to define the current mitosis and apoptosis
%  probabilities of the PCa cells whose internal clock reached the 24h.
%  Specifically, several conditions will define these probabilities:
%
%  1) Tumor Growth under control Regimen
%  2) Tumor Response to Rad223 Therapy
%  3) Tumor Response to Cabo   Therapy

if flag_cabo == 1 
    [row, col] = find(bone == site.tumor_edge & obs_influenced_cells == 0);
    if ~isempty(row)  % Check if there are any such cells
        for i = 1:size(row, 1)  % Loop through all identified cells
            bone(row(i), col(i)) = site.tumor;
        end
    end
    % Find Tumor Cells Influeced by the Obs
    % [row, col] = find(bone == site.tumor) & bone_marrow_edge);
    [row, col] = find(bone == site.tumor & obs_influenced_cells == 1);
    % Cycle through all those cells...
    for i = 1 : size(row, 1)
        % ... And Assign a different Value to the Resistant cells
        bone(row(i), col(i)) = site.tumor_edge;        
    end
    clear row col
end

% To spot mitotic and apoptotic cells
mitotic_cells = bone;
apoptotic_cells = bone;

% Initialize the matrix that will update mitotic/apoptotic events 
change_cell = zeros(rows, columns);

% Define Random access order to the ABM grid
[random_grid_order] = random_grid_access(rows, columns);

overall_mit_prob = ones(rows, columns) * alpha.p_mit_min;

% Now Investigate Each Grid Site Randomly
for agent = 1 : (rows * columns)

    % Define Coordinates (row and column) of the Current Site
    [row_agent, col_agent] = current_grid_coordinates(rows, columns, ...
        agent, random_grid_order);
    % Check if the Agent is a PCa cell
    if (bone(row_agent, col_agent) == site.tumor || bone(row_agent, col_agent) == site.tumor_edge)

        % Check if the Agent is in a Potential Mitotic/Apoptotic Stage
        if internal_clock(row_agent, col_agent) == T_cell_division

            % if the tumor is not a microtumor then normal dynamics,
            % otherwise set probabilities for apoptosis and mitosis
            if sum(sum(bone == site.tumor)) > 5
                % CASE 1: Tumor Growth Under Control Condition
                if (Rad.activity == 0 && flag_cabo == 0)
                    [p_mit, p_apo] = events_probability_control(rows, columns, ...
                        X, Y, site_dim, ...
                        bone, site, ...
                        row_agent, col_agent, ...
                        center_vessels, flag_min_distance, flag_cabozantinib, papo_cell_near_obs, obs_influenced_cells);
                end


                % CASE 2: Rad223 Therapy Activated
                if (Rad.activity == 1 && flag_cabo == 0)
                    [p_mit, p_apo] = events_probability_rad(rows, columns, X, Y, site_dim, bone, ...
                        site, row_agent, col_agent, B_mit, B_apo, ...
                        start_therapy_rad, hour, tau, rad_max, center_vessels, flag_min_distance, flag_cabozantinib, papo_cell_near_obs, obs_influenced_cells);
                end


                % CASE 3: Cabozantinib Therapy is Activated
                if (Rad.activity == 0 && flag_cabo == 1)
                    [p_mit, p_apo] = events_probability_cabo(rows, columns, X, Y, site_dim, ...
                        bone, alpha, site, row_agent, col_agent, ...
                        center_vessels, flag_min_distance, ...
                        vessel_retard, vessels_time, ...
                        pmit_cell_near_obs, papo_cell_near_obs);
                end

                % Combinatorial Therapy
                if (Rad.activity == 1 && flag_cabo == 1)
                    [p_mit_rad, p_apo_rad] = events_probability_rad(rows, columns, X, Y, site_dim, bone, ...
                        site, row_agent, col_agent, B_mit, B_apo, ...
                        start_therapy_rad, hour, tau, rad_max, center_vessels, flag_min_distance, papo_cell_near_obs, obs_influenced_cells);
                    [p_mit_cabo, p_apo_cabo] = events_probability_cabo(rows, columns, X, Y, site_dim, ...
                        bone, alpha, site, row_agent, col_agent, ...
                        center_vessels, flag_min_distance, ...
                        vessel_retard, vessels_time, ...
                        pmit_cell_near_obs, papo_cell_near_obs);
                    p_mit = min(p_mit_rad, p_mit_cabo);
                    p_apo = max(p_apo_rad, p_apo_cabo);
                end
            else

                p_mit = (alpha.p_mit_max + alpha.p_mit_min)/2;
                p_apo = (alpha.p_apo_max + alpha.p_apo_min)/4;

            end
      
            % MonteCarlo Algorithm
            [change_cell, mitotic_cells, apoptotic_cells] = monte_carlo_algorithm(p_mit, p_apo, row_agent, col_agent, change_cell, mitotic_cells, apoptotic_cells, site);

        end
    end
    
end

clear agent row_agent col_agent

UpdateFixedMatrix; % Defined in ABM_Initialization -> FixedMatrix





