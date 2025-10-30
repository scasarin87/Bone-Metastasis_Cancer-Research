%% Events Probability Cabo %%

%  This function is used to compute the mitosis and apoptosis probabilities
%  for each PCa cell at the current 'hour' when cabo therapy is active. 

%  Input  -> rows, columns : Y and X size of the grid
%            X, Y          : hexagonal grid
%            site_dim      : ABM scaling factor 
%            bone          : main ABM matrix, it labels each model agent
%            alpha         : struct with the 6 model driving parameters 
%            site          : struct with the ABM agents names
%            row/col_agent : current site coordinates (row/column)
%            center_vessels: n.vessels x 4 matrix (contains vessels
%                            coordinates and semi-axis values)
%            flag_min_dist : if any vessels within the alpha.vess_influence
%                            area is set to 1 -declared in FixedParametes.m
%            vessel_retard : check FixedParameters.m in ABM_Initialization
%            vessels_time  : check FixedParameters.m in ABM_Initialization
%
%  Output -> p_mit, p_apo  : mitosis and apoptosis probability of the
%                            current tumor 'agent'. These probs will be the
%                            input of the M.Carlo algorithm to define
%                            whether the cell will divide or die
%            vessels_time  : keeps track of the vessels targeted by cabo



function [p_mit, p_apo, vessels_time] = events_probability_cabo(rows, columns, X, Y, site_dim, ...
                                                                bone, alpha, site, row_agent, col_agent, ...
                                                                center_vessels, flag_min_distance, ...
                                                                vessel_retard, vessels_time, ...
                                                                pmit_cell_near_obs, papo_cell_near_obs)
    
    % Until line 89 is equal to the first part of CASE 1(control)
    clear cont_index_vessels distance_matrix
    
    % Variables Initialization
    p_apo               = 1;   % Put it to high
    dist_matrix_init    = 100; % Initialize it to high
    distance_matrix     = ones(size(center_vessels, 1), 25) * dist_matrix_init;     
    cont_index_vessels  = ones(size(center_vessels,1), 1); % vector with 
    
    % Define Again Random access to the ABM grid
    [random_grid_order_1] = random_grid_access(rows, columns);
    
    for agent_1 = 1 : (rows * columns)
        % Define Coordinates (row and column) of the Current Site 
        [row_agent_1, col_agent_1] = current_grid_coordinates(rows, columns, ...
                                                               agent_1, random_grid_order_1);
        % Check if the current site is a vessel
        if bone(row_agent_1, col_agent_1) == site.vessel || bone(row_agent_1, col_agent_1) == site.vessel_cabo
            
            % Compute the distance btw tumor site (X(row_agent, col_agent), Y(row_agent, col_agent))
            % and a random point of a vessel X(row_agent_1, col_agent_1), Y(row_agent_1, col_agent_1))
            [agent_vessel_dist] = compute_distance(X, Y, row_agent, col_agent, ...
                                                   row_agent_1, col_agent_1);       
                                                 
            % Check if this distance is < site.vess_influence...   
            if agent_vessel_dist * site_dim < alpha.vess_influence  
                
               % ... and compute the distance between the site and the vessel and store it in distance_matrix
               for k = 1:size(center_vessels, 1) 
                   
                   % First we define to which vessel the site (row_agent_1, col_agent_1) belongs 
                   dist_vess_cenvess = sqrt((Y(row_agent_1, col_agent_1) - center_vessels(k, 1)) ^ 2 + ...
                                            (X(row_agent_1, col_agent_1) - center_vessels(k, 2)) ^ 2);
                   Ellipse_Radius = ellipse_radius(ceil(center_vessels(k, 3) / site_dim), ...
                                                   ceil(center_vessels(k ,4) / site_dim), ...
                                                   abs(X(row_agent_1, col_agent_1) - center_vessels(k, 2)), ... 
                                                   abs(Y(row_agent_1, col_agent_1) - center_vessels(k, 1)));
                                                                           
                   if dist_vess_cenvess <= Ellipse_Radius 
                       
                      % Then we add the computed distance at the k-th row of distance_matrix, 
                      % where k is the index of the vessel
                      [distance_matrix(k, cont_index_vessels(k))] = compute_distance(X, Y, row_agent, ...
                                                                                     col_agent, row_agent_1, col_agent_1);
                      cont_index_vessels(k) = cont_index_vessels(k) + 1;                                             
                   end       
               end   
               
               clear k
               
            end              
        end 
    end 
    
    clear agent_1
    
    % We need the minimum distance between each vessel within alpha.vess_influence um and the current tumor site
    min_distance_matrix = min(distance_matrix, [], 2);
    p_div = 0;
    
    for i = 1 : size(min_distance_matrix, 1)
        
        % I consider only the vessels that actually give contribution
        if min_distance_matrix(i) ~= dist_matrix_init 
            
           % If the vessel was not target by cabo, it provides full
           % contribution to the mitosis probability
           if center_vessels(i, 5) == 0 
              p_div = p_div + (alpha.p_mit_max - alpha.p_mit_min) * logsig(-alpha.slope_prob * min_distance_matrix(i) * site_dim + 4.5) + alpha.p_mit_min; 
              flag_min_distance = 1;
              
           % Otherwise, we linearly re-scale its contribution over time
           elseif center_vessels(i, 5) == 1 
                  
                  %I check how much time has passed since it was targeted...
                  for i_time = 2 : size(vessels_time, 1)  
                      if vessels_time(i_time, 1) == center_vessels(i, 1) && vessels_time(i_time, 2) == center_vessels(i, 2) && vessels_time(i_time, 3) <= vessel_retard 
                         time_target_cabo = vessels_time(i_time, 3);
                      end 
                  end 
                 
                  
                  %... And I compute its mitosis contribution considering a
                  % linear decreasing effect of the vessel wrt to the time
                  p_div = p_div + (1 - time_target_cabo / vessel_retard) * ((alpha.p_mit_max - alpha.p_mit_min) * logsig(-alpha.slope_prob * min_distance_matrix(i) * site_dim + 4.5) + alpha.p_mit_min); 
                  flag_min_distance = 1;                                               
           end                                           
        end   
    end   
    
    clear i
                                                        
   p_mit = min(alpha.p_mit_max, p_div); %Upper bound
   p_mit = max(p_mit, alpha.p_mit_min);  %Lower bound  
   p_mit = max(p_mit, pmit_cell_near_obs(row_agent, col_agent)); % Comparison with Obs induced resistance probabilities
    
   % If there are vessels within alpha.vess_influence [um] area I'll compute
   % p_apo with respect to the min distance of the closest vessel 
   if flag_min_distance == 1 
      p_apo = min(p_apo, (alpha.p_apo_max - alpha.p_apo_min) * logsig(alpha.slope_prob * min(min_distance_matrix) * site_dim - 4.5) + alpha.p_apo_min);
      p_apo = min(p_apo, papo_cell_near_obs(row_agent, col_agent));
      flag_min_distance = 0;
   
   % Otherwise I'll just consider the center of the closest vessels
   elseif flag_min_distance == 0 
          for p = 1 : size(center_vessels,1)
              x_dist(p) = sqrt((X(row_agent, col_agent) - center_vessels(p, 2)) ^ 2 + ...
                               (Y(row_agent, col_agent) - center_vessels(p, 1)) ^ 2);                          
              p_apo = min(p_apo, (alpha.p_apo_max - alpha.p_apo_min) * logsig(alpha.slope_prob * x_dist(p) * site_dim - 4.5) + alpha.p_apo_min);
              p_apo = min(p_apo, papo_cell_near_obs(row_agent, col_agent));
          end  
          
          clear p
   end      
   
   Rad.mitosis = p_mit;
   Rad.apoptosis = p_apo;
   Rad.quiescence = 1 - Rad.mitosis - Rad.apoptosis;
    
   
end