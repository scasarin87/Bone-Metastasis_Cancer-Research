%% Define Apoptotic Site %%

%  This function is used to create the path between the cell which is
%  currently undergoing mitosis or apoptosis and the first free bone_marrow
%  site without crossing any bone/OBs/OCs

%  Input  -> bone           : current bone matrix
%            site           : struct defining the ABM sites
%            X, Y           : hexagonal grid
%            row/col_mitosis: coordinates of the mitotic site
%            directionx/y   : defined in ABM_Initialization
%            alpha          : struct containing model driving parameters
%            site_dim       : pix -> um (20.8333 um)

%  Output -> row/col_apoptosis: apoptotic site

function [row_apoptosis, col_apoptosis] = define_apoptotic_site(bone, site, X, Y, row_mitosis, col_mitosis, ...
                                                                directionx, directiony, flag_cabo, alpha, ...
                                                                center_vessels, vessels_time, vessel_retard, site_dim)
    apo_site_found = 0;
    
%     % When Cabozantinib Therapy is on
%     if flag_cabo == 1
%         % Find PCa Cells that won't resist to cabo (site.tumor)
%         [row_tum, col_tum] = find(bone == site.tumor | bone == site.tumor_edge);
%         % Random Access to these Cells
%         rand_tum = randperm(length(row_tum));
%         % Find Distance Between Cells and all the Vessels Centers
%         for cell = 1 : size(row_tum, 1)
%             apo_site_found = 1;
%             % Scan all the vessels ...
%             for vess = 1 : size(center_vessels, 1)
%                 % distance = compute_distance(X, Y, row_tum(rand_tum(cell)), col_tum(rand_tum(cell)), row_ves(vess), col_ves(vess));
%                 % ... and compute the Cell - Vessel distance
%                 distance = compute_distance(X, Y, row_tum(rand_tum(cell)), col_tum(rand_tum(cell)), center_vessels(vess, 1), center_vessels(vess, 2));
%                 % Scale distance to um
%                 distance = distance * site_dim;
%                 
%                 % Check if the current vessel was targeted by Cabo
%                 if center_vessels(vess, 5) == 1
%                     % If so I find the vessel index in the vessels_time matrix which contains the hours spent since targeting
%                     vess_index = find(ismember(vessels_time(:, 1:2), [center_vessels(vess, 1) center_vessels(vess, 2)], 'rows'));
%                     % Find the number of hours spent since cabo targeting
%                     time_from_target = vessels_time(vess_index, 3);
%                     % The longer the time since targeting the lower the vessel influence
%                     vess_influence_distance = alpha.vess_influence * (1 - time_from_target / vessel_retard);
%                     % If current cell is within vessel influcence ... 
%                     if distance < vess_influence_distance
%                         apo_site_found = 0;
%                         break
%                     end
%                 end
%                 % If current vessel was not targeted by cabo
%                 if center_vessels(vess, 5) == 0
%                    % If current cell is within stardard vessel influcence ... 
%                    if distance < alpha.vess_influence
%                         % ... No apoptosis
%                         apo_site_found = 0;
%                         break
%                     end
%                 end 
%             end
%         
%             % If a cell far from every vessel is found it will go apoptosis
%             if apo_site_found == 1 
%                 
%                 % These 2 lines select a tumor cell within the tumor mass,
%                 % chosen accordingly to the conditions above
%                 % row_apoptosis = row_tum(rand_tum(cell));
%                 % col_apoptosis = col_tum(rand_tum(cell));
%                 
%                 % These code instead select the closest tumor edge to the
%                 % select cell chosen with the conditions above
%                 % Define Whole Tumor and Complementary Tumor Mask Matrix
%                 % Define Tumor Mask Matrix
%                 [curr_tumor, ~] = def_tumor_masks(bone, site);
%                 % Define Current Tumor Internal Boundary 
%                 tumor_edge = bwperim(curr_tumor);
%                 % Get Tumor Edges Coordinates
%                 [row_tum_edge, col_tum_edge] = find(tumor_edge == 1);
%                 for cell_edge = 1 : size(row_tum_edge, 1)        
%                     distance(cell_edge) = compute_distance(X, Y, row_tum(rand_tum(cell)), col_tum(rand_tum(cell)), row_tum_edge(cell_edge), col_tum_edge(cell_edge));
%                 end
% 
%                 [~, n_target] = min(distance); 
%                 clear distance
% 
%                 row_apoptosis = row_tum_edge(n_target);
%                 col_apoptosis = col_tum_edge(n_target);                
%                 
%                 break
%             end           
%         end            
%     end
%     
%     % If cabo is on but all tumor cells resist to the therapy
%     if flag_cabo == 1 && apo_site_found == 0
%         row_apoptosis = 0;
%         col_apoptosis = 0;
%     end
%     
    % When Cabozantinib Therapy is Off (Control - Rad Condition)
    if flag_cabo == 0  || flag_cabo == 1
        % Explore the neighbors of the target site (row_target, col_target)
        if mod(col_mitosis, 2) == 0
           j_k2 = 2;
        else
           j_k2 = 1;
        end     

        dist_temp = zeros(1, 6);
        for neighbour_index = 1 : length(dist_temp) 
            % Get the Neighbour Coordinates
            row_neighbour = row_mitosis + directionx(neighbour_index, j_k2);
            col_neighbour = col_mitosis + directiony(neighbour_index);
            % If the current neighbour (among the 6) is a PCa cell ...
            if bone(row_neighbour, col_neighbour) == site.tumor || bone(row_neighbour, col_neighbour) == site.tumor_edge            
                % ... Get the neighbour - mitosis site distance ...
                dist_temp(neighbour_index) = compute_distance(X, Y, row_mitosis, col_mitosis, row_neighbour, col_neighbour); 
            else
                % ... Otherwise we set it to high (won't be considered)
                dist_temp(neighbour_index) = 500;
            end                              
        end     
        clear neighbour_index

        % Get the index of the closest tumor neighbour to the mitotic site ...
        [~, neigbour_target] = min(dist_temp); 
        clear dist_temp

        % ... The closest site to the mitotic site will be the apoptotic site
        row_apoptosis = row_mitosis + directionx(neigbour_target, j_k2);
        col_apoptosis = col_mitosis + directiony(neigbour_target); 
    end

end
         