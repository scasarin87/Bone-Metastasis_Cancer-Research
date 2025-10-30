%% Redistribute Vessels %%

%  This function is used to re distribute the vessels in the PC3 tumor, 
%  so that when it grows they are still localized on its edges

% Define global variables
global CxT CyT

% Find tumor perimeter
tumor = (bone == site.tumor) | (bone == site.vessel);
tumor = imfill(tumor, 'holes'); % If any holes remained
tumor_per = bwperim(tumor);

% Extract tumor and tumor perimeter coordinates
[row_tum, col_tum] = find(tumor);
[row_tum_per, col_tum_per] = find(tumor_per); 

% Extract vessel center coordinates
row_vess_cen = center_vessels(:, 1);
col_vess_cen = center_vessels(:, 2);

% Init distance cells
curr_vess_tumperim_distances = {};
indexes = {};

% For loop to store the distances btw any tumor perim site and vess centers
for vess = 1 : size(row_vess_cen, 1)   
    
    % Define row and col of the curr vessel
    curr_vess_row_cen = row_vess_cen(vess);
    curr_vess_col_cen = col_vess_cen(vess);
    vess_tumperim_distances = [];
    
    % For loop across the sites at the tumor perimeter
    for tum = 1 : size(row_tum_per, 1)
        
        % Define row and col of the curr tumor perim coordinates
        curr_row_tum_per = row_tum_per(tum);
        curr_col_tum_per = col_tum_per(tum);
        
        % Compute closest distance between current vessel and tumor perimeter
        vess_tumperim_distance = compute_distance(X, Y, ...
                                                  curr_row_tum_per, ...
                                                  curr_col_tum_per, ...
                                                  curr_vess_row_cen, ...
                                                  curr_vess_col_cen);
        % Add current distance to the array
        vess_tumperim_distances = [vess_tumperim_distances, vess_tumperim_distance];
    end
    
    % Sort the distance vectors
    [curr_vess_tumperim_distances{vess}, indexes{vess}] = sort(vess_tumperim_distances);
end   


for vess = 1 : size(curr_vess_tumperim_distances, 2)
    
    % If I see overlaps -> increase this number to generate different paths
    % with the second, third, fourth, ... closest perim sites
    ind_dist = 1;
    
    % While loop that checks for possible overlaps before moving the vess
    while true
        
        % With too many vessels it is impossible to avoid the overlap
        try
            % Find closest perimeter distance and to current vessel and its index
            closest_dist = curr_vess_tumperim_distances{vess}(ind_dist); 
            closest_ind = indexes{vess}(ind_dist);
        
        % I'd just skip the redistribution if this happens
        catch exception
            break
        end

        % Arbitrary threshold not to keep all vessels at the same distance from border
        if closest_dist < 5
            break
        end
        
        % I give just the 33% probability to a vessel to be moved (to avoid
        % every vessel moving at the same time)
        if rand() > 1
            break
        end

        % Extract the current vessel binary mask
        curr_vess = single_vessel_masks{vess};
        
        % Remove the current vessel from the 'vessels' mask 
        vessels = vessels - curr_vess;

        % Generate a path between vessel center and closest perimeter point
        diff_row = abs(row_vess_cen(vess) - row_tum_per(indexes{vess}(ind_dist)));
        diff_col = abs(col_vess_cen(vess) - col_tum_per(indexes{vess}(ind_dist)));
        max_diff = max(diff_row, diff_col);
        row_coords = round(linspace(row_vess_cen(vess), ...
                                    row_tum_per(indexes{vess}(ind_dist)), ...
                                    max_diff + 1));
        col_coords = round(linspace(col_vess_cen(vess), ...
                                    col_tum_per(indexes{vess}(ind_dist)), ...
                                    max_diff + 1));
                                
        % For loop across the new possible center coordinates
        for coord = 1 : size(row_coords, 2)

            % When I find the first coords closer to the perim than the actual vessel center
            if row_coords(coord) == row_vess_cen(vess) && col_coords(coord) == col_vess_cen(vess)
                continue
            else
                % I define new vess center coordinates
                new_vess_row = row_coords(coord);
                new_vess_col = col_coords(coord);    
                
                % Coords found, exit the loop
                break
            end
            
        end
            
        % Define temp vessel variable to check possible overlaps
        temp_vess = zeros(rows, columns);

        % I define the new vessel with the new center
        for tum_cell = 1 : size(row_tum, 1)

            % I compute the site-vessel center distance
            dist_centro_vaso = compute_distance(X, Y, ...
                                                row_tum(tum_cell), col_tum(tum_cell), ...
                                                new_vess_row, new_vess_col);
            vessel_ellipse_area = ellipse_radius(ceil(center_vessels(vess, 3) / site_dim), ...
                                                 ceil(center_vessels(vess, 4) / site_dim), ...
                                                 abs(X(row_tum(tum_cell), col_tum(tum_cell)) - new_vess_row), ...
                                                 abs(Y(row_tum(tum_cell), col_tum(tum_cell)) - new_vess_col));

            % Define if the current tumor site belongs to the vessel
            if dist_centro_vaso <= vessel_ellipse_area
                temp_vess(row_tum(tum_cell), col_tum(tum_cell)) = 1;
            end

        end 

        % Check if there was a vessel overlapping
        temp_vess_dilate = imdilate(temp_vess, strel("disk", 1));
        vess_check = vessels + temp_vess_dilate; % current vess + new vessel
        vess_overlap = any(vess_check(:) == 2);

        % If vess_overlap = 1 means there was an overlap
        if vess_overlap == 0
            
            % If no overlap ... correct 'center_vessels' matrix
            center_vessels(vess, 1) = new_vess_row;
            center_vessels(vess, 2) = new_vess_col;
            
            % Update 'vessels' matrix ...
            vessels = vessels + temp_vess;
            
            % ... the cell that stores all unique vessel masks ...
            single_vessel_masks{vess} = temp_vess;
            
            % ... And the main bone matrix
            bone = bone - curr_vess;            
            bone = bone + temp_vess;
            
            %Change also the vessel time variable if cabo therapy is on
            if flag_cabo == 1
                for i_time = 1 : size(vessels_time, 1)  
                    if vessels_time(i_time, 1) == center_vessels(vess, 1) && vessels_time(i_time, 2) == center_vessels(vess, 2) 
                        vessels_time(i_time, 1) = new_vess_row;
                        vessels_time(i_time, 2) = new_vess_col;
                    end  
                end   
            end
            
            % Finally, exit the while loop 
            break
            
        else
            
            % Restore the 'vessels' matrix 
            vessels = vessels + curr_vess;
            
            % Increment the ind_dist variable and compute a new path if I had overlaps
             ind_dist = ind_dist + 1;
        end        
            
    end 
    
end  