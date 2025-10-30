%% Angio PC3 6 Hours %%

global ratio_pc3_vess_tum ratio_pc3_vess_tum_cabo;
load VesselsPropertiesPC3

% Find tumor coordinates 
[row_tum, col_tum] = find(bone == site.tumor);

% starting_vessels = vessels; 

% While loop to generate the vessel mask
while true
    
    % Every six hours i perform the standard angiogenesis
    if mod(hour, 6) == 0
        % Compute the ration btw the number of tumor and vessel sites
        curr_tumor_sites = sum(sum(bone == site.tumor));
        curr_vess_sites = sum(sum(bone == site.vessel));
        ratio = curr_vess_sites / curr_tumor_sites;
        % Exit condition: if the curr ratio is below the in vivo one exit
        if ratio > ratio_pc3_vess_tum
            break
        end
    end
    
    % Otherwise put the next vessel in the most empty area
    [most_empty_quadrant] = vessel_number_per_quadrant(tumor, vessels);
    
    % Erode this area to avoid placing vess center at the border
        
    % Define the morph operator
    se = strel('disk', 1);
    
    % Apply erosion
    most_empty_quadrant = imerode(most_empty_quadrant, se);
    
    current_tumor_area = sum(sum(bone == site.tumor)) + sum(sum(site.vessel));
    
    % If tumor grows I choose a more precise empty area to generate the vessel
    if current_tumor_area > 10 * tumor_area_start        
        [most_empty_quadrant] = vessel_number_per_quadrant(most_empty_quadrant, vessels);
    end

    % Find tumor coordinates of the most empty area
    [row_tum_area, col_tum_area] = find(most_empty_quadrant == 1);

    % Get the length of the vectors
    n = length(row_tum_area);

    % Generate a random permutation of indices
    perm_indices = randperm(n);

    % Shuffle both vectors and take the first value
    row_tum_area = row_tum_area(perm_indices);
    row_tum_vess = row_tum_area(1);
    col_tum_area = col_tum_area(perm_indices);
    col_tum_vess = col_tum_area(1);

    %Randomly select new vessel semi-axis from CS data
    a_vessel = max(random(pd_a_kernels), 10);
    b_vessel = min(random(pd_b_InverseGaussian), 60);

    % Define temp vessel variable to check possible overlaps
    temp_vess = zeros(rows, columns);

    % Find all vessel sites
    for tum_cell = 1 : size(row_tum, 1)

        % I compute the site-vessel center distance
        dist_centro_vaso = compute_distance(X, Y, ...
                                            row_tum(tum_cell), ...
                                            col_tum(tum_cell), ...
                                            row_tum_vess, ...
                                            col_tum_vess);

        % Find ellipse defining vess shape centered in the vessel center
        vessel_ellipse_area = ellipse_radius(ceil(a_vessel / site_dim), ...
                                             ceil(b_vessel / site_dim), ...
                                             abs(X(row_tum(tum_cell), col_tum(tum_cell)) - row_tum_vess), ...
                                             abs(Y(row_tum(tum_cell), col_tum(tum_cell)) - col_tum_vess));

        % Verify proximity condition
        if dist_centro_vaso <= vessel_ellipse_area
            temp_vess(row_tum(tum_cell), col_tum(tum_cell)) = 1;
        end  

    end   

    % Check if there was a vessel overlapping
    vess_check = vessels + temp_vess; % current vess + new vessel
    vess_overlap = any(vess_check(:) == 2);

    % If vess_overlap = 1 means there was an overlap
    if vess_overlap 
        continue
    end

    % If there was no overlap i add the vessel to the main matrix ...
    vessels = vessels + temp_vess;

    % ... to the cell that stores all unique vessel masks ... 
    single_vessel_masks{end+1} = temp_vess;
    
    % ... and to the bone matrix
    bone = bone + temp_vess;

    % New vessel properties
    new_vessel_prop = [row_tum_vess, col_tum_vess, a_vessel, b_vessel];
    % Update the main vessels matrix
    center_vessels  = [center_vessels; new_vessel_prop]; 
    % Clear temp variable
    clear new_vessel_prop  
    
end 