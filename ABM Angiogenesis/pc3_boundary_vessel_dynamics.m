%% PC3 Boundary Vessels Dynamics %%

% This function is used to:
% 1) Guarantee that pc3 boundary vessels keep belonging to the tumor edge 
% throughout the simulation. 
% 2) Place new boundary vessels (angiogenesis) until reaching the in vivo
% coverage.

% Every hour I find the new edge of the growing tumor, if a vessel center no
% more belongs to it, I choose the new vessel center as the closest edge 
% site to the previous center. Finally I move to whole vessel to the boundary
% and update both center_vessels and single_vessel_masks variable.

% Then I pick a random tumor edge site and I make it the center of a new
% boundary vessel, updating all the related variables.


% Load PC3 Vessel Data
% load 'ABM Main'/VesselData/VesselsProperties_PC3.mat

%% 1) Shift boundary vessels to make them stick to the tumor boundary
%% Find new edge on growing tumor

if flag_cabozantinib 
    if enter_condition
        pc3_vess_boundary_coverage = pc3_vess_boundary_coverage / 1.9; % 1.9 taken from x Luca Cabo.xlsx sheet_name = 'PC3, blood vess area'
        pc3_central_vess_coverage = pc3_central_vess_coverage / 1.9;
        enter_condition = 0; % Defined in ABM Initialization
    end
end

% Get tumor area with vessels
tumor = (bone == site.tumor | bone == site.tumor_edge | bone == site.vessel);

% Get the largest tumor portion (some tumor cells might detach from its core)
labeled_mask = bwlabel(tumor);
props = regionprops(labeled_mask, 'Area');
areas = [props.Area];
[~, max_id] = max(areas);
tumor = labeled_mask == max_id;

% Fill tumor area to remove holes, if any
tumor = imfill(tumor, 'holes');

% Erode tumor edge to place vessels not on the boundary 
se = strel('disk', 1);
tumor_eroded = imerode(tumor, se);
if sum(sum(tumor_eroded == 1)) < 6
    tumor_eroded = tumor;
end
tumor_bound = bwperim(tumor_eroded); 

% Get the largest tumor boundary area (some tumor cells might detach from its core)
labeled_mask = bwlabel(tumor_bound);
props = regionprops(labeled_mask, 'Area');
areas = [props.Area];
[max_area, max_id] = max(areas);
tumor_bound = labeled_mask == max_id;

%% Check if vessel centers are on the tumor edge

% Find current vessel number
vess_number = size(center_vessels, 1);

% Extract boundary coordinates
[row_tum_bound, col_tum_bound] = find(tumor_bound == 1);

% Put them in the same array
coord_tum_bound = [row_tum_bound col_tum_bound];

% Extract the coordinates of center vessels
coord_vess_cent = center_vessels(:, 1 : 2);

% Flag for vess belonging to the tumor boundary
flag_vess_on_boun = 0;

% Save index for vessels that do not belong to the boundary
vess_index = [];

% Find if vess centers belong to tumor boundary
for i_vess = 1 : size(coord_vess_cent, 1)
    
    % Check if current vessel is a boundary vessel
    if strcmp(single_vessel_masks{i_vess}{5}, 'boundary_vess') == 0
        continue
    end
    
    % Check among tumor boundary coords
    for i_boun = 1 : size(coord_tum_bound, 1)
        
        % Check if vess center belongs to boundary
        if isequal(coord_vess_cent(i_vess, :), coord_tum_bound(i_boun, :))
            
            % Raise the flag
            flag_vess_on_boun = 1;            
            % Break the inner for loop
            break
            
        end
    end
    
    % If vess center belongs to tumor edge
    if flag_vess_on_boun
        
        % Re initialize the flag
        flag_vess_on_boun = 0;
    
    % If not
    else
        
        % Save the curren vessel center index
        vess_index = [vess_index; i_vess];        
        
    end
end

%% Find the closest tumor edge site to vessels center not on the boundary   

% Initialize to zero the variable storing the vessels coordinates to be moved
coord_vess_cent_to_move = [];
    
% If I found a/some vessel/s center not on tumor boundary
if size(vess_index, 1) > 0
    
    % I extract their center coordinates
    coord_vess_cent_to_move = coord_vess_cent(vess_index, :);
    
    % Variable that stores all vess center - tumor boundary site dist
    distance = [];
    
    % Array with coord indexes of the closest tum boundary site 
    all_index = [];
    
    % For loop to find the closest tumor boundary to the vess center
    for i_vess = 1 : size(coord_vess_cent_to_move, 1)        
        % For loop over each tumor boundary site
        for i_boun = 1 : size(coord_tum_bound, 1)
            
            % Compute distance
            curr_dist = compute_distance(X, Y, ...
                                         coord_vess_cent_to_move(i_vess, 1), ...
                                         coord_vess_cent_to_move(i_vess, 2), ...
                                         coord_tum_bound(i_boun, 1), ...
                                         coord_tum_bound(i_boun, 2));
                                     
            % Append to the distance variable
            distance = [distance; curr_dist];
            
        end
        
        % Get the index of the closest tumor boundary site
        [~, index] = min(distance);
        
        % Update the all_index variable
        all_index = [all_index; index];
        
        % Re initialize the distance variable
        distance = [];
        
    end
    
end
    
%% Move the vessel center to the new tumor boundary site

% For loop across all the vessels that should be moved
for i_vess = 1 : size(coord_vess_cent_to_move, 1)
    
    % Get current vess center coordinates
    curr_row_vess_cent = coord_vess_cent_to_move(i_vess, 1);
    curr_col_vess_cent = coord_vess_cent_to_move(i_vess, 2);
    
    % Get the new vess_center coordinates
    new_row_vess_cent = coord_tum_bound(all_index(i_vess), 1);
    new_col_vess_cent = coord_tum_bound(all_index(i_vess), 2);
    % Extract axis info
    major_axis = single_vessel_masks{vess_index(i_vess)}{2};
    
    % Change center vessel variable
    center_vessels(vess_index(i_vess), 1) = new_row_vess_cent;
    center_vessels(vess_index(i_vess), 2) = new_col_vess_cent;
    
    % Initialize temp_vess variable to store current vess binary matrix
    temp_vess = zeros(rows, columns);
    % Initialize temp_vess_final variable to store final vess binary matrix
    temp_vess_final = zeros(rows, columns);

    % Find all vessel sites
    for tum_cell = 1 : size(row_tum_bound, 1)

        % I compute the site-vessel center distance
        site_vess_dist = compute_distance(X, Y, ...
                                          row_tum_bound(tum_cell), ...
                                          col_tum_bound(tum_cell), ... 
                                          new_row_vess_cent, ...
                                          new_col_vess_cent);

        % Verify proximity condition
        if site_vess_dist * site_dim <= major_axis / 2
            temp_vess(row_tum_bound(tum_cell), col_tum_bound(tum_cell)) = 1; 
        end 
        
        % Verify proximity condition
        if site_vess_dist * site_dim <= single_vessel_masks{vess_index(i_vess)}{4} / 2
            temp_vess_final(row_tum_bound(tum_cell), col_tum_bound(tum_cell)) = 1; 
        end 

    end
    
    % Extract the original vessel coords
    [row_orig_vess, col_orig_vess] = find(single_vessel_masks{vess_index(i_vess)}{1} == 1);
    
    % Extract the shifted vessel coords
    [row_shif_vess, col_shif_vess] = find(temp_vess == 1);    
    
    % Remove former vessel from bone matrix
    for i = 1 : length(row_orig_vess)
        bone(row_orig_vess(i), col_orig_vess(i)) = site.tumor;
    end
    
    % Remove former vessel from bone matrix
    for i = 1 : length(row_shif_vess)
        bone(row_shif_vess(i), col_shif_vess(i)) = site.vessel;
    end
        
    % Add the current vessel size to the cell
    single_vessel_masks{vess_index(i_vess)}{1} = temp_vess;
    % Add the final vessel size to the cell
    single_vessel_masks{vess_index(i_vess)}{3} = temp_vess_final;
        
end

%% Make boundary vessel larger, until it reaches its final size

% For loop across all vessels
for i_vess = 1 : size(single_vessel_masks, 2)
    
    % It must be a boundary vessel
    if strcmp(single_vessel_masks{i_vess}{5}, 'boundary_vess') == 0
        continue
    end
    
    % Each vessel has 33% of chance of growing at the current iteration
    if rand(1) < 0.33
        continue
    end
    
    % Get vessel center row and column
    curr_row_vess_cent = center_vessels(i_vess, 1);
    curr_col_vess_cent = center_vessels(i_vess, 2);
    
    % Get current major axis
    major_axis = single_vessel_masks{i_vess}{2};
    
    % Update major axis
    major_axis = min((major_axis + site_dim), single_vessel_masks{i_vess}{4});
    
    % Adjust its value
    if major_axis > single_vessel_masks{i_vess}{4} - site_dim
        major_axis = single_vessel_masks{i_vess}{4};
    end        
            
    % Initialize temp_vess matrix and for vess final dimension
    temp_vess = zeros(rows, columns);
    temp_vess_final = zeros(rows, columns);    

    % Find all vessel sites
    for tum_cell = 1 : size(row_tum_bound, 1)

        % I compute the site-vessel center distance
        site_vess_dist = compute_distance(X, Y, ...
                                          row_tum_bound(tum_cell), ...
                                          col_tum_bound(tum_cell), ... 
                                          curr_row_vess_cent, ...
                                          curr_col_vess_cent);

        % Verify proximity condition
        if site_vess_dist * site_dim <= major_axis / 2
            temp_vess(row_tum_bound(tum_cell), col_tum_bound(tum_cell)) = 1; 
        end  
        
        % Verify proximity condition
        if site_vess_dist * site_dim <= single_vessel_masks{i_vess}{4} / 2
            temp_vess_final(row_tum_bound(tum_cell), col_tum_bound(tum_cell)) = 1; 
        end  

    end
    
    % Define flag for vessel overlap
    flag_ovrl = 0;
    
    % Check for overlapping between the enlarged vessel and previous ones
    for j_vess = 1 : size(single_vessel_masks, 2)
        
        % It must be a boundary vessel
        if strcmp(single_vessel_masks{j_vess}{5}, 'boundary_vess') == 0
            continue
        end
        
        % I must not compare it to the same vessel
        if j_vess == i_vess
            continue
        end
        
        % Get matrix of the current vessel
        curr_vess_matrix = single_vessel_masks{j_vess}{1};
        
        % Check for overlap
        vess_check = curr_vess_matrix + temp_vess;             
        vess_overlap = any(vess_check(:) == 2);
        
        % If overlap, I exit the loop and put off the vessel size change
        if vess_overlap
            
            % Raise ovrl flag
            flag_ovrl = 1;
            
            % Exit the loop
            break
        
        end
        
    end
    
    % If I found an overlap, I do not change vessel size rn
    if flag_ovrl
        continue
    end
    
    % Update major axis variables
    single_vessel_masks{i_vess}{2} = major_axis;
    center_vessels(i_vess, 3) = major_axis;
    
    % Extract the original vessel coords
    [row_orig_vess, col_orig_vess] = find(single_vessel_masks{i_vess}{1} == 1);
    
    % Extract the shifted vessel coords
    [row_shif_vess, col_shif_vess] = find(temp_vess == 1);        
    
    % Remove former vessel from bone matrix
    for i = 1 : length(row_orig_vess)
        bone(row_orig_vess(i), col_orig_vess(i)) = site.tumor;
    end
    
    % Remove former vessel from bone matrix
    for i = 1 : length(row_shif_vess)
        bone(row_shif_vess(i), col_shif_vess(i)) = site.vessel;
    end
        
    % Add the vessel to the cell
    single_vessel_masks{i_vess}{1} = temp_vess;  
    single_vessel_masks{i_vess}{3} = temp_vess_final;  
    
            
end

%% 2) Boundary Vessels Angiogenesis
%% Get the coordinates of boundary vessels

% Variable init
row_bound_vess = [];
col_bound_vess = [];

row_bound_vess_fin = [];
col_bound_vess_fin = [];

% I fill the two array with the row and col coords of boundary vessels
for i_vess = 1 : size(single_vessel_masks, 2)
    
    % I consider only boundary vessels
    if strcmp(single_vessel_masks{i_vess}{5}, 'boundary_vess') == 0
        continue
    end
    
    % Get vessel coord
    [row_vess, col_vess] = find(single_vessel_masks{i_vess}{1} == 1);
    [row_vess_fin, col_vess_fin] = find(single_vessel_masks{i_vess}{3} == 1);
    
    % Append them to the arrays
    row_bound_vess = [row_bound_vess; row_vess];
    col_bound_vess = [col_bound_vess; col_vess];
    
    row_bound_vess_fin = [row_bound_vess_fin; row_vess_fin];
    col_bound_vess_fin = [col_bound_vess_fin; col_vess_fin];

end

% Put them in the same array
coord_bound_vess = [row_bound_vess col_bound_vess];

%% Boundary Vessels Angiogenesis

% Every 18 hours angiogenesis only if ratio is below the threshold
if mod(hour, 12) == 0

    % Compute the boundary vessel coverage ratio
    ratio = length(row_bound_vess_fin) / size(coord_tum_bound, 1);

    % Add vessels only if ratio condition is satisfied
    while (ratio < pc3_vess_boundary_coverage)
        
        % Distance between a single tumor boundary site and all boundary vessels
        distance = [];

        % Minimum distances between each tumor boundary site and boundary vessels
        min_distances = [];

        % Array with coord indexes of the closest tum boundary site 
        all_index = [];

        % Re initialize temp_vess variable
        temp_vess = zeros(rows, columns);

        % For loop to find the furthest tumor boundary site from blood vessels
        for i_boun = 1 : size(coord_tum_bound, 1)        
            % For loop over each boundary vess site
            for i_vess = 1 : size(coord_bound_vess, 1)

                % Compute distance
                curr_dist = compute_distance(X, Y, ...
                                             coord_tum_bound(i_boun, 1), ...
                                             coord_tum_bound(i_boun, 2), ...
                                             coord_bound_vess(i_vess, 1), ...
                                             coord_bound_vess(i_vess, 2));

                % Append to the distance variable
                distance = [distance; curr_dist];

            end

            % Get the index of the closest tumor boundary site
            [min_dist, ~] = min(distance);

            % Update the min_distances variable
            min_distances = [min_distances; min_dist];

            % Re initialize the distance variable
            distance = [];

        end

        % Get the index of the furthest tumor_boundary site from vessels
        [~, tum_bound_index] = max(min_distances);

        % Extract corresponding coordinates
        new_row_vess_cent = coord_tum_bound(tum_bound_index, 1);
        new_col_vess_cent = coord_tum_bound(tum_bound_index, 2);
        
        % Add vessel to the bone matrix
        bone(new_row_vess_cent, new_col_vess_cent) = site.vessel;

        % Add vessel to the temp vess variable
        temp_vess(new_row_vess_cent, new_col_vess_cent) = 1;
        
        % Update single_vessel_maks cell
        single_vessel_masks{end+1}{1} = temp_vess;
        % Vessel initialization major axis dim
        single_vessel_masks{end}{2} = site_dim; 

        % Extract major and minor axis from distributions
        major_axis = random(pd_maj_axis_boun);
        minor_axis = max(random(pd_min_axis_boun), 5); % Avoid < 0um axis
        
        % Place the vessel with its final dimensions at the current hour
        for tum_cell = 1 : size(row_tum_bound, 1)

            % I compute the site-vessel center distance
            site_vess_dist = compute_distance(X, Y, ...
                                              row_tum_bound(tum_cell), ...
                                              col_tum_bound(tum_cell), ... 
                                              new_row_vess_cent, ...
                                              new_col_vess_cent);

            % Verify proximity condition
            if site_vess_dist * site_dim <= major_axis / 2
                temp_vess(row_tum_bound(tum_cell), col_tum_bound(tum_cell)) = 1;
            end  

        end
        
        % Add the final vessel size to the cell       
        single_vessel_masks{end}{3} = temp_vess;
        % Add the final major axis size
        single_vessel_masks{end}{4} = major_axis;
        % Define the vessel as boundary vessel
        single_vessel_masks{end}{5} = 'boundary_vess';

        % Update center vessels variable
        new_vessel_prop = [new_row_vess_cent, new_col_vess_cent, site_dim, minor_axis];
        center_vessels  = [center_vessels; new_vessel_prop]; 
        vessels_number = size(center_vessels, 1);

        clear new_vessel_prop
        
        % Variable init
        row_bound_vess = [];
        col_bound_vess = [];

        row_bound_vess_fin = [];
        col_bound_vess_fin = [];

        % I fill the two array with the row and col coords of boundary vessels
        for i_vess = 1 : size(single_vessel_masks, 2)

            % I consider only boundary vessels
            if strcmp(single_vessel_masks{i_vess}{5}, 'boundary_vess') == 0
                continue
            end

            % Get vessel coord
            [row_vess, col_vess] = find(single_vessel_masks{i_vess}{1} == 1);
            [row_vess_fin, col_vess_fin] = find(single_vessel_masks{i_vess}{3} == 1);

            % Append them to the arrays
            row_bound_vess = [row_bound_vess; row_vess];
            col_bound_vess = [col_bound_vess; col_vess];

            row_bound_vess_fin = [row_bound_vess_fin; row_vess_fin];
            col_bound_vess_fin = [col_bound_vess_fin; col_vess_fin];

        end

        % Put them in the same array
        coord_bound_vess = [row_bound_vess col_bound_vess];
        
        % Compute the boundary vessel coverage ratio
        ratio = length(row_bound_vess_fin) / size(coord_tum_bound, 1);
    end
end


%% Central Vessels Angiogenesis

%% Angiogenesis that fills the empty tumor areas

% Every 12 hours, just to fill empty spaces
if mod(hour, 12) == 0
    
    % Re-define the morph operator
    se = strel('disk', 2);
    
    % Erode tumor
    tumor_eroded = imerode(tumor_eroded, se);
    
    % Get tumor central coords
    [row_tum_cent, col_tum_cent] = find(tumor_eroded == 1);
    
    % Generate a random permutation of indices
    perm_indices = randperm(length(row_tum_cent));
    row_tum_cent = row_tum_cent(perm_indices);
    col_tum_cent = col_tum_cent(perm_indices);
    
    % Put them in the same array
    coord_tum_cent = [row_tum_cent col_tum_cent];
    
    % Initialize number of central vessel
    n_vess_sites = 0;
    
    % Define the number of sites occupied by central vessels
    for i_vess = 1 : size(single_vessel_masks, 2)
        
        % Check if vessel belongs to the center
        if strcmp(single_vessel_masks{i_vess}{5}, 'central_ellipse') == 0
            continue
        end
        
        % Get the coordinates of the current vessel
        [row_curr_vess, col_curr_vess] = find(single_vessel_masks{i_vess}{1} == 1);
        % Get the number of sites
        n_vess_sites = n_vess_sites + length(row_curr_vess);
        
    end
    
    % Compute the ratio
    ratio = n_vess_sites / length(row_tum_cent);
    
    % Add just one vessel
    if ratio < pc3_central_vess_coverage

        % Create tumor and vessels mask
        tumor = (bone == site.tumor | bone == site.tumor_edge);
        vessels = bone == site.vessel;

        % Erode tumor
        tumor_eroded = imerode(tumor_eroded, se);

        % Find the most empty area (for vessels)
        [most_empty_quadrant] = vessel_number_per_quadrant(tumor_eroded, vessels);
        [most_empty_quadrant] = vessel_number_per_quadrant(most_empty_quadrant, vessels);

        % Get tumor coordinates of this area
        [row_tum, col_tum] = find(most_empty_quadrant == 1);

        % Generate a random permutation of indices
        perm_indices = randperm(length(row_tum));
        row_tum = row_tum(perm_indices);
        col_tum = col_tum(perm_indices);

        % Get vess center coords
        row_vess_cent = row_tum(1);
        col_vess_cent = col_tum(1);

        % Get Major and minor axis
        major_axis = random(pd_maj_axis_cent_ellipse) / 2;
        minor_axis = random(pd_min_axis_cent_ellipse) / 2;

        % Initialize temp_vess variable
        temp_vess = zeros(rows, columns);

        % Place vessel in the tumor
        for i_tum = 1 : size(coord_tum_cent, 1)

            % Compute tumor vess center distance
            curr_dist = compute_distance(X, Y, row_vess_cent, col_vess_cent, ...
                                         coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2));
            % Get ellipse area
            vess_area = ellipse_radius(ceil(major_axis / site_dim), ...
                                       ceil(minor_axis / site_dim), ...
                                       abs(X(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) - row_vess_cent), ...
                                       abs(Y(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) - col_vess_cent));

           % Check distance condition
           if curr_dist <= vess_area

               % Add vessel site to temp vess variable
               temp_vess(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) = 1; 

           end        

        end

        % At least the vessel center must be shown
        temp_vess(row_vess_cent, col_vess_cent) = 1;

        % Get coordinates of the current vessel
        [row_curr_vess, col_curr_vess] = find(temp_vess == 1);

        % Otherwise upate variable
        single_vessel_masks{end + 1}{1} = zeros(rows, columns);
        single_vessel_masks{end}{2} = major_axis;
        single_vessel_masks{end}{3} = minor_axis;
        single_vessel_masks{end}{4} = '';
        single_vessel_masks{end}{5} = 'central_ellipse';

        % Update center vessels variable
        new_vessel_prop = [row_vess_cent, col_vess_cent, major_axis, minor_axis];
        center_vessels  = [center_vessels; new_vessel_prop]; 
        clear new_vessel_prop
        vessels_number = size(center_vessels, 1);

        % Update vessel variables
        for i_vess = 1 : size(row_curr_vess, 1)

            % Fill the binary masks
            vessels(row_curr_vess(i_vess), col_curr_vess(i_vess)) = 1;
            single_vessel_masks{end}{1}(row_curr_vess(i_vess), col_curr_vess(i_vess)) = 1;
            bone(row_curr_vess(i_vess), col_curr_vess(i_vess)) = site.vessel;

        end

        % Update the number of vess sites
        n_vess_sites = n_vess_sites + length(row_curr_vess);    
    end
end



%% 1st angiogenesis, depending on the in vivo data (sprouting of existing vessels)
% Every 6 hours angiogenesis only if ratio is below the threshold
if mod(hour, 12) == 0
    
    % Generate a random permutation of indices
    perm_indices = randperm(length(row_tum_cent));
    row_tum_cent = row_tum_cent(perm_indices);
    col_tum_cent = col_tum_cent(perm_indices);
    
    % Put them in the same array
    coord_tum_cent = [row_tum_cent col_tum_cent];
    
    % Initialize number of central vessel
    n_vess_sites = 0;
    
    % Define the number of sites occupied by central vessels
    for i_vess = 1 : size(single_vessel_masks, 2)
        
        % Check if vessel belongs to the center
        if strcmp(single_vessel_masks{i_vess}{5}, 'central_ellipse') == 0
            continue
        end
        
        % Get the coordinates of the current vessel
        [row_curr_vess, col_curr_vess] = find(single_vessel_masks{i_vess}{1} == 1);
        % Get the number of sites
        n_vess_sites = n_vess_sites + length(row_curr_vess);
        
    end
    
    % Compute the ratio
    ratio = n_vess_sites / length(row_tum_cent);
    
    % If ratio is below the vessel coverage
    while (ratio < pc3_central_vess_coverage)        
        
        % I'll create a new vessel considering the in vivo distance distribution
        % but it can happen that new coords rest outside tumor
        flag_vess_within_tumor = 0;
        
        % I need to find vess coords within tumor       
        while (flag_vess_within_tumor == 0)
            
            % Variable listing central vessel indexes
            cent_vess_index = [];
            
            % Find current central vessels
            for i_vess = 1 : size(single_vessel_masks, 2)
               if strcmp(single_vessel_masks{i_vess}{5}, 'central_ellipse')
                   cent_vess_index = [cent_vess_index; i_vess];
               end
            end
            
            % Choose a random vessel to perform angiogenesis
            rand_vess_index = randi(size(cent_vess_index, 1));
            % Define its center coords
            start_point = [center_vessels(cent_vess_index(rand_vess_index), 1), ...
                           center_vessels(cent_vess_index(rand_vess_index), 2)];
            % Define the distance from pd (and scale to pixel)
            % 10 is added to 
            vess_distance = (random(pd_distances) * 10) / site_dim; 
            % Define angle
            vess_angle = deg2rad(random(pd_angles_vessels));
            % Get the the new vess center
            new_center = [start_point(1) - ceil(sin(vess_angle) * vess_distance), ...
                          start_point(2) + ceil(cos(vess_angle) * vess_distance)];
                      
            % Extract cordinates
            row_vess_cent = new_center(1);
            col_vess_cent = new_center(2);
            
            % Get axis 
            major_axis = random(pd_maj_axis_cent_ellipse) / 2;
            minor_axis = random(pd_min_axis_cent_ellipse) / 2;
            
            % If new center is inside the tumor center I exit the loop
            if any(all(coord_tum_cent == new_center, 2))
                
                % Raise the flag to one
                flag_vess_within_tumor = 1;
                
            end
        end

        % Initialize temp_vess variable
        temp_vess = zeros(rows, columns);

        % Place vessel in the tumor
        for i_tum = 1 : size(coord_tum_cent, 1)

            % Compute tumor vess center distance
            curr_dist = compute_distance(X, Y, row_vess_cent, col_vess_cent, ...
                                         coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2));
            % Get ellipse area
            vess_area = ellipse_radius(ceil(major_axis / site_dim), ...
                                       ceil(minor_axis / site_dim), ...
                                       abs(X(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) - row_vess_cent), ...
                                       abs(Y(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) - col_vess_cent));

           % Check distance condition
           if curr_dist <= vess_area

               % Add vessel site to temp vess variable
               temp_vess(coord_tum_cent(i_tum, 1), coord_tum_cent(i_tum, 2)) = 1; 

           end        

        end

        % At least the vessel center must be shown
        temp_vess(row_vess_cent, col_vess_cent) = 1;

        % Get coordinates of the current vessel
        [row_curr_vess, col_curr_vess] = find(temp_vess == 1);
        
        % Get current vessel matrix
        vessels = bone == site.vessel;
        
        % Check for possible overlaps
        vess_check = temp_vess + vessels;
        vess_overlap = any(vess_check(:) == 2);

        % If I find vessel overlapping, I go to the next tumor site
        if vess_overlap
            continue
        end
        
        % Otherwise upate variable
        single_vessel_masks{end + 1}{1} = zeros(rows, columns);
        single_vessel_masks{end}{2} = major_axis;
        single_vessel_masks{end}{3} = minor_axis;
        single_vessel_masks{end}{4} = '';
        single_vessel_masks{end}{5} = 'central_ellipse';

        % Update center vessels variable
        new_vessel_prop = [row_vess_cent, col_vess_cent, major_axis, minor_axis];
        center_vessels  = [center_vessels; new_vessel_prop]; 
        clear new_vessel_prop
        vessels_number = size(center_vessels, 1);

        % Update vessel variables
        for i_vess = 1 : size(row_curr_vess, 1)

            % Fill the binary masks
            vessels(row_curr_vess(i_vess), col_curr_vess(i_vess)) = 1;
            single_vessel_masks{end}{1}(row_curr_vess(i_vess), col_curr_vess(i_vess)) = 1;
            bone(row_curr_vess(i_vess), col_curr_vess(i_vess)) = site.vessel;
            
        end
        
        % Update the number of vess sites
        n_vess_sites = n_vess_sites + length(row_curr_vess);

        % Compute ratio between vessel and tumor sites
        ratio = n_vess_sites / size(coord_tum_cent, 1);                
        
    end 
end

