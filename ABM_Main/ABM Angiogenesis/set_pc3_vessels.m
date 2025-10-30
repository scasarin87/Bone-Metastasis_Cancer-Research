% Load PC3 Vessel Data
load 'ABM Main'/VesselData/VesselsProperties_PC3.mat

%% Initialization of PC3 vessels on tumor boundaries

% Initialize cell that contains each vessel info:
% Cell 1 -> Binary mask of the current blood vessel
% Cell 2 -> Size of the current major axis
% Cell 3 -> Binary mask of the final size vessel
% Cell 4 -> Size of the final major axis
% Cell 5 -> Vessel type (boundary, ellipse, elongated)
single_vessel_masks = cell(1, 0);
% Initialize temporary matrix for vessels
temp_vess = zeros(rows, columns);
% Initialize vessel matrix
vessels = zeros(rows, columns);
% Initialize the number of vessel sites to 0
curr_vess_sites = 0;
% Overlap flag
flag_ovrl = 0;
% Initialize matrix keeping track of the properties of the vessels
center_vessels = [];
vessels_number = 0;

% vessels are created just if the tumor is larger than a certain number
% of cells
% Define the morph operator
se = strel('disk', 1);

% Apply erosion
tumor_eroded = imerode(tumor, se);

% Find circular crown at tumor boundary
tumor_bound = bwperim(tumor_eroded);

% Get the number of boundary sites
boundary_site = sum(sum(tumor_bound == 1));

% Extract boundary coordinates
[row_tum_bound, col_tum_bound] = find(tumor_bound == 1);

% Get the length of the tumor boundary site vector
n_site_bound = length(row_tum_bound);

% Generate a random permutation of indices
perm_indices = randperm(n_site_bound);
row_tum_bound = row_tum_bound(perm_indices);
col_tum_bound = col_tum_bound(perm_indices);

% For loop to generate boundary vessels
for i_vess = 1 : n_site_bound

    % Check if I've already placed to many vessels
    curr_ratio = curr_vess_sites / boundary_site;

    % For loop exit condition
    if curr_ratio > pc3_vess_boundary_coverage
        break
    end

    % At the first iteration I choose randomly
    if i_vess == 1

        % Define current vess center row and colums
        row_vess_cent = row_tum_bound(i_vess);
        col_vess_cent = col_tum_bound(i_vess);

        % Otherwise I pick a new vessel in the furthest tum boundary site
    else

        % Distance between a single tumor boundary site and all boundary vessels
        distance = [];

        % Minimum distances between each tumor boundary site and boundary vessels
        min_distances = [];

        % Array with coord indexes of the closest tum boundary site
        all_index = [];

        % Re initialize temp_vess variable
        temp_vess = zeros(rows, columns);

        % For loop to find the furthest tumor boundary site from blood vessels
        for i_boun = 1 : size(row_tum_bound, 1)
            % For loop over each boundary vess center
            for i_vess = 1 : size(center_vessels, 1)

                % Compute distance
                curr_dist = compute_distance(X, Y, ...
                    row_tum_bound(i_boun), ...
                    col_tum_bound(i_boun), ...
                    center_vessels(i_vess, 1), ...
                    center_vessels(i_vess, 2));

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

        % Extract corresponsing coordinates
        row_vess_cent = row_tum_bound(tum_bound_index);
        col_vess_cent = col_tum_bound(tum_bound_index);

    end

    % Extract major and minor axis from distributions
    major_axis = min(random(pd_maj_axis_boun), 150); %Avoid too large at init
    minor_axis = max(random(pd_min_axis_boun), 5); % Avoid < 0um axis

    % Add vessel to the temp vess variable
    temp_vess(row_vess_cent, col_vess_cent) = 1;
    % Add to vessels matrix
    vessels(row_vess_cent, col_vess_cent) = 1;
    % Update single_vessel_maks cell
    single_vessel_masks{end+1}{1} = temp_vess;
    % Vessel initialization major axis dim
    single_vessel_masks{end}{2} = site_dim;

    % Find all vessel sites
    for tum_cell = 1 : size(row_tum_bound, 1)

        % I compute the site-vessel center distance
        site_vess_dist = compute_distance(X, Y, ...
            row_tum_bound(tum_cell), ...
            col_tum_bound(tum_cell), ...
            row_vess_cent, ...
            col_vess_cent);

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

    % Update number of vessel sites
    curr_vess_sites = curr_vess_sites + sum(sum(temp_vess == 1));

    % Reinitialize temporary vessel matrix
    temp_vess = zeros(rows, columns);

    % Update center vessels variable
    new_vessel_prop = [row_vess_cent, col_vess_cent, site_dim, minor_axis];
    center_vessels  = [center_vessels; new_vessel_prop];
    vessels_number = size(center_vessels, 1);

    clear new_vessel_prop
end

%% Initialization of PC3 vessel in tumor center
% I erode again the tumor_eroded variable to exclude tumor boundaries
% from the computation of the tumor central area. I have an in vivo
% defined treshold that counts the coverage of central pc3 vessel wrt
% the tumor area. If this ratio < threshold I generate a vessel in the
% space following the distance distribution manually computed.

if sum(sum(tumor)) > 50
    % Re-define the morph operator
    se = strel('disk', 2);
elseif sum(sum(tumor)) > 10
    se = strel('disk', 1);
else
    se = strel('disk', 0);
end

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

% Get current vess center after randomization
row_vess_cent = row_tum_cent(1);
col_vess_cent = col_tum_cent(1);

% Get major and minor axis for the first central vessel
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

% Curr vess sites
n_vess_sites = length(row_curr_vess);

% Compute ratio between vessel and tumor sites
ratio = n_vess_sites / size(coord_tum_cent, 1);

% Put one center_vessel anyway
flag_first_vessel = 1;

% If lower that the n vivo coverage value
while (ratio < pc3_central_vess_coverage || flag_first_vessel == 1)

    % Set to zero this variable
    flag_first_vessel = 0;

    % Init single vessel mask
    % Cell 1 -> Binary mask of the current blood vessel
    % Cell 2 -> Size of the current major axis
    % Cell 3 -> Size of the current minor axis
    % Cell 4 -> empty ''
    % Cell 5 -> Vessel type (boundary, ellipse, elongated)
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

    end

    % Add other vessel/s
    for i_vess = 2 : length(row_tum_cent)

        % Define new vess center
        row_vess_cent = coord_tum_cent(i_vess, 1);
        col_vess_cent = coord_tum_cent(i_vess, 2);

        % Get major and minor axis for the first central vessel
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

        % Check for possible overlaps
        vess_check = temp_vess + vessels;
        vess_overlap = any(vess_check(:) == 2);

        % If I find vessel overlapping, I go to the next tumor site
        if vess_overlap
            continue
        end

        % If I have no overlaps, I exit the loop
        break

    end

    % Get coordinates of the current vessel
    [row_curr_vess, col_curr_vess] = find(temp_vess == 1);

    % Update the number of vess sites
    n_vess_sites = n_vess_sites + length(row_curr_vess);

    % Compute ratio between vessel and tumor sites
    ratio = n_vess_sites / size(coord_tum_cent, 1);
end