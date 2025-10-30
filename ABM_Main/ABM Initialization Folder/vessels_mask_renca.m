%% Vessels Mask Renca

% This function initializes both elongated and elliptical renca blood 
% vessels following in vivo data.

% Input  -> rows, columns   : are x and y size of the processed image matrix
%           tumor           : is the tumor mask (= 1 if site is PCa cell)
%           X, Y            : hexagonal grid sites
%           site_dim        : pixel -> um conversion rate (set to 500/24)
%
% Output -> vessels         : is the vessel mask (= 1 if site is vessel)                                                                 
%           center_vessels  : vessels' center matrix with axis dimension
%           vessels_number  : number of vessels after initialization
%           single_vessel_masks: array of cells. Each one has the binary
%                             mask for each vessel and its tag ('boundary',
%                             'central ellipse', 'central elongated')

function [vessels, single_vessel_masks, center_vessels, vessels_number, vess_degree_range] = vessels_mask_renca(rows, columns, tumor, site_dim, X, Y)

    % Load Renca Vessel Data
    load 'ABM Main'/VesselData/VesselsProperties_RENCA.mat
%     
% %     % Definire il vettore x da 0 a 400
%     x = 1:1:follow_up;
%    
%     % Definire la funzione log-sigmoide manualmente
%     logsig2 = @(x) 1 ./ (1 + exp(-x));
% 
%     % Applicare la funzione logsig sul vettore x
%     scale_factor = 0.02;  % Fattore di scala per rendere la transizione più smooth
%     y = logsig2(scale_factor * (x - follow_up / 2));  % Centrare e scalare la funzione logsig
%     y = 6 + 17 * y;
%     y = y / 100;
%     
%     hold on
%     plot(x(1:360), y(1:360))
    
    
    renca_bv_coverage = renca_bv_coverage / 100;
    
    % Choose preferential orentation
    vess_degree = randi(360);
    vess_degree_range = [vess_degree - 15, vess_degree + 15];
    
    % Define empty vessel output matrix
    vessels = zeros(rows, columns);
    
    % Initialize cell that contains each vessel info:
    % For ELONGATED VESSELS
    % Cell 1 -> Binary mask of the current blood vessel
    % Cell 2 -> Size of the current major axis
    % Cell 3 -> empty ''
    % Cell 4 -> Size of the final major axis
    % Cell 5 -> Vessel type (boundary, ellipse, elongated)
    % For ELLIPTICAL VESSELS
    % Cell 1 -> Binary mask of the current blood vessel
    % Cell 2 -> Size of the current major axis
    % Cell 3 -> Size of the current minor axis
    % Cell 4 -> empty ''
    % Cell 5 -> Vessel type (boundary, ellipse, elongated)
    single_vessel_masks = cell(1, 0);
    
    % Define empty variable to store vessels info
    % Each row corresponds to a vessel
    % 1° and 2° columns: vessel center coordinates 
    % 3° and 4° columns: major and minor vessel semi-axis
    center_vessels = [];
    
    % Tumor area in squared um
    area_tumor = sum(sum(tumor == 1)) * site_dim * site_dim;  
    
    % Define the number of starting vessels
    vessels_number = ceil((Mean_Vessel_Density) * area_tumor);
    
    % Temp vessel center mask initialization
    temp_center_mask = zeros(rows, columns);
    
    % Poisson disc sampling: The vessel distance is at least equal to the spacing variable
    size_mask = [rows, columns]; 
    spacing = 2.5;  % between poisson points 
    
    % Poissons points (sampling to define vessels position)
    vess_coords = poissonDisc(size_mask, spacing);
    
    % Fill the temp_center_mask with the Poisson center points
    for i = 1 : size(vess_coords, 1)
        temp_center_mask(round(vess_coords(i,1)), round(vess_coords(i,2))) = 1;
    end 
    
    % Vessel centers are the Poisson points within the tumor mask.
    [row_vess_cen, col_vess_cen] = find(tumor == 1 & temp_center_mask);
    
    % Random permutation of the points list 
    rng('shuffle')
    
    % Define the number of vessels to be generated
    n_vess = numel(row_vess_cen); 
    
    % Get a random index shuffle
    rand_vess_ind = randperm(n_vess);
    
    % Shuffle center rows and columns
    row_rand = row_vess_cen(rand_vess_ind);
    col_rand = col_vess_cen(rand_vess_ind);
    
    % Define a list of possible vessel centers
    possible_centers = [row_rand, col_rand]; 
    
    % Compare the theoretical vess number to the possible centers found
    vessels_number = min(vessels_number, size(possible_centers, 1));
  
    % For loop to generate the vessels 
    for i = 1 : vessels_number      
        
        % Define temporary vessel variable
        temp_vess = zeros(rows, columns);
        
        % Randomly select new vessel center
        row_cen = possible_centers(i, 1); 
        col_cen = possible_centers(i, 2);
        
        % Generate random number between 0 and 1
        rand_numb = rand(1);
        
        % if 0 <= rand_numb <= renca_ellips_bv -> Elliptical vessel
        if rand_numb <= 0 %renca_ellipse_bv
            
            % Randomly select new vessel semi-axis from renca ellipse data
            major_axis = random(pd_maj_axis_ellipse) / 2; 
            minor_axis = random(pd_min_axis_ellipse) / 2;
            
            % Place new elliptical RENCA blood vessel
            for i_vess = 1 : rows
                 for  j_vess = 1 : columns
                     
                     % It must be a tumor site
                     if tumor(i_vess, j_vess) == 1 
                         
                         % Compute the vess_center - curr_site distance
                         dist_centro_vaso = compute_distance(X, Y, row_cen, col_cen, i_vess, j_vess); 
                         
                         % Get the vessel elliptical area
                         vess_area = ellipse_radius(ceil(major_axis / site_dim), ...
                                                    ceil(minor_axis / site_dim), ...
                                                    abs(X(i_vess, j_vess) - row_cen), ...
                                                    abs(Y(i_vess, j_vess) - col_cen));
                         
                         % If current point is within vessel area
                         if dist_centro_vaso <= vess_area
                             
                             % Update vessel variable
                             vessels(i_vess, j_vess) = 1; 
                             
                             % Update temporary vessel variable
                             temp_vess(i_vess, j_vess) = 1;
                         end   
                         
                     end 
                     
                 end  
            end 
            
            % Clear loop variables
            clear i_vess j_vess
            
            % Minimum vessel dimension is 1 pixel
            vessels(row_cen, col_cen) = 1;  
            temp_vess(row_cen, col_cen) = 1;
            
            % Update single_vessel_masks variable with elliptical vess info
            single_vessel_masks{end + 1}{1} = temp_vess;
            single_vessel_masks{end}{2} = major_axis;
            single_vessel_masks{end}{3} = minor_axis;
            single_vessel_masks{end}{4} = '';
            single_vessel_masks{end}{5} = 'central_ellipse';
        
        % If renca_ellips_bv < rand_numb <= 1 -> Elongated vessel
        elseif rand_numb > 0 %renca_ellipse_bv
            
            % mask for the final vessel
            temp_vess_final = zeros(rows, columns);
            
            % Select a random preferential growth direction
            vess_direction = deg2rad(randi(vess_degree_range));
            
            % Randomly select new vessel axis from renca ellipse data
            major_axis = random(pd_maj_axis_elongat);
            vessel_sites = round(major_axis / site_dim);
            minor_axis = random(pd_min_axis_elongat);
            
            % Initialize just the center in the vessel masks
            vessels(row_cen, col_cen) = 1;
            temp_vess(row_cen, col_cen) = 1;
            temp_vess_final = temp_vess;
            
            % Incremental growth based on the growth direction
            dx = cos(vess_direction);
            dy = sin(vess_direction);
            
            % For loop to generate the vessel
            for j = 1 : round(vessel_sites / 2)
                for k = 0 : 1
                
                    if k == 0
                        % Get the new coords
                        new_vess_row = round(row_cen + j * dx);
                        new_vess_col = round(col_cen + j * dy);

                    elseif k == 1 

                        new_vess_row = round(row_cen - j * dx);
                        new_vess_col = round(col_cen - j * dy);

                    end
                
                    % Update the temp_vess_final mask
                    temp_vess_final(new_vess_row, new_vess_col) = 1;
                end
            end
            
            % Update single_vessel_masks variable with elongated vess info
            single_vessel_masks{end + 1}{1} = temp_vess;
            single_vessel_masks{end}{2} = site_dim;
            single_vessel_masks{end}{3} = temp_vess_final;
            single_vessel_masks{end}{4} = major_axis;
            single_vessel_masks{end}{5} = 'elongated';  
            single_vessel_masks{end}{6} = [row_cen, col_cen];
           
        end
        
        new_vessel_prop = [row_cen, col_cen, major_axis, minor_axis];
        center_vessels  = [center_vessels; new_vessel_prop]; %center vessels matrix update
        clear new_vessel_prop  
        
        area_vessels = 0;
        
        % Control of a possibile vessel over-generation 
        for i_mask = 1 : size(single_vessel_masks, 2)
            area_vessels = area_vessels + sum(sum(single_vessel_masks{i_mask}{3}))* site_dim * site_dim;
        end
        
        % Compute the vess tumor ratio
        ratio = area_vessels / area_tumor;
        
        % Check vessel and tumor area
        if ratio > renca_bv_coverage
            break 
        end  
        
    end
    
    % Define vessels number 
    vessels_number = i;
    
    % Clear iterative variable
    clear i    

end