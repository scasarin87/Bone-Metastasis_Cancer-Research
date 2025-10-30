%% Vessels Mask

% Input  -> rows, columns   : are x and y size of the processed image matrix
%           tumor           : is the tumor mask (= 1 if site is PCa cell)
%           X, Y            : hexagonal grid sites
%           site_dim        : pixel -> um conversion rate (set to 500/24)
%
% Output -> vessels         : 1 -> vessel site
%                             0 -> other                                                                 
%           center_vessels_0: vessels' center matrix with axis dimension
%           vessels_number  : number of vessels after initialization


function [vessels, center_vessels, vessels_number, single_vessel_masks] = vessels_mask_c42b(rows, columns, tumor, vess_prob, site_dim, X, Y)
    global cell_line CxT CyT ratio_pc3_vess_tum 
    
    if strcmp(cell_line, 'c42b')
        load VesselsProperties_C42B
    elseif strcmp(cell_line, 'pc3')
        load VesselsPropertiesPC3
    else
        error('Wrong cell line selected')
    end
    
    % Define empty output matrix
    vessels = zeros(rows, columns);
    % Each array of this empty vector has the binary mask for each vessel
    single_vessel_masks = cell(1, 0);
    % Define empty variable to store vessels info
    center_vessels = [];
    
    %I compute the number of vessels that have to be created and I randomly select their size and centers  
    Area_tumor        = sum(sum(tumor == 1)) * site_dim * site_dim; %Area in squared um      
    vessels_number = ceil((Mean_Vessel_Density) * Area_tumor); %Number of initial vessels to be created
    mask = zeros(rows, columns);
    
    % C42B Vessels
    if strcmp(cell_line, 'c42b')
        % Poisson disc sampling: The distance between vessels must be at least of a number of pixel determined by spacing variable
        sizeI   = [rows, columns]; 
        vessels = zeros(rows, columns);
        spacing = 3;  
        vess_coords = poissonDisc(sizeI, spacing); % It Determines the poissons points (sampling to define vessels position)
        for i = 1 : size(vess_coords, 1)
            mask(round(vess_coords(i,1)), round(vess_coords(i,2))) = 1;
        end 
        clear i

        %I will draw the centers of the vessels among the points belonging to tumor and being poisson points.
        %I perform a random permutation of the points list to determine the centers of the vessels.
        [xx, yy] = find(tumor == 1 & mask); 
        rng('shuffle')
        n       = numel(xx); 
        ii      = randperm(n);
        xx_rand = xx(ii);
        yy_rand = yy(ii);
        possible_centers = [xx_rand, yy_rand]; %List of possibile center coordinates
        center_vessels   = []; % Each row corresponds to a vessel
                               % 1° and 2° columns: cohordinate of the vessel center 
                               % 3° and 4° columns: Semiaxis of the vessel dimension
        %Cycle that will create the vessels
        for iiii = 1 : vessels_number                

            %Defining randomly vessels center and semiaxis
            %I select the first number of vessels center from NewVesselCenter Flag_Vess_at_edge = 1;

            %Randomly select new vessel center
            Cyy = possible_centers(iiii,1); 
            Cxx = possible_centers(iiii,2);
            %Randomly select new vessel semi-axis from CS data
            a_vessel = random(pd_a_kernels);
            b_vessel = random(pd_b_InverseGaussian);
            
            new_vessel_prop = [Cyy, Cxx, a_vessel, b_vessel];
            center_vessels  = [center_vessels; new_vessel_prop]; %center vessels matrix update
            clear new_vessel_prop  

            %Vessel Mask will be updated with respect to the vessels that are being created                            
            %The vessels mask will be = 1 where a new vessel has been created (0 otherwise)

            for jj = 1 : rows
                 for  kk = 1 : columns
                      if tumor(jj,kk) == 1 
                         dist_centro_vaso = sqrt((X(jj, kk) - Cxx)^2 + (Y(jj, kk) - Cyy)^2); %I compute the site-vessel center distance
                         if dist_centro_vaso <= ellipse_radius(ceil(a_vessel / site_dim), ceil(b_vessel / site_dim), abs(X(jj, kk) - Cxx), abs(Y(jj, kk) - Cyy))
                            vessels(jj, kk) = 1; %Vessels mask
                         end  
                      end 
                 end 
            end

            clear jj kk

            vessels(Cyy, Cxx) = 1; %minimum dimension of a vessel is 1 pixel 

            %Control of a possibile Over-creation of vessels 
            Area_vessels = sum(sum(vessels))*site_dim*site_dim; % in µm^2        
            if(Area_vessels * 15 > Area_tumor) %Parameter of control: Atum>Avess*15 *from experimental images*
               vessels_number = iiii; %I stop iterations if the area criterion is not met
               break 
            end      
        end
        clear iiii
    
    % PC3 Vessels
    else 
        
        % Define the morph operator
        se = strel('disk', 1);

        % Apply erosion
        tumor_eroded = imerode(tumor, se);
        
        % Find circular crown at tumor boundary
        circular_crown = bwperim(tumor_eroded);   
        
        % Find circular crown tumor coordinates 
        [row_tum_cc, col_tum_cc] = find(circular_crown == 1);
        
        % Get the length of the vectors
        n = length(row_tum_cc);
        
        % Find number of PC3 vessel sites according to their average density
        pc3_boundary_vess_sites = round(n * 0.45);
        
        while true
            
            % Exit condition
            % Check if I already added to many vessels
            if sum(sum(vessels == 1)) > pc3_boundary_vess_sites
                vessels_number = size(center_vessels, 1);
                % Exit the while loop
                break                
            end
            
            % Find circular crown tumor coordinates 
            [row_tum_cc, col_tum_cc] = find(circular_crown == 1);
            
            % Get the length of the vectors
            n = length(row_tum_cc);

            % Generate a random permutation of indices
            perm_indices = randperm(n);

            % Shuffle both vectors and take the first value
            row_tum_cc = row_tum_cc(perm_indices);
            row_tum_vess = row_tum_cc(1);
            col_tum_cc = col_tum_cc(perm_indices);
            col_tum_vess = col_tum_cc(1);
            
            % Randomly select the semi-axis
            a_vessel = min(random(pd_b_InverseGaussian), 80); %max(random(pd_a_kernels), 10);
            
            % Define temp vessel variable to check possible overlaps
            temp_vess = zeros(rows, columns);
            
            % Find all vessel sites
            for tum_cell = 1 : size(row_tum_cc, 1)

                % I compute the site-vessel center distance
                dist_centro_vaso = compute_distance(X, Y, ...
                                                    row_tum_cc(tum_cell), ...
                                                    col_tum_cc(tum_cell), ...
                                                    row_tum_vess, ...
                                                    col_tum_vess);
                
                % Verify proximity condition
                if dist_centro_vaso <= a_vessel
                    temp_vess(row_tum_cc(tum_cell), col_tum_cc(tum_cell)) = 1;
                end  
                
            end   
            
            % Check if there was a vessel overlapping
            vess_check = vessels + temp_vess; % current vess + new vessel
            vess_overlap = any(vess_check(:) == 2);
            
            % If vess_overlap = 1 means there was an overlap
            if vess_overlap 
                continue
            end
              
            % If there was no overlap i add the vessel to the main matrix
            vessels = vessels + temp_vess;
            
            
            
            
        end 
        
%         % Find tumor coordinates 
%         [row_tum, col_tum] = find(tumor == 1);
%         
%         % Get the length of the vectors
%         n = length(row_tum);
%         
%         % Find number of PC3 vessel sites according to their average density
%         pc3_vess_sites = round(n * ratio_pc3_vess_tum);
%         
%         % While loop to generate the vessel mask
%         while true
%             
%             % Put the next vessel in the most empty area
%             [most_empty_quadrant] = vessel_number_per_quadrant(tumor, vessels);
%             
%             % Find tumor coordinates of the most empty area
%             [row_tum_area, col_tum_area] = find(most_empty_quadrant == 1);
%             
%             % Get the length of the vectors
%             n = length(row_tum_area);
% 
%             % Generate a random permutation of indices
%             perm_indices = randperm(n);
% 
%             % Shuffle both vectors and take the first value
%             row_tum_area = row_tum_area(perm_indices);
%             row_tum_vess = row_tum_area(1);
%             col_tum_area = col_tum_area(perm_indices);
%             col_tum_vess = col_tum_area(1);
%             
%             % If cpu probability < than prob to find vessel at that site
%             if rand(1) > vess_prob(row_tum_vess, col_tum_vess)
%                 continue
%             end 
%             
%             %Randomly select new vessel semi-axis from CS data
%             a_vessel = max(random(pd_a_kernels), 10);
%             b_vessel = min(random(pd_b_InverseGaussian), 80);
%             
%             % Define temp vessel variable to check possible overlaps
%             temp_vess = zeros(rows, columns);
%             
%             % Find all vessel sites
%             for tum_cell = 1 : size(row_tum, 1)
% 
%                 % I compute the site-vessel center distance
%                 dist_centro_vaso = compute_distance(X, Y, ...
%                                                     row_tum(tum_cell), ...
%                                                     col_tum(tum_cell), ...
%                                                     row_tum_vess, ...
%                                                     col_tum_vess);
%                                                 
%                 % Find ellipse defining vess shape centered in the vessel center
%                 vessel_ellipse_area = ellipse_radius(ceil(a_vessel / site_dim), ...
%                                                      ceil(b_vessel / site_dim), ...
%                                                      abs(X(row_tum(tum_cell), col_tum(tum_cell)) - row_tum_vess), ...
%                                                      abs(Y(row_tum(tum_cell), col_tum(tum_cell)) - col_tum_vess));
%                 
%                 % Verify proximity condition
%                 if dist_centro_vaso <= vessel_ellipse_area
%                     temp_vess(row_tum(tum_cell), col_tum(tum_cell)) = 1;
%                 end  
%                 
%             end   
%             
%             % Check if there was a vessel overlapping
%             vess_check = vessels + temp_vess; % current vess + new vessel
%             vess_overlap = any(vess_check(:) == 2);
%             
%             % If vess_overlap = 1 means there was an overlap
%             if vess_overlap 
%                 continue
%             end
%               
%             % If there was no overlap i add the vessel to the main matrix
%             vessels = vessels + temp_vess;
%             
%             % And to the cell that stores all unique vessel masks
%             single_vessel_masks{end+1} = temp_vess;
%             
%             % New vessel properties
%             new_vessel_prop = [row_tum_vess, col_tum_vess, a_vessel, b_vessel];
%             % Update the main vessels matrix
%             center_vessels  = [center_vessels; new_vessel_prop]; 
%             % Clear temp variable
%             clear new_vessel_prop
%             
%             % Exit condition
%             % Check if I already added to many vessels
%             if sum(sum(vessels == 1)) > pc3_vess_sites
%                 vessels_number = size(center_vessels, 1);
%                 % Exit the while loop
%                 break                
%             end        
%             
%         end 
            
            
  
        
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CODE TO ADD LONGITUDINAL VESSEL TO THE TUMOR %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Or a longitudinal vessel
% add_vess = 0;
% vess_type = 1; % 1 for longitudinal vessels
% %Randomly select new vessel semi-axis from Long data
% a_vessel = min(max(random(pd_a_kernels_par),10), 200);
% b_vessel = max(random(pd_b_InverseGaussian_par), 10);
% % a_vessel must be the major axis
% if a_vessel < b_vessel
%     c_vessel = a_vessel;
%     a_vessel = b_vessel;
%     b_vessel = c_vessel; 
% end
% % Get ellipse radius centered in CxT, CyT passing through row_tum(vess), col_tum(vess)
% radius = compute_distance(X, Y, CxT, CyT, row_tum(vess), col_tum(vess));
% % Expand radius according to vessel dimension
% int_radius = radius - (b_vessel / site_dim);
% ext_radius = radius + (b_vessel / site_dim);
% % Shuffle again tumor indexes
% perm_indices = randperm(n);
% % Shuffle both vectors based on the same permutation
% new_row_tum = row_tum(perm_indices);
% new_col_tum = col_tum(perm_indices);
% % Find all vessel sites
% for tum_cell = 1 : size(new_row_tum, 1)
%     % Distance between tumor center and current vess site
%     dist_centro_tumsite = compute_distance(X, Y, new_row_tum(tum_cell), new_col_tum(tum_cell), CxT, CyT);
%     % Distance between tumor site and the vessel center
%     dist_centrovess_tumsite = compute_distance(X, Y, new_row_tum(tum_cell), new_col_tum(tum_cell), row_tum(vess), col_tum(vess));
%     % I go on only if tum site is between the ellipse boundaries
%     if dist_centro_tumsite >= int_radius && dist_centro_tumsite <= ext_radius
%         % I go on only if tum site is not further than the major semi-axis value
%         if dist_centrovess_tumsite < a_vessel / site_dim
%             add_vess = 1;
%             vessels(new_row_tum(tum_cell), new_col_tum(tum_cell)) = 1;
%             curr_vess(new_row_tum(tum_cell), new_col_tum(tum_cell)) = 1;                                
%             % I put to zero the probability to find choose a new vessel center where it was already added
%             vess_prob(new_row_tum(tum_cell), new_col_tum(tum_cell)) = 0;
%         end
%     end  
% end                         
