%% Renca Vessel Dynamics %%

% This function is used to:
% 1) Guarantee that renca elongated vessels increase their size up to reaching
%    the major axis value. 
% 2) Place new elongated and elliptical vessels (angiogenesis) up to reaching 
%    the in vivo coverage.


% UN IDEA POTREBBE ESSERE QUELLA DI INIZIARE A GENERARE TUTTO IL VASO IN
% UNA DIREZIONE CHE MI VA MOLTO BENE, E POI FARLO APPARIRE VOLTA PER VOLTA

%% 1) Vessel Growth

% % Definire il vettore x da 0 a 400
% x = 1:1:504;
% 
% % Definire la funzione log-sigmoide manualmente
% logsig = @(x) 1 ./ (1 + exp(-x));
% 
% % Applicare la funzione logsig sul vettore x
% scale_factor = 0.02;  % Fattore di scala per rendere la transizione più smooth
% y = logsig(scale_factor * (x - 252));  % Centrare e scalare la funzione logsig
% 
% % Adattare i valori di y per scalare tra 5 e 20
% y = 6 + 16 * y; % 17 was the original value 14 was the second best
% % 
% % Visualizzare la curva
% hold on
% plot(x(1:360), y(1:360), 'r-', 'LineWidth', 1.5);
% xlabel('x');
% ylabel('y');
% title('Curva con Funzione Logsig e Transizione Smooth');
% grid on;


% 
% if hour < follow_up / 2
%     renca_bv_coverage = 0.17;
% else 
%     renca_bv_coverage = 0.27;
% end

global tumor_dim alpha

if hour <= 504
    renca_bv_coverage = renca_bv_coverage_growth(hour) / 100;
else 
    renca_bv_coverage = renca_bv_coverage_growth(end) / 100;
end

if flag_cabozantinib
    if strcmp(cabo_modulation, 'linear')
        cabo_factor = 1 + (alpha.cabo_max - 1)*(cabo_percentage/100);
    elseif strcmp(cabo_modulation, 'quadratic')
        cabo_factor = 1 + (alpha.cabo_max - 1)*(cabo_percentage/100)^2;
    elseif strcmp(cabo_modulation, 'cubic')
        cabo_factor = 1 + (alpha.cabo_max - 1)*(cabo_percentage/100)^3;
    elseif strcmp(cabo_modulation, 'sigmoid')
        x0 = 0.5;
        k = 10;
        sigmoid = 1 / (1 + exp(-k * ((cabo_percentage/100) - x0)));
        cabo_factor = 1 + sigmoid * (alpha.cabo_max - 1);
    end
    renca_bv_coverage = renca_bv_coverage/cabo_factor;
end

% For loop across all the unique renca vessels generated
for i = 1 : size(single_vessel_masks, 2)
    
    % I need the elongated vessels
    if strcmp(single_vessel_masks{i}{5}, 'elongated') == 0
        continue
    end
    
    % Generate a random number btw 0 and 1
    rand_numb = rand(1);
    
    % The current elongated vessel will growth just the 75% of times
    if rand_numb > 0.5
        continue
    end
    
    % Get the current elongated vessel
    curr_vess = single_vessel_masks{i}{1};

    % Get its coords
    [r_vess, c_vess] = find(curr_vess);
    % Need the number of sites
    n_sites = length(r_vess);
    
    % Get the current major axis
    major_axis = single_vessel_masks{i}{2};
    
    % Get the final major axis
    final_major_axis = single_vessel_masks{i}{4};
    
    % If vessel has already fully growth ... skip
    if major_axis >= final_major_axis
        continue
    end  
    
    % Get the final size blood vessel
    final_vess = single_vessel_masks{i}{3};
    % final_vess = final_vess - curr_vess;
    
    % Get its coordinates
    [row_fin_vess, col_fin_vess] = find(final_vess);
    
    % Get the centroid row and columns
    row_cent = single_vessel_masks{i}{6}(1);
    col_cent = single_vessel_masks{i}{6}(2);
    
    % Matrix for the vess centroid - new vess site distance
    all_distances = [];
    
    % For loop across the coordinates
    for j = 1 : size(row_fin_vess, 1)
        
        % Get current coords
        curr_row = row_fin_vess(j);
        curr_col = col_fin_vess(j);
        
        % Compute the distance
        dist = compute_distance(X, Y, curr_row, curr_col, row_cent, col_cent);
        
        % Update the matrix
        all_distances = [all_distances; dist];          
        
    end

    % Get the minimum distance
    for kk = 1 : n_sites
        [~, min_ind] = min(all_distances);
        all_distances(min_ind) = 10000000000;
    end
    
    [~, min_ind] = min(all_distances);
    
    % It must be a tumor site
    if bone(row_fin_vess(min_ind), col_fin_vess(min_ind)) ~= site.tumor && bone(row_fin_vess(min_ind), col_fin_vess(min_ind)) ~= site.tumor_edge
        continue
    end
    
    % Update the matrixes with the vess growth
    vessels(row_fin_vess(min_ind), col_fin_vess(min_ind)) = 1;
    curr_vess(row_fin_vess(min_ind), col_fin_vess(min_ind)) = 1; 
    bone(row_fin_vess(min_ind), col_fin_vess(min_ind)) = site.vessel;
    single_vessel_masks{i}{1} = curr_vess;
    
    % Update major axis value
    single_vessel_masks{i}{2} = min((major_axis + site_dim), final_major_axis);
    
    
end

%% Angiogenesis - Distance Based and Tumor Growing Based

% Check every 15 hours for angiogenesis
if mod(hour, 3) == 0
    
    % To evaluate the tumor region growing I need to get tumor area on two
    % consecutive time steps
    if hour == 3 
            area_tumor_t0 = sum(sum(bone == site.tumor | bone == site.tumor_edge));
            tumor_mask_t0 = (bone == site.tumor | bone == site.tumor_edge);
    else 
        area_tumor_t1 = sum(sum(bone == site.tumor | bone == site.tumor_edge));
        tumor_mask_t1 = (bone == site.tumor | bone == site.tumor_edge);
    end
    
    area_vessels = 0;
    
    % Control of a possibile vessel over-generation 
    for i_mask = 1 : size(single_vessel_masks, 2)
        area_vessels = area_vessels + sum(sum(single_vessel_masks{i_mask}{1}));
    end
    
    % Find tumor area
    area_tumor = sum(sum(bone == site.tumor | bone == site.tumor_edge));
    
    % Compute the vess tumor ratio
    ratio = area_vessels / area_tumor;

    % Generate new vessels only if there is no overgeneration
    while ratio < renca_bv_coverage
        
        % If < 0.5 I generate a vessel with the distance based criteria -
        % In this code I'll compute the vessel boundaries an avoid them to
        % be placed in a distance < 20-40um. Otherwise I'll boost the
        % vessel growth in the region of tumor expansion
        rand_numb = rand(1);
        
        % Distance based vessel generation
        % If rand < 0.5
        % The first angiogenesis if forced to be a distance based (since I
        % do not have the tumor growth info yet)
        % if tumor size is decreasing there is no tumor growth, so I'll go
        % with the distance based 
        if rand_numb < 0.25 || hour == 3 || area_tumor_t1 <= area_tumor_t0 
        
            % Apply dilation to the vessels mask
            dilate_vess = boundarymask(vessels);
            dilate_vess = imdilate(dilate_vess, strel('disk', 1));

            % Tumor and vessels mask
            tum_vess_mask = bone == site.tumor | bone == site.vessel | bone == site.tumor_edge;

            % Mask difference to get new possible vess sites
            new_vess_mask = tum_vess_mask - dilate_vess;
            % Correct possible errors
            new_vess_mask(new_vess_mask == -1) = 0;

            % Get row and column from the new vess mask
            [row_vess, col_vess] = find(new_vess_mask);
        
        % Tumor growth
        else
            
            % Find the tumor growth region
            tumor_growth_mask = tumor_mask_t1 - tumor_mask_t0;
            % Get the most empty quadrant
            % [tumor_growth_mask] = vessel_number_per_quadrant(tumor_growth_mask, vessels);
            % Get row and column from the new vess mask
            [row_vess, col_vess] = find(tumor_growth_mask  == 1);     

        end
        
        if hour > 3
            % Update tumor growth info
            area_tumor_t0 = area_tumor_t1;
            tumor_mask_t0 = tumor_mask_t1;
        end 
        
        % Shuffle coordinates with a permutation
        rng('shuffle')
        perm_indices = randperm(length(row_vess));
        row_vess = row_vess(perm_indices);
        col_vess = col_vess(perm_indices);
        
        for i = 1 : size(row_vess, 1)
            
            curr_row = row_vess(i);
            curr_col = col_vess(i);
            
            % mask for the final vessel
            temp_vess = zeros(rows, columns);
            temp_vess_final = zeros(rows, columns);
            
            % Initialize just the center in the vessel masks            
            temp_vess(curr_row, curr_col) = 1;
            temp_vess_final = temp_vess;
            
            % Select a random preferential growth direction
            if rand(1) < 0.55
                vess_direction = deg2rad(randi(vess_degree_range));
            else
                vess_direction = deg2rad(randi(vess_degree_range) + 90);
            end
            
            % Randomly select new vessel axis from renca ellipse data
            major_axis = random(pd_maj_axis_elongat);
            vessel_sites = round(major_axis / site_dim);
            minor_axis = random(pd_min_axis_elongat);
            
            % Incremental growth based on the growth direction
            dx = cos(vess_direction);
            dy = sin(vess_direction);
            
            % For loop to generate the vessel
            for j = 0 : round(vessel_sites / 2)
                for k = 0 : 1
                
                    if k == 0
                        % Get the new coords
                        new_vess_row = round(curr_row + j * dx);
                        new_vess_col = round(curr_col + j * dy);

                    elseif k == 1 

                        new_vess_row = round(curr_row - j * dx);
                        new_vess_col = round(curr_col - j * dy);

                    end
                
                    % Update the temp_vess_final mask
                    temp_vess_final(new_vess_row, new_vess_col) = 1;
                end
            end

            temp = zeros(rows, columns);
            for jj = 1 : size(single_vessel_masks, 2)
                temp = temp + single_vessel_masks{jj}{3};
            end 
            
            % Check for overlap
            vess_check = temp_vess_final + temp;             
            vess_overlap = any(vess_check(:) == 2);
            
            % If overlap, try with the next coords
            if vess_overlap
                continue
            % Otherwise, if no overlap, exit the loop
            else
                break
            end
            
        end
            
        % Update single_vessel_masks variable with elongated vess info
        single_vessel_masks{end + 1}{1} = temp_vess;
        single_vessel_masks{end}{2} = site_dim;
        single_vessel_masks{end}{3} = temp_vess_final;
        single_vessel_masks{end}{4} = major_axis;
        single_vessel_masks{end}{5} = 'elongated';  
        single_vessel_masks{end}{6} = [curr_row, curr_col];
        
        % Update vessels mask
        vessels = vessels + temp_vess;
        
        % Update bone variable
        bone(curr_row, curr_col) = site.vessel;
        
        % Update center vessels with the new info
        new_vessel_prop = [curr_row, curr_col, major_axis, minor_axis];
        center_vessels  = [center_vessels; new_vessel_prop]; %center vessels matrix update
        clear new_vessel_prop 
        
        % Re Initialize area vessels
        area_vessels = 0; 
        
        % Control of a possibile vessel over-generation 
        for i_mask = 1 : size(single_vessel_masks, 2)
            area_vessels = area_vessels + sum(sum(single_vessel_masks{i_mask}{1}));
        end 
        
        % Compute the vess tumor ratio
        ratio = area_vessels / area_tumor;
           
    end 
          
end 
% 

% if hour > follow_up / 2 && mod(hour, 9) == 0
%     
%     area_vessels = 0;
%     
%     % Control of a possibile vessel over-generation 
%     for i_mask = 1 : size(single_vessel_masks, 2)
%         area_vessels = area_vessels + sum(sum(single_vessel_masks{i_mask}{3}));
%     end
%     
%     % Find tumor area
%     area_tumor = sum(sum(bone == site.tumor));
%     
%     % Compute the vess tumor ratio
%     ratio = area_vessels / area_tumor;
% 
%     % Generate new vessels only if there is no overgeneration
%     if ratio < renca_bv_coverage + 0.1
%         
%         tumor = bone == site.tumor;
%         
%         % Find the most empty area (for vessels)
%         [most_empty_quadrant] = vessel_number_per_quadrant(tumor, vessels);
%         [most_empty_quadrant] = vessel_number_per_quadrant(most_empty_quadrant, vessels);     
%         
%         % Get row and column from the new vess mask
%         [row_vess, col_vess] = find(most_empty_quadrant);
% 
%         % Shuffle coordinates with a permutation
%         rng('shuffle')
%         perm_indices = randperm(length(row_vess));
%         row_vess = row_vess(perm_indices);
%         col_vess = col_vess(perm_indices);
%         
%         for i = 1 : size(row_vess, 1)
%             
%             curr_row = row_vess(i);
%             curr_col = col_vess(i);
%             
%             % mask for the final vessel
%             temp_vess = zeros(rows, columns);
%             temp_vess_final = zeros(rows, columns);
%             
%             % Initialize just the center in the vessel masks            
%             temp_vess(curr_row, curr_col) = 1;
%             temp_vess_final = temp_vess;
%             
%             % Select a random preferential growth direction
%             if rand(1) < 0.50
%                 vess_direction = deg2rad(randi(vess_degree_range));
%             else
%                 vess_direction = deg2rad(randi(vess_degree_range) + 90);
%             end
%             
%             % Randomly select new vessel axis from renca ellipse data
%             major_axis = random(pd_maj_axis_elongat);
%             vessel_sites = round(major_axis / site_dim);
%             minor_axis = random(pd_min_axis_elongat);
%             
%             % Incremental growth based on the growth direction
%             dx = cos(vess_direction);
%             dy = sin(vess_direction);
%             
%             % For loop to generate the vessel
%             for j = 0 : round(vessel_sites / 2)
%                 for k = 0 : 1
%                 
%                     if k == 0
%                         % Get the new coords
%                         new_vess_row = round(curr_row + j * dx);
%                         new_vess_col = round(curr_col + j * dy);
% 
%                     elseif k == 1 
% 
%                         new_vess_row = round(curr_row - j * dx);
%                         new_vess_col = round(curr_col - j * dy);
% 
%                     end
%                 
%                     % Update the temp_vess_final mask
%                     temp_vess_final(new_vess_row, new_vess_col) = 1;
%                 end
%             end
% 
%             temp = zeros(rows, columns);
%             for jj = 1 : size(single_vessel_masks, 2)
%                 temp = temp + single_vessel_masks{jj}{3};
%             end 
%             
%             % Check for overlap
%             vess_check = temp_vess_final + temp;             
%             vess_overlap = any(vess_check(:) == 2);
%             
%             % If overlap, try with the next coords
%             if vess_overlap
%                 continue
%             % Otherwise, if no overlap, exit the loop
%             else
%                 break
%             end
%             
%         end
%             
%         % Update single_vessel_masks variable with elongated vess info
%         single_vessel_masks{end + 1}{1} = temp_vess;
%         single_vessel_masks{end}{2} = site_dim;
%         single_vessel_masks{end}{3} = temp_vess_final;
%         single_vessel_masks{end}{4} = major_axis;
%         single_vessel_masks{end}{5} = 'elongated';  
%         single_vessel_masks{end}{6} = [curr_row, curr_col];
%         
%         % Update vessels mask
%         vessels = vessels + temp_vess;
%         
%         % Update bone variable
%         bone(curr_row, curr_col) = site.vessel;
%         
%         % Update center vessels with the new info
%         new_vessel_prop = [curr_row, curr_col, major_axis, minor_axis];
%         center_vessels  = [center_vessels; new_vessel_prop]; %center vessels matrix update
%         clear new_vessel_prop 
%         
%         % Re Initialize area vessels
%         area_vessels = 0; 
%         
%         % Control of a possibile vessel over-generation 
%         for i_mask = 1 : size(single_vessel_masks, 2)
%             area_vessels = area_vessels + sum(sum(single_vessel_masks{i_mask}{3}));
%         end 
%         
%         % Compute the vess tumor ratio
%         ratio = area_vessels / area_tumor;
%            
%     end 
%           
% end 
% 


%% 1) Make the elongated vessel grow

% L'idea di questa sezione di codice è la seguente. Per ogni vaso definito
% come 'elongated', nel caso in cui non si sia ancora raggiunta la saturazione
% di vasi per la dimensione corrente del tumore, vado a vedere se ha
% raggiunto la sua dimensione finale (final_major_axis) e, in caso
% contrario, decido il sito in cui andare ad appendere un solo nuovo sito di
% vaso per permettergli mano a mano di raggiungere la completa maturazione.

% Per trovare il sito di espansione prima calcolo il centroide del vaso
% stesso, dopo di che parto dal sito, sempre appartenente al vaso, il più  
% lontano dal centroide e provo a vedere se uno dei suoi neighbour è papabile 
% come nuovo sito. Le condizioni perché ciò avvenga sono: 
% che deve appartenere al tumore, e che non deve avere più di due siti a 
% contatto con il vaso stesso (altrimenti andrei a creare un vaso che si 
% attorciglia su se stesso). Nel caso in cui questo nuovo sito non vada
% bene passo al neighbour successivo. Nel caso in cui nessuno dei neighbour
% vada bene vado ad a cercare tra i neighbour del secondo sito più lontano
% dal centroide, e così via. Se nemmeno in questo caso trovo un sito che
% vada bene aspetto.

% % For loop across all the unique renca vessels generated
% for i = 1 : size(single_vessel_masks, 2)
%     
%     % Check if I had an overdeposition of blood vessels
%     % Compute the current tumor and bv areas
%     curr_tumor_area = sum(sum(bone == site.tumor));
%     curr_vess_area = sum(sum(bone == site.vessel));
%     
%     % Compute the ratio
%     ratio = curr_vess_area / curr_tumor_area;
%     
%     % Check against experimental data
%     if ratio > renca_bv_coverage
%         break
%     end
%     
%     % Initialize flag vess growth
%     flag_vess_growth = 0;
%     
%     % I need the elongated vessels
%     if strcmp(single_vessel_masks{i}{5}, 'elongated') == 0
%         continue
%     end
%     
%     % Generate a random number btw 0 and 1
%     rand_numb = rand(1);
%     
%     % The current elongated vessel will growth just the 50% of times
%     if rand_numb > 0.50
%         continue
%     end
%     
%     % Get the current elongated vessel
%     curr_vess = single_vessel_masks{i}{1};
%     
%     % Get the current major axis
%     major_axis = single_vessel_masks{i}{2};
%     
%     % Get the final major axis
%     final_major_axis = single_vessel_masks{i}{4};
%     
%     % If vessel has already fully growth ... skip
%     if major_axis >= final_major_axis
%         continue
%     end
%     
%     % Get the vessel centroid
%     stats = regionprops(curr_vess, 'Centroid');
%     
%     % Get the centroid row and columns
%     row_cent = round(stats.Centroid(2));
%     col_cent = round(stats.Centroid(1));
%     
%     % Find curr_vess coordinates
%     [row_vess_orig, col_vess_orig] = find(curr_vess);
%     
%     % Shuffle coordinates with a permutation
%     rng('shuffle')
%     perm_indices = randperm(length(row_vess_orig));
%     row_vess = row_vess_orig(perm_indices);
%     col_vess = col_vess_orig(perm_indices);
%     
%     % To keep track of how many tumor neighbour each site has
%     all_distances = [];
%     
%     % For loop across each vessel site
%     for j = 1 : size(row_vess, 1)
%         
%         % Compute vessel site - centroid distance
%         cent_vess_dist = compute_distance(X, Y, row_vess(j), col_vess(j), row_cent, col_cent);
%         
%         % Append variable to the all distance 
%         all_distances = [all_distances; cent_vess_dist];        
%         
%     end
%     
%     [~, max_indexes] = sortrows(all_distances, 'descend');
%     
%     % Sort the vector starting from the furthest point 
%     row_vess = row_vess(max_indexes);
%     col_vess = col_vess(max_indexes);
%     
%     % For loop among the vess sites to check the current site neighbours
%     for j = 1 : size(row_vess, 1)
%         
%         if flag_vess_growth == 1
%             break
%         end
%         
%         % Define a temporary site mask
%         temp_mask_site = zeros(rows, columns);
%         
%         % Fill it with the current site
%         temp_mask_site(row_vess(j), col_vess(j)) = 1;
%         
%         % Get its neighbours
%         dilate_temp_mask_site = boundarymask(temp_mask_site);
%         
%         % Remove the current site with mask subtraction
%         temp_mask_site = dilate_temp_mask_site - temp_mask_site;
%         
%         % Get the coordinates of the neighbours
%         [row_neigh, col_neigh] = find(temp_mask_site);
%         
%         rng('shuffle')
%         perm_indices = randperm(length(row_neigh));
%         row_neigh = row_neigh(perm_indices);
%         col_neigh = col_neigh(perm_indices);
%         
%         % For loop across the neighbours
%         for k = 1 : size(row_neigh, 1)
%             
%             % Check if the neighbours is a tumor site
%             if bone(row_neigh(k), col_neigh(k)) ~= site.tumor
%                 continue
%             end 
%             
%             % Define a second temp_mask to study the neighbour neighbours
%             temp_mask_neigh = zeros(rows, columns);
%             
%             % Fill it with the current neighbour
%             temp_mask_neigh(row_neigh(k), col_neigh(k)) = 1;
% 
%             % Get its neighbours
%             dilate_temp_mask_neigh = boundarymask(temp_mask_neigh);
%             
%             % Remove the current site with mask subtraction
%             temp_mask_neigh = dilate_temp_mask_neigh - temp_mask_neigh;
%             
%             % Get the coordinates of the neighbour neighbours
%             [row_neigh_n, col_neigh_n] = find(temp_mask_neigh);
%             
%             % Counter (if >= 2 skip to the next cell)
%             sum_curr_vess = 0;
%             
%             % Put coords in the same array
%             vess_coord = [row_vess_orig, col_vess_orig];
%             neigh_coord = [row_neigh_n, col_neigh_n];
%             
%             % For loop to find same rows
%             for z1 = 1 : size(neigh_coord, 1)
%                 for z2 = 1 : size(vess_coord, 1)
%                     if isequal(neigh_coord(z1, :), vess_coord(z2, :))
%                         sum_curr_vess = sum_curr_vess + 1;
%                     end
%                 end
%             end
%             
%             % If two or more vessel neighbours are found I need another one 
%             if sum_curr_vess > 1
%                 continue
%             end    
%             
%             % If I found an overlap do not add the vessel in this site
%             if bone(row_neigh(k), col_neigh(k)) == site.vessel
%                 continue
%                        
%             % Otherwise the neighbour site will become a vess
%             else                 
%                 % Update vessels mask
%                 vessels(row_neigh(k), col_neigh(k)) = 1;
%                 
%                 % Update curr_vess mask
%                 curr_vess(row_neigh(k), col_neigh(k)) = 1;
%                 
%                 % Update single_vessel_masks mask and major axis
%                 single_vessel_masks{i}{1} = curr_vess;
%                 single_vessel_masks{i}{2} = min((major_axis + site_dim), final_major_axis);
%                 
%                 % Update bone matrix
%                 bone(row_neigh(k), col_neigh(k)) = site.vessel;
%                 
%                 % Raise a flag
%                 flag_vess_growth = 1;
%                 
%                 % Exit the loop
%                 break
%             end
%             
%         end 
%             
%     end         
%         
% end 

% Check every 18 hours for angiogenesis
% if mod(hour, 600) == 0
%     
%     % Check if I had an overdeposition of blood vessels
%     % Compute the current tumor and bv areas
%     curr_tumor_area = sum(sum(bone == site.tumor));
%     curr_vess_area = sum(sum(bone == site.vessel));
%     
%     % Compute the ratio
%     ratio = curr_vess_area / curr_tumor_area;
%     
%     % Check against experimental data
%     while ratio < renca_bv_coverage
%     
%         % Get tumor area with vessels
%         tumor = (bone == site.tumor) | (bone == site.vessel);
% 
%         % Get the largest tumor portion (some tumor cells might detach from its core)
%         labeled_mask = bwlabel(tumor);
%         props = regionprops(labeled_mask, 'Area');
%         areas = [props.Area];
%         [~, max_id] = max(areas);
%         tumor = labeled_mask == max_id;
% 
%         % Fill tumor area to remove holes, if any
%         tumor = imfill(tumor, 'holes');
% 
%         % Erode tumor edge to place vessels not on the boundary 
%         tumor_eroded = imerode(tumor, strel('disk', 3));
%         
%         % Find the most empty area (for vessels)
%         [most_empty_quadrant] = vessel_number_per_quadrant(tumor_eroded, vessels);
%         [most_empty_quadrant] = vessel_number_per_quadrant(most_empty_quadrant, vessels);
% 
%         % Get tumor coordinates of this area
%         [row_tum, col_tum] = find(most_empty_quadrant == 1);
%         
%         % Generate a random permutation of indices
%         perm_indices = randperm(length(row_tum));
%         row_tum = row_tum(perm_indices);
%         col_tum = col_tum(perm_indices);
% 
%         % Get vess center coords
%         row_vess_cent = row_tum(1);
%         col_vess_cent = col_tum(1);
%         
%         % Initialize temp_vess mask
%         temp_vess = zeros(rows, columns);
%         
%         % Generate random number between 0 and 1
%         rand_numb = rand(1);
%         
%         % if 0 <= rand_numb <= renca_ellips_bv -> Elliptical vessel
%         if rand_numb <= renca_ellipse_bv
%             
%             % Randomly select new vessel semi-axis from renca ellipse data
%             major_axis = random(pd_maj_axis_ellipse) / 2; 
%             minor_axis = random(pd_min_axis_ellipse) / 2;
%             
%             % Place new elliptical RENCA blood vessel
%             for i_vess = 1 : rows
%                  for  j_vess = 1 : columns
%                      
%                      % It must be a tumor site
%                      if bone(i_vess, j_vess) == site.tumor 
%                          
%                          % Compute the vess_center - curr_site distance
%                          dist_centro_vaso = compute_distance(X, Y, row_vess_cent, col_vess_cent, i_vess, j_vess); 
%                          
%                          % Get the vessel elliptical area
%                          vess_area = ellipse_radius(ceil(major_axis / site_dim), ...
%                                                     ceil(minor_axis / site_dim), ...
%                                                     abs(X(i_vess, j_vess) - row_vess_cent), ...
%                                                     abs(Y(i_vess, j_vess) - col_vess_cent));
%                          
%                          % If current point is within vessel area
%                          if dist_centro_vaso <= vess_area
%                              
%                              % Update temporary vessel variable
%                              temp_vess(i_vess, j_vess) = 1;
%                              
%                          end   
%                          
%                      end 
%                      
%                  end  
%             end 
%             
%             % Minimum vessel dimension is 1 pixel
%             temp_vess(row_vess_cent, col_vess_cent) = 1;
%             
%             % Clear loop variables
%             clear i_vess j_vess
%             
%             % Get current vessel matrix
%             vessels = bone == site.vessel;
%             
%             % Check for overlap
%             vess_check = vessels + temp_vess;             
%             vess_overlap = any(vess_check(:) == 2);
% 
%             % If overlap
%             if vess_overlap
%                 % Try again
%                 continue
%             end
%             
%             % Update single_vessel_masks variable with elliptical vess info
%             single_vessel_masks{end + 1}{1} = temp_vess;
%             single_vessel_masks{end}{2} = major_axis;
%             single_vessel_masks{end}{3} = minor_axis;
%             single_vessel_masks{end}{4} = '';
%             single_vessel_masks{end}{5} = 'central_ellipse';
%         
%         % if renca_ellips_bv < rand_numb <= 1 -> Elongated vessel
%         elseif rand_numb > renca_ellipse_bv
%             
%             % Randomly select new vessel axis from renca ellipse data
%             major_axis = random(pd_maj_axis_elongat);
%             minor_axis = random(pd_min_axis_elongat);
%             
%             % Initialize just the center in the vessel masks
%             vessels(row_vess_cent, col_vess_cent) = 1;
%             temp_vess(row_vess_cent, col_vess_cent) = 1;
%             
%             % Update single_vessel_masks variable with elongated vess info
%             single_vessel_masks{end + 1}{1} = temp_vess;
%             single_vessel_masks{end}{2} = site_dim;
%             single_vessel_masks{end}{3} = '';
%             single_vessel_masks{end}{4} = major_axis;
%             single_vessel_masks{end}{5} = 'elongated';            
%             
%             % Update bone matrix
%             bone(row_vess_cent, col_vess_cent) = site.vessel;
%         end 
%         
%         new_vessel_prop = [row_vess_cent, col_vess_cent, major_axis, minor_axis];
%         center_vessels  = [center_vessels; new_vessel_prop]; %center vessels matrix update
%         clear new_vessel_prop  
%         
%         % Check if I had an overdeposition of blood vessels
%         % Compute the current tumor and bv areas
%         curr_tumor_area = sum(sum(bone == site.tumor));
%         curr_vess_area = sum(sum(bone == site.vessel));
% 
%         % Compute the ratio
%         ratio = curr_vess_area / curr_tumor_area;
% 
%         
%     end 
%         
% end         