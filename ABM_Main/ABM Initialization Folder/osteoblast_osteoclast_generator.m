%% Osteoblast Osteoclast Generator %%

%  This function computes the internal edge of the cortical bone and then
%  it generates OBs and OCs matrices accordingly. In particular, at first
%  OCs are generated in the area close to the tumor. Then the
%  remaining edge is covered by OBs. Finally, the mitotis/apoptosis
%  probabilities gradients are generated, to implement the OBs induced
%  resistance to cabozantinb.

% Input  -> bone: -2 -> outer site; -1 -> cortical bone
%                  0 -> bone marrow  2 -> tumor site
%                  3 -> vessel site  4 -> near bone tumor site
%                  5 -> OBs site     6 -> OCs site
%           site: struct containing ABM sites names
%           X, Y: hexagonal grid matrices
%           site_dim: ABM scaling factor (20.833 um)
%
% Output -> bone: updated bone matrix, by adding obs and ocs
%           BONE: records bone matrix over follow-up time
%           obs_mit_gradient: mit gradient for obs induced cabo resistance
%           obs_apo_gradient: apo gradient for obs induced cabo resistance

% function [bone, BONE, pmit_cell_near_obs, papo_cell_near_obs, ...
%     obs_influenced_cells, obs_number_start] = osteoblast_osteoclast_generator(bone, site, X, Y, rows, columns, site_dim)      

function [bone, BONE, MITOTIC, APOPTOTIC, papo_cell_near_obs, ...
    obs_influenced_cells, obs_number_start, cortical_bone_edge, no_bone_marrow, method, internal_clock_ocs] = osteoblast_osteoclast_generator(bone, site, X, Y, rows, columns, site_dim, T_ocs_resorption) 
   
    global alpha
    internal_clock_ocs  = zeros(rows, columns);
    % 1) Find CorticalBone-BoneMarrow Edge Coordinates
    % Thanks to https://www.youtube.com/watch?v=5yaj7QWurgM
    method = strel('disk', 1);
    no_bone_marrow = (bone == site.outer | bone == site.cortical_bone); 
    cortical_bone_edge = imerode(no_bone_marrow, method);
    cortical_bone_edge = (no_bone_marrow - cortical_bone_edge);    
    [row_cortbone_edge, col_cortbone_edge] = find(cortical_bone_edge == 1);
    
    % 2) Find Tumor Perimeter Coordinates
    tumor_edge = bwperim(bone == site.tumor | bone == site.vessel);
    [row_tumor_edge, col_tumor_edge] = find(tumor_edge == 1);
    
    % 3) Place Baseline OCs
    % According to in-house experimental data we know that the majority of
    % OCs is located on the bone epiphysis and drastically reduces towards
    % the diaphysis. Moreover, OCs generation is induced by tumor proximity
    
    % Define bone marrow mask 
    bone_marrow = bone == site.bone_marrow;
    % Remove smaller area from the bone marrow mask
    % bone_marrow = keep_largest_region(bone_marrow);
    % Find max and min row of bone marrow mask
    [bm_rows, bm_cols] = find(bone_marrow);
    min_row = min(bm_rows) - 1; 
    max_row = max(bm_rows) + 1;
    avg_row = ceil((min_row + max_row) / 2);
    % Get OCs distribution function
    [ocs_probability, bone_section, bone_type] = get_ocs_probability(rows, min_row, avg_row, max_row);
   
    % Place OCs according to their spatial probability
    for cb_site = 1 : size(row_cortbone_edge,1)
        % Skip iteration if current site is already part of an OCs
        if bone(row_cortbone_edge(cb_site), col_cortbone_edge(cb_site)) == site.osteoclast
            continue
        end
        % Place OCs only if the event probability is matched for current row with longitudinal geometry
        if strncmp(bone_section, 'ls', numel('ls')) && (strncmp(bone_type, 'tibia', numel('tibia')) || strncmp(bone_type, 'femur', numel('femur')))
            if randi(100) < ocs_probability(row_cortbone_edge(cb_site))
                % flag to check if at least one site was already an OCs
                flag_ocs = 0;
                % These lines are used to increase the size of the single OC so that it will occupy from 3 to 5 ABM sites
                osteoclasts = get_osteoclast_matrix(bone, cortical_bone_edge, ...
                    row_cortbone_edge, col_cortbone_edge, cb_site, method);
                [row_ocs, col_ocs] = find(osteoclasts);
                ocs_hour = round(rand(1) * T_ocs_resorption);
                % Check if the n-sites are already occupied by OCs ...
                for oc_site = 1 : size(row_ocs, 1)
                    if bone(row_ocs(oc_site), col_ocs(oc_site)) == site.osteoclast
                        flag_ocs = 1;
                    end
                end
                % ... If so, go to the next site
                if flag_ocs == 1
                    continue
                end
                % ... Otherwis I Place the OC cell and provide the internal clock 
                for oc_site = 1 : size(row_ocs, 1)
                    bone(row_ocs(oc_site), col_ocs(oc_site)) = site.osteoclast;
                    internal_clock_ocs(row_ocs(oc_site), col_ocs(oc_site)) = ocs_hour;                
                end                               
            end
        % Consider Transversal Geometry
        elseif (strncmp(bone_section, 'cs', numel('cs')) || strncmp(bone_type, 'calvaria', numel('calvaria')) || strncmp(bone_type, 'vertebra', numel('vertebra'))) && rand(1) < ocs_probability
            
            % Place only the in vivo ocs amount
            if randi(100) > ocs_probability
                continue
            end            
            % flag to check if at least one site was already an OCs
            flag_ocs = 0;
            % These lines are used to increase the size of the single OC so that it will occupy from 3 to 5 ABM sites
            osteoclasts = get_osteoclast_matrix(bone, cortical_bone_edge, ...
                row_cortbone_edge, col_cortbone_edge, cb_site, method);
            [row_ocs, col_ocs] = find(osteoclasts);
            ocs_hour = round(rand(1) * T_ocs_resorption);
            % Check if the n-sites are already occupied by OCs ...
            for oc_site = 1 : size(row_ocs, 1)
                if bone(row_ocs(oc_site), col_ocs(oc_site)) == site.osteoclast
                    flag_ocs = 1;
                end
            end
            % ... If so, go to the next site
            if flag_ocs == 1
                continue
            end
            % ... Otherwis I Place the OC cell and provide the internal clock 
            for oc_site = 1 : size(row_ocs, 1)
                bone(row_ocs(oc_site), col_ocs(oc_site)) = site.osteoclast;
                internal_clock_ocs(row_ocs(oc_site), col_ocs(oc_site)) = ocs_hour;                
            end     
        end      
    end
           
    % 4) Place OBs
    for ob_site = 1 : size(row_cortbone_edge, 1)
        if bone(row_cortbone_edge(ob_site), col_cortbone_edge(ob_site)) ~= site.osteoclast
            if strncmp(bone_type, 'tibia', numel('tibia'))
                if rand(1) > 0.05 %0.95 is the average OBs distribution in the tibia
                    bone(row_cortbone_edge(ob_site), col_cortbone_edge(ob_site)) = site.osteoblast;
                end
            elseif strncmp(bone_type, 'femur', numel('femur'))
                if rand(1) > 0.20 % 80% is the average OBs distribution in the femur 
                    bone(row_cortbone_edge(ob_site), col_cortbone_edge(ob_site)) = site.osteoblast;
                end
            elseif strncmp(bone_type, 'calvaria', numel('calvaria'))
                if rand(1) > 0.32 % 68% is the average OBs distribution in the calvaria 
                    bone(row_cortbone_edge(ob_site), col_cortbone_edge(ob_site)) = site.osteoblast;
                end
            elseif strncmp(bone_type, 'vertebra', numel('vertebra'))
                if rand(1) > 0.18 % 78% is the average OBs distribution in the vertebra 
                    bone(row_cortbone_edge(ob_site), col_cortbone_edge(ob_site)) = site.osteoblast;
                end
            end
        end            
    end
    
    obs_number_start = sum(sum(bone == site.osteoblast));
    
    % 5) Compute the Mit/Apo Probability Gradient to Make PCa Cell Resist to Cabo 
    
    % Define Costants
    
    % Max Pixel distance which will make cells cabo resistant
    pix_dist = 6;
    % Cells Influenced by Obs resistance
    influenced_cells = pix_dist * site_dim;
    % Closest PMit - PApo - Step 
    papo = alpha.p_apo_min; 
    step = 0.05;
    
    % Define Matrices
    obs_influenced_cells = zeros(rows, columns);
    papo_cell_near_obs = 1 * ones(rows, columns); 

    % Find X and Y coordinates for both Obs and BM sites
    [obs_row, obs_col] = find(bone == site.osteoblast);
    [bm_row, bm_col] = find(bone == site.bone_marrow | bone == site.tumor | bone == site.vessel);

    % Find Distance Btw each OBs site and BM
    for ob = 1 : size(obs_row, 1)
        for bm = 1 : size(bm_row, 1)
            % Compute Ob - Bm Distance
            distance = compute_distance(X, Y, bm_row(bm), bm_col(bm), obs_row(ob), obs_col(ob));
            % Transform Pixel Distance to um Distance
            distance = distance * site_dim;
            % Check if BM site is within the OBs influece region
            if distance <= influenced_cells
                obs_influenced_cells(bm_row(bm), bm_col(bm)) = 1;
            end  
            % Set Fixed Mit/Apo Probabilities for cells within the region
            if distance <= (influenced_cells / pix_dist) * 2
                if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + step
                    papo_cell_near_obs(bm_row(bm), bm_col(bm)) = papo + rand(1) * step;   
                end
            end
            % A Probability Gradient is applied
            if distance > (influenced_cells / pix_dist) * 2 && distance <= (influenced_cells / pix_dist) * 4
                if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + 2 * step
                    papo_cell_near_obs(bm_row(bm), bm_col(bm)) = (papo + step) + rand(1) * step;    
                end
            end
            if distance > (influenced_cells / pix_dist) * 4 && distance <= influenced_cells
                if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + 3 * step
                    papo_cell_near_obs(bm_row(bm), bm_col(bm)) = (papo + 2 * step) + rand(1) * step;   
                end
            end

        end    
    end
    
    % Initialize Matrix to Keep Track of the 'bone' development over time
    BONE(:, :, 1) = bone;
    MITOTIC(:,:,1) = bone;
    APOPTOTIC(:,:,1) = bone;
   
end









    % 3) Place OCs According to tumor distance (https://pubmed.ncbi.nlm.nih.gov/30068572/)
    % Compute CorticalBoneEdge-TumorPerimeter Distance
%     for cb_site = 1 : size(row_cortbone_edge, 1)
%         for tum_site = 1 : size(row_tumor_edge, 1)
%             distance(tum_site) = compute_distance(X, Y, row_tumor_edge(tum_site), col_tumor_edge(tum_site), ...
%                                                         row_cortbone_edge(cb_site), col_cortbone_edge(cb_site));
%         end
%         flag_ocs = 0;
%         % Compute Current CBsite - TumorCell minimum distance
%         cb_tum_mindistance = min(distance);
%         % Compute the probability of generating a new OCs according to the tumor distance from experimental data
%         p_ocs = (0.98623).^cb_tum_mindistance;
%         p_ocs = p_ocs/3.57;
%         % Compare the computed probability with a random generated btw [0,1]
%         if rand(1) < p_ocs
%             % These 5 lines are used to increase the size of the single OC,
%             % so that it will occupy 3 ABM sites at least (60 um)
%             osteoclasts = ones(size(bone, 1), size(bone, 2));
%             osteoclasts(row_cortbone_edge(cb_site), col_cortbone_edge(cb_site)) = 0;
%             osteoclasts = imerode(osteoclasts, method);
%             osteoclasts = imerode(osteoclasts, method);
%             osteoclasts = not(osteoclasts);
%             osteoclasts = (osteoclasts & cortical_bone_edge);
%             [row_ocs, col_ocs] = find(osteoclasts);
%             ocs_hour = round(rand(1) * T_ocs_resorption);
%             % Check if the n-sites are already occupied by ocs 
%             for oc_site = 1 : size(row_ocs, 1)
%                 if bone(row_ocs(oc_site), col_ocs(oc_site)) == site.osteoclast
%                     flag_ocs = 1;
%                 end
%             end
%             if flag_ocs == 1
%                 continue
%             end
%             % Here I Place the OC cell 
%             for oc_site = 1 : size(row_ocs, 1)
%                 bone(row_ocs(oc_site), col_ocs(oc_site)) = site.osteoclast;
%                 internal_clock_ocs(row_ocs(oc_site), col_ocs(oc_site)) = ocs_hour;                
%             end    
%         end      
%     end
    
    % 3A: Internal Clock for OCs
    % internal_clock_ocs  = zeros(rows, columns);
    % Define Random access order to the ABM grid
    % [random_grid_order] = random_grid_access(rows, columns);

    % Now Investigate Each Grid Site Randomly
    % for agent = 1 : (rows * columns)

        % Define Coordinates (row and column) of the Current Site 
    %     [row_agent, col_agent] = current_grid_coordinates(rows, columns, ...
    %                                                       agent, random_grid_order);
        % Check if the Agent is a PCa cell
    %     if bone(row_agent, col_agent) == site.osteoclast
    %         internal_clock_ocs(row_agent, col_agent) = round(rand(1) * T_ocs_resorption); 
    %     end
    % end