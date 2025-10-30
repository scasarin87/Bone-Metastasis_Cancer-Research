%% Change the OBs influence matrix 
%
% this function chnages the influence matrix associated to the OBs and
% their resistance to cabo 
%
% Input -->  bone                              : the bone matrix
%                 obs_influenced_matrix : the binary matrix that point to
%                                                       the cells influenced by Obs
%                 pmit_cell_near_cells    : the matrix that continains the
%                                                       probabilities of mitosis altered by OBs influence
%                 papo_cell_near_obs      : the matrix that continains the
%                                                       probabilities of apoptosis altered by OBs influence
%
% Output --> obs_inflienced_matrix : updated matrix of OBs influence
%                  pmit_cell_near_obs      : updated matrix with probabilites
%                  papo_cell_near_obs      : updated matrix with probabilites
                                          
function [obs_influenced_cells, papo_cell_near_obs] = change_influence_matrix(bone, site, site_dim, X, Y, rows, columns)

    global alpha 

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
    [bm_row, bm_col] = find(bone == site.bone_marrow | bone == site.tumor | bone == site.vessel | bone == site.tumor_edge);

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
    