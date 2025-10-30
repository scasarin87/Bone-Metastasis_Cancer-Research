%% Ocs_activity%%

% This function is used to implement the resorption activity of the
% cortical bone by the action of osteoclasts

% When ZA therapy start I compute the resorption area from hour=0 to
% hour=start therapy to get the initial reabsorbed area
if hour == 72
    curr_bm_area = sum(sum(bone == site.bone_marrow));
    curr_tumor_area = sum(sum(bone == site.tumor | bone == site.tumor_edge));
    curr_vessel_area = sum(sum(bone == site.vessel | bone == site.vessel_cabo)); 
    total_bm_area = curr_bm_area + curr_tumor_area + curr_vessel_area;
    res_area_3_days = total_bm_area - bone_marrow_start_area; 
end

% Initialized the cortical bone edge interface
cortical_bone_edge_new = cortical_bone_edge * (-1);
% Find tumor perimeter  
tumor_edge = bwperim(bone == site.tumor | bone == site.tumor_edge | bone == site.vessel | bone == site.vessel_cabo);
% Get tumor perimeter coordinates
[row_tumor_edge, col_tumor_edge] = find(tumor_edge == 1);    
% Get Ocs coordinates
[row_ocs, col_ocs] = find(bone == site.osteoclast);
for ocs = 1:length(row_ocs)
    if bone(row_ocs(ocs), col_ocs(ocs)) == site.osteoclast  
        if internal_clock_ocs(row_ocs(ocs), col_ocs(ocs)) == T_ocs_resorption
            internal_clock_ocs(row_ocs(ocs), col_ocs(ocs)) = 0;
            % Investigate the neighborhood
            if mod(col_ocs(ocs), 2) == 0
                j_k = 2;
            else
                j_k = 1;
            end        
            % Inizialization of matrix liste
            liste = zeros(4, 6);
            % look for neighbors
            for n = 1 : 6
                % (j3,k3) --> coordinates of the temporary (j,k) neighbor
                j3 = row_ocs(ocs)+directionx(n,j_k);
                k3 = col_ocs(ocs)+directiony(n);
                % stop the propagation of tumor once out of the bone
                if bone(j3,k3) ~= site.outer                    
                    % liste(1,:) --> neighbors of (j,k) indexes (1 if it is cortical bone, 0 if not)
                    % liste(2,:) --> row of current position
                    % liste(3,:) --> column of current position
                    % liste(4,:) --> (1 if it is not part of the cortical bone edge interface) 
                    if (bone(j3,k3) == site.cortical_bone) || (bone(j3,k3) == site.cortical_bone_induced)
                        liste(1,n) = 1;
                        liste(2,n) = j3;
                        liste(3,n) = k3;
                        if bone(j3,k3) ~= cortical_bone_edge_new(j3,k3)
                            liste(4,n) = 1;
                        end 
                    end 
                end 
            end 

            % Random choose the site of resorption ensuring the stochasticity of the model
            new_liste = zeros(4,6);        
            rng('shuffle');
            new_liste(:,:)= liste(:, randperm(6));
            
            % Compute distance btw tumor and each ocs
            distance = [];
            for tum_site = 1 : size(row_tumor_edge, 1)
                distance(tum_site) = compute_distance(X, Y, row_tumor_edge(tum_site), col_tumor_edge(tum_site), ...
                                                        row_ocs(ocs), col_ocs(ocs));
            end

            % Defining the minimum distance nedeed to initiate the resorption activity
            % Resorption occurs when the distance between tumor and ocs site is less than 
            % 1 site and only cortical bone is resorbed. We also impose the condition that 
            % this site of cortical bone must not belong to the edge to ensure the propagation 
            % toward external sites
            cb_tum_mindistance = min(distance);
            if cb_tum_mindistance < 1.12
                for n2 = 1:6
                    if new_liste(1,n2) == 1 && new_liste(4,n2) == 1
                        bone(new_liste(2,n2),new_liste(3,n2)) = site.osteoclast;
                        internal_clock_ocs(new_liste(2,n2),new_liste(3,n2)) = 0;
                        bone(row_ocs(ocs),col_ocs(ocs)) = site.bone_marrow;
                        break
                    end  
                end 
            end 


        end  
    end
end

% Update the cortical bone edge after the resorption
no_bone_marrow = (bone == site.outer | bone == site.cortical_bone | bone == site.cortical_bone_induced | bone == site.osteoblast | bone == site.osteoclast);
cortical_bone_edge = imerode(no_bone_marrow, method);
cortical_bone_edge = (no_bone_marrow - cortical_bone_edge);

% Ensure the presence of Ocs at CB-Tumor interface
[row_int,col_int] = find(cortical_bone_edge == 1 & bone == site.cortical_bone);
for int = 1:length(row_int)
    bone(row_int(int),col_int(int)) = site.osteoclast;
end

% If pc3 tumor is close to the cortical bone
if strcmp(cell_line, 'pc3') || strcmp(cell_line, 'renca')
            
    % Find cortical bone sites
    [row_cb, col_cb] = find(cortical_bone_edge);
       
    % Get tumor perim sites
    tumor = (bone == site.tumor) | (bone == site.vessel) | (bone == site.tumor_edge);
    tumor_perim = bwperim(tumor);
    [row_tp, col_tp] = find(tumor_perim);
    
    % Get distance between cb and tumor
    for site_cb = 1 : size(row_cb, 1)
        for site_tp = 1 : size(row_tp, 1)
            
            % Don't do anything if it is already an OC
            if bone(row_cb(site_cb), col_cb(site_cb)) == site.osteoclast
                continue
            end                
            
            % Get distance
            curr_dist = compute_distance(X, Y, row_cb(site_cb), ...
                col_cb(site_cb), row_tp(site_tp), col_tp(site_tp)) * site_dim;
            
            % I just need distances below site_dim
            if curr_dist > site_dim
                continue
            end
            
            % Give stocasticity
            rand_numb = rand(1);
            if rand_numb > 0.1
                continue
            end
            
            bone(row_cb(site_cb), col_cb(site_cb)) = site.osteoclast;                                      
        
        end
    end
        
    
end




 
