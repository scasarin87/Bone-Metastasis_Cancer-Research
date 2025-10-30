%% ABM Vessels Response To Cabo %%

%  This function is used to define the tumor vessels to be eliminated at 
%  each iteration hour according to the Varkaris' tumor vessels response 
%  curve. 
%  NB: Each cabo targeted vessel, even if its targeted, will influence the
%  events probabilities in a linearly decreasing fashio for the following 
%  4.5 days.

global ratio_pc3_vess_tum_cabo;

% Vessels are eliminated only under Cabo Regimen
if  flag_cabo == 1
    
    % Update the time since vessels were targeted by cabo
    if size(vessels_time, 1) > 0 
            vessels_time(:, 3) = vessels_time(:, 3) + 1;
    end
    
    % Find index of the Vessels Targeted for 4.5 days (vessel_retard) in vessels_time_matrix ...
    vess_coord_1 = find(vessels_time(:, 3) == (vessel_retard) & vessels_time(:, 1) ~= 0); 
    
    % ... If Some Vessel was targeted exactly 4.5 days ago ...
    if size(vess_coord_1, 1) > 0
       for i = 1 : size(vess_coord_1, 1) 
           % ... Identify the corresponding index in the center vessel matrix ...
           vess_coord_2 = find(center_vessels(:, 1) == vessels_time(vess_coord_1(i), 1) & center_vessels(:, 2) == vessels_time(vess_coord_1(i), 2)); 
           
           % ... Get Vessels Info ...
           row_cvess = center_vessels(vess_coord_2(1, 1), 1); % vessel center row
           col_cvess = center_vessels(vess_coord_2(1, 1), 2); % vessel center column
           a_vessel = center_vessels(vess_coord_2(1, 1), 3); % vessel major axis dim
           b_vessel = center_vessels(vess_coord_2(1, 1), 4); % vessel minor axis dim
           
           % ... And Eliminate the Vessel ...
            for jj = 1 : rows
                for  kk = 1 : columns
                    if bone(jj, kk) == site.vessel_cabo
                        dist_centro_vaso = sqrt((X(jj, kk) - col_cvess)^2 + (Y(jj, kk) - row_cvess)^2);
                        if dist_centro_vaso <= ellipse_radius(ceil(a_vessel / site_dim), ceil(b_vessel / site_dim), abs(X(jj, kk) - col_cvess), abs(Y(jj, kk) - row_cvess))
                            bone(jj, kk) = site.bone_marrow;
                        end 
                    end 
                end 
            end            
           
           % Finally we delete the data from the center_vessel matrix
           center_vessels(vess_coord_2(1, 1), :) = []; 
           clear vess_coord_2 row_cvess col_cvess a_vessel b_vessel
       end 
    end  
    clear vess_coord_1
    
    % Vessel number to be eliminated with c42b cell line
    if strcmp(cell_line, 'c42b')
        % Update Time Since Cabo Therapy Started
        time_cabo = hour - start_therapy_cabo + 1; 
        % Update % of Vessels alive from Cabo CD31 curve
        percentual_alive = fix(CURVACD31(time_cabo) * 100);
        % Compute the Vessels Number to be Eliminated at current time
        n_vess_tobe_eliminated = fix((n_vess_start_cabo) * ((100 - percentual_alive) / 100)) - n_vess_eliminated; 
    
    % Vessel number to be eliminated with pc3 cell line
    elseif strcmp(cell_line, 'pc3')
        n_vess_tobe_eliminated = 0;
        if mod(hour, 12) == 0          
            curr_vess_sites = sum(sum(bone == site.vessel));
            curr_tumor_sites = sum(sum(bone == site.tumor));
            curr_tumor_edge_sites = sum(sum(bone == site.tumor_edge));
            ratio = (curr_vess_sites) / (curr_tumor_sites + curr_tumor_edge_sites);
            % If ratio 3% I can remove a vessel, otherwise, exit the loop
            if ratio > ratio_pc3_vess_tum_cabo
                n_vess_tobe_eliminated = 1;
            end 
        end         
    end
    
    % I randomly select the vessel/s to be destroyed
    vess_number = size(center_vessels, 1);
    rng('shuffle')
    random_order = randperm(vess_number); 
    % Random permutation of the Vessels Matrix
    center_vessels2 = center_vessels(random_order, :); 
    center_vessels = center_vessels2;
    clear center_vessels2 
    
    % Check if I need to Eliminate 1 or more vess at the current hour
    if n_vess_tobe_eliminated > 0 
        % For every vessel to be eliminated ...
        for vess = 1 : n_vess_tobe_eliminated
            % ... I check that it hasn't already been targeted by cabo ...
            flag_cabo_control = 0; 
            while flag_cabo_control == 0
                % ... And Randomly Choose the Vessel to be Deleted
                random_number = randi(size(center_vessels, 1)); 
                if center_vessels(random_number, 5) == 0 
                    flag_cabo_control = 1;
                end 
            end 
            
            % Get Vessels Info
            row_cvess = center_vessels(random_number, 1); % vessel center row
            col_cvess = center_vessels(random_number, 2); % vessel center column
            a_vessel = center_vessels(random_number, 3); % vessel major axis dim
            b_vessel = center_vessels(random_number, 4); % vessel minor axis dim
            center_vessels(random_number, 5) = 1; % Raise Cabo Flag for chosen vess

            % Target the Vessel
            for jj = 1 : rows
                for  kk = 1 : columns
                    if bone(jj, kk) == site.vessel
                        dist_centro_vaso = sqrt((X(jj, kk) - col_cvess)^2 + (Y(jj, kk) - row_cvess)^2);
                        if dist_centro_vaso <= ellipse_radius(ceil(a_vessel / site_dim), ceil(b_vessel / site_dim), abs(X(jj, kk) - col_cvess), abs(Y(jj, kk) - row_cvess))
                            bone(jj, kk) = site.vessel_cabo;
                        end 
                    end 
                end 
            end 

           % Vessels_time matrix Update
           vessels_time = [vessels_time; center_vessels(random_number, 1), center_vessels(random_number, 2), 0];

        end   
    end   

    n_vess_eliminated = n_vess_eliminated + n_vess_tobe_eliminated;
    
end  


 