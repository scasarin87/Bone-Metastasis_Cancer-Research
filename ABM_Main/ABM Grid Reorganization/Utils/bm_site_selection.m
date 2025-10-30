%% Bone Marrow Site Selection %%

%  This function is used to define the bone marrow site which will undergo
%  mitosis or adjacent to the PCa cell that will suffer apoptosis

%  Input  -> marrow_edge       : BM sites adjacent to tumor 
%            row_tum, col_tum  : tumor cell undergoing mit/apo
%            X, Y              : hexagonal grid 
%
%  Output -> row/col mitosis     : chosen site for mitosis/apoptosis event

function [row_mitosis, col_mitosis] = bm_site_selection(marrow_edge, row_tum, col_tum, X, Y, Rad, hour, start_rad, end_hour, flag_cabo)
    
    % Get Tumor Edges Coordinates
    [row_bm, col_bm] = find(marrow_edge == 1);
    
    if size(row_bm, 1) == 0
        
        row_mitosis = [];
        col_mitosis = [];
        
    else
        
    
%%  UNCOMMENT TO GET THE CLOSEST BM SITE    
%     for cell = 1 : size(row_bm, 1)        
%         distance(cell) = compute_distance(X, Y, row_bm(cell), col_bm(cell), row_tum, col_tum);
%     end
%     
%     [~, n_target] = min(distance); 
%     clear distance
% 
%     row_mitosis = row_bm(n_target);
%     col_mitosis = col_bm(n_target);

%%  UNCOMMENT TO GET RANDOM BM SITE AMONG THE N CLOSEST 
    
        % When no therapy is active
        if Rad.activity == 0 
            parameter = 0.25; %0.30 standard for control
        end
        % If rad therapy is on
        if Rad.activity == 1 && flag_cabo == 0 
            parameter = 0.001 + (hour - start_rad) * ((0.30 - 0.001) / (end_hour - start_rad));
            parameter = max(parameter, 0.001);
            parameter = min(parameter, 0.30);
        end
        % If cabo therapy is on
        if flag_cabo == 1 || (Rad.activity == 1 && flag_cabo == 1)
            parameter = 0.01; 
        end

        n_closest_sites = ceil(parameter * size(row_bm, 1));

        for cell = 1 : size(row_bm, 1)        
            distance(cell) = compute_distance(X, Y, row_bm(cell), col_bm(cell), row_tum, col_tum);
        end

        for cell = 1 : n_closest_sites
            [~, n_target(cell)] = min(distance); 
            distance(n_target(cell)) = 500; % Put distance to high
        end 
        clear distance

        n_target = n_target(randperm(length(n_target)));

        row_mitosis = row_bm(n_target(1));
        col_mitosis = col_bm(n_target(1));
    end
 
%%  UNCOMMENT TO GET A RANDOM BM SITE

%     % Attempt Vector Will be Randomized To Access Randomly to the Sites
%      attempt = linspace(1, size(row_bm,  1), size(row_bm,  1));
%      attempt = attempt(randperm(length(attempt)));
%     
%     % Assign Value
%      row_marrow = row_bm(attempt(1));
%      col_marrow = col_bm(attempt(1));


end