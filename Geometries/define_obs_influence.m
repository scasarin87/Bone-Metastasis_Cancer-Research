%% Define OBs Influence %%

%  This script generates the mask of the PCa cell which will resist to
%  cabozantibib effect influenced by the OBs

clear all

% ABM scaling factor (1 pix = 20.833 um).
site_dim = 500/24; 
% Max Pixel distance which will make cells cabo resistant
pix_dist = 6;
% Cells Influenced by Obs resistance
influenced_cells = pix_dist * site_dim;
% Closest PMit - PApo - Step 
pmit = 0.50;
papo = 0.44;
step = 0.05;

% Desired Geometry
current_geometry = 'femur_1'; % 'long_bone_1'; % long_bone_2, long_bone_3

% Define Directories
codeDirectory = pwd; 
outputDirectory = strcat(codeDirectory, "\Geometries\", current_geometry, "\obs_influence_region");
mitmaskDirectory = strcat(codeDirectory, "\Geometries\", current_geometry, "\pmit_cell_near_obs");
apomaskDirectory = strcat(codeDirectory, "\Geometries\", current_geometry, "\papo_cell_near_obs");

% Load the Desired Bone Geometry
[cortical_bone, bone_marrow, osteoblasts, ~ , ~, ~, ...
    rows, columns] = load_geometry(current_geometry);

% ABM Hexagonal Grid Building
[X, Y, ax, ay, bx, by, ...
    central_col, central_row, directionx, directiony] = hexagonal_grid(rows, columns);

obs_influenced_cells = zeros(rows, columns);
pmit_cell_near_obs = zeros(rows, columns); 
papo_cell_near_obs = 1 * ones(rows, columns); 

% Find X and Y coordinates for both Obs and BM sites
[obs_row, obs_col] = find(osteoblasts == 1);
[bm_row, bm_col] = find(bone_marrow == 1);

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
            if pmit_cell_near_obs(bm_row(bm), bm_col(bm)) < pmit
                pmit_cell_near_obs(bm_row(bm), bm_col(bm)) = pmit + rand(1) * step;
            end
            if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + step
                papo_cell_near_obs(bm_row(bm), bm_col(bm)) = papo + rand(1) * step;   
            end
        end
        % A Probability Gradient is applied
        if distance > (influenced_cells / pix_dist) * 2 && distance <= (influenced_cells / pix_dist) * 4
            if pmit_cell_near_obs(bm_row(bm), bm_col(bm)) < pmit - step
                pmit_cell_near_obs(bm_row(bm), bm_col(bm)) = (pmit - step) + rand(1) * step;
            end
            if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + 2 * step
                papo_cell_near_obs(bm_row(bm), bm_col(bm)) = (papo + step) + rand(1) * step;    
            end
        end
        if distance > (influenced_cells / pix_dist) * 4 && distance <= influenced_cells
            if pmit_cell_near_obs(bm_row(bm), bm_col(bm)) < pmit - 2 * step
                pmit_cell_near_obs(bm_row(bm), bm_col(bm)) = (pmit - 2 * step) + rand(1) * step;
            end
            if papo_cell_near_obs(bm_row(bm), bm_col(bm)) >= papo + 3 * step
                papo_cell_near_obs(bm_row(bm), bm_col(bm)) = (papo + 2 * step) + rand(1) * step;   
            end
        end
        
    end    
end

imagesc(obs_influenced_cells)
colorbar
figure
imagesc(osteoblasts)
figure
imagesc(pmit_cell_near_obs)
colorbar
figure
imagesc(papo_cell_near_obs)
colorbar
video_answer = input("\n\nDo you want to save the current results? \n 1 = yes \n 2 = no\n\nAnswer: ");

if video_answer == 1
    message = sprintf('Saving to: %s', outputDirectory);
    disp(message);
    save(outputDirectory, 'obs_influenced_cells')
    save(mitmaskDirectory, 'pmit_cell_near_obs')
    save(apomaskDirectory, 'papo_cell_near_obs')
elseif video_answer == 2
    disp(["Ok... Your simulation is complete"]);
else
    disp(["Error ...Close the application!"]);
end

