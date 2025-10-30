%% MAIN %%

%  This is the main function of the ABM.
%  Each epoch simulates the tumor dynamics and bone tissue plasticity under
%  in a defined follow-up time under several conditions:
%  1) Control Regimen: tumor growth without any therapy implementation
%  2) Rad-223 Therapy: tumor response to Rad-223 therapy
%  3) Cabozantinib Th: tumor responde to Cabozantinib therapy

if exist('epoch', 'var') == 0
    clear all 

    % Experimental In-Vivo Control Data
    CONTROL_IN_VIVO_DATA;
     
    epoch     = 0; % Starting epoch
    tot_epoch = 6000; % Final epoch (+1)

    tumor_area = [];      % N. Rows    -> number of epochs
                      % N. Columns -> normalized tumor area at each follow-up hour      
    all_tumor_areas = [];      

end

while (epoch < tot_epoch) 
    
    clearvars -except control_curve_renca time_points_renca tic toc cicle_number mean_mape_array mean_rmse_array all_tumor_areas epoch tot_epoch alpha control_data_interpolated time_points MEAN STD CURVACD31 CONTROL_MEDIA CONTROL_STD control_curve_pc3 time_points_pc3
    close all
       
    %% Select Model Features:
    % 1) Geometry: choose between femur_ls_1, femur_ls_1_lr, tibia_ls_1, tibia_ls_1_lr, tibia_cs_1, tibia_cs_2, tibia_cs_3  
    global current_geometry cell_line CxT CyT tumor_dim
    current_geometry = 'femur_ls_1_lr'; 
    cell_line = 'renca';
    
    % 2) Tumor Dimension Range from [1, 10] and Tumor Center Coordinates
    tumor_dim = 3;
    CxT = 200; CyT = 70;    
    
    % 3) Start Therapies
    start_therapy_cabo = 30 * 24; % Cabozantinib Therapy starting hour.
    start_therapy_rad  = 30 * 24; % Rad223 Therapy starting hour.
    start_therapy_za   = 30 * 24;
    end_therapy_cabo   = 30 * 24; % Withdrawal hour of Cabozantinib.
    flag_resorption    = 0;       % 1 means OCs resorption is activates
    
    show_plot = 0;
    
    end_cicle = 3;
   
    % 3) Define ABM Driving Parameters
    global alpha    
    
    % Update driving parameters every 5 simulations
    if epoch == 0 || mod(epoch, end_cicle) == 0
        
        % Compute the time to process a single set of parameters
        tic 
        
        clear mean_mape mean_rmse 

        % Initialize variables
        alpha1 = 0; alpha2 = 1;
        alpha4 = 0; alpha5 = 1;
        
        % I want Pmitmax > Pmitmin && Papomax > Papomin && Pmitmax + Papomin <= 1
        while (alpha2 > alpha1) || (alpha5 > alpha4) || (alpha1 + alpha5 > 1) 
            rng('shuffle')
            alpha1 = 0.5 + rand(1) * 0.48; %max probability of mitosis between 0.5 and 0.98  
            alpha2 = 0.05 + rand(1) * 0.40; %min probability of mitosis between 0.05 and 0.45
            alpha4 = 0.10 + rand(1) * 0.50; %max probability of apoptosis between 0.10 and 0.60
            alpha5 = 0.01 + rand(1) * 0.29; %min probability of apoptosis between 0.01 and 0.30 
        end  
     
        alpha3 = 0.01 + rand(1) * 0.44; %slope within 0.01 and 0.45
        alpha6 = 110 + rand(1) * 50; %distance between 120 and 200 
        
        alpha = struct(...
            'p_mit_max'     , alpha1, ...  % Maximum probability of mitosis
            'p_mit_min'     , alpha2, ...  % Minimum probability of mitosis
            'slope_prob'    , alpha3, ...  % Tumor Rate Growth
            'p_apo_max'     , alpha4, ...  % Maximum probability of apoptosis
            'p_apo_min'     , alpha5, ...  % Minimum probability of mitosis
            'vess_influence', alpha6);     % Maximum cell-vessel distance mitosis contribution
        
        clear alpha1 alpha2 alpha3 alpha4 alpha5 alpha6
        
        % A cicle re-starts every end_cicle epochs
        cicle_number = 0;
        
        % Define array to store the rmse and mape mean values
        mean_mape_array = [];
        mean_rmse_array = [];
    end 
    
    disp(['Epoch: ', num2str(epoch), ', Mit Max: ', num2str(alpha.p_mit_max), ', Mit Min: ', num2str(alpha.p_mit_min), ...
          ', Slope: ', num2str(alpha.slope_prob), ', Apo Max: ', num2str(alpha.p_apo_max), ...
          ', Apo Min: ', num2str(alpha.p_apo_min), ', Vess Dist: ', num2str(alpha.vess_influence)]);
    
    % 4) Display ABM and Area Plot
    show_abm  = 0; % if = 1 -> ABM is shown at the end of each iteration hour
    
    %% ABM Initialization         
    ABM_Initialization;
    
%     available_area = sum(sum(bone == site.bone_marrow));
    
    %% Main Loop: ABM Dynamics Development
    for hour = 1 : follow_up % Time Step: 1 Hours
                           
        % 2) MAIN LOOP BODY
        
        % 2a) Update PCa Cells Internal Clock Every Hour
        [internal_clock, internal_clock_ocs] = update_internal_clock(internal_clock, internal_clock_ocs, rows, ...
                                                                    columns, bone, site, T_cell_division, T_ocs_resorption);
        
        % 2b) Define Mitosis/Apoptosis Probabilities for each PCa agent
        ABM_Tumor_Dynamics;
        
        % 2c) Grid Re-Organization ~ Tissue Plasticity
        ABM_Grid_Reorganization;
        
        % 2d) Angiogenesis Process
        ABM_Angiogenesis;
                        
        % 3) VARIABLE UPDATE AND UTILS  
      
        % Record the bone variable at each simulation hour
        BONE(:, :, hour) = bone;
        %MITOTIC(:, :, hour) = mitotic_cells;
        %APOPTOTIC(:, :, hour) = apoptotic_cells;
        
        % Update the temporal dynamic matrix
        cells_matrix(1, hour) = length(find(bone == site.tumor | bone == site.tumor_edge));
        
        % 3a) Show ABM at the end of every simulation hour
        [abm_matrix] = display_abm(rows, columns, bone, show_abm);
        
        % Check tumor excessive growth condition
%         actual_tumor_area = sum(sum(bone == site.tumor)) + sum(sum(bone == site.vessel));
%         if actual_tumor_area > 0.75 * available_area
%             % fill arbitrarily BONE variable
%             for remaining_hours = hour : follow_up
%                 BONE(:, :, remaining_hours) = bone;
%             end
%             % Exit the simulation
%             break
%         end
    end 
      
    % Compute mape and rmse for pc3 calibration
    if strcmp(cell_line, 'pc3')
        
        % Get in silico area growth and plot it with in vivo tumor growth
        [tumor_area] = show_area_plot(CONTROL_MEDIA, CONTROL_STD, control_curve_pc3, time_points_pc3, ...
                                      BONE, follow_up, tumor_area_start, site, show_plot);
                                  
        % Compute errors
        mape = 100 * mean(abs(control_curve_pc3 - tumor_area) ./ control_curve_pc3);
        rmse = sqrt(mean((control_curve_pc3 - tumor_area) .^ 2));
    
    % Compute mape and rmse for renca calibration
    elseif strcmp(cell_line, 'renca')
        
        % Get in silico area growth and plot it with in vivo tumor growth
        [tumor_area] = show_area_plot(CONTROL_MEDIA, CONTROL_STD, control_curve_renca, time_points_renca, ...
                                      BONE, follow_up, tumor_area_start, site, show_plot);
                                  
        % Compute errors
        mape = 100 * mean(abs(control_curve_renca - tumor_area) ./ control_curve_renca);
        rmse = sqrt(mean((control_curve_renca - tumor_area) .^ 2));
        
    end
    
    % Collect all tumor areas
    all_tumor_areas = [all_tumor_areas; tumor_area]; 
    if mod(epoch, 50) == 0 && epoch ~= 0
        save("tumor_areas.mat", 'all_tumor_areas');
    end
        
    % Create excel file with error values for each simulation
    output_row = [alpha.p_mit_max, alpha.p_mit_min, alpha.slope_prob, ...
                  alpha.p_apo_max, alpha.p_apo_min, alpha.vess_influence, ...
                  mape, rmse];
    % Convert data to a matlabtable
    output_row = array2table(output_row);
    
    % The first simulation I have to create the excel
    if epoch == 0
        
        % Write the table to an Excel file
        writetable(output_row, 'error_single_table.xlsx');
        
    % Otherwise I just append values to the created one
    else
        
        % Load existing Excel file
        previous_data = readtable('error_single_table.xlsx');
        % Append new data to existing data
        combined_data = [previous_data; output_row];
        % Write combined data back to Excel file
        writetable(combined_data, 'error_single_table.xlsx');
    end
    
    % Append to the mean variables
    mean_mape_array = [mean_mape_array; mape];
    mean_rmse_array = [mean_rmse_array; rmse];
    
    % Increase cicle number variable
    cicle_number = cicle_number + 1;
    
    % Create excel file
    if cicle_number == end_cicle
        
        % Compute the mean of the two errors
        mean_mape = mean(mean_mape_array);
        mean_rmse = mean(mean_rmse_array);
        
        % Create excel file with error values for each simulation
        output_row = [alpha.p_mit_max, alpha.p_mit_min, alpha.slope_prob, alpha.p_apo_max, ...
                      alpha.p_apo_min, alpha.vess_influence, mean_mape, mean_rmse];
        % Convert data to a matlabtable
        output_row = array2table(output_row);
        
        % The first simulation I have to create the excel
        if epoch == end_cicle - 1

            % Write the table to an Excel file
            writetable(output_row, 'error_average_table.xlsx');

        % Otherwise I just append values to the created one
        else

            % Load existing Excel file
            previous_data = readtable('error_average_table.xlsx');
            % Append new data to existing data
            combined_data = [previous_data; output_row];
            % Write combined data back to Excel file
            writetable(combined_data, 'error_average_table.xlsx');
        end
        
        elapsedTime = toc;
        disp(['Elapsed time: ', num2str(elapsedTime), ' seconds']);
        
    end      
        
    epoch = epoch + 1; % update the number of epochs    
        
end

