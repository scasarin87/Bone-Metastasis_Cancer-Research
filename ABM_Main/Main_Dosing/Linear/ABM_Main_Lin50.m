%% MAIN %%

%  This is the main function of the ABM.
%  Each epoch simulates the tumor dynamics and bone tissue plasticity under
%  in a defined follow-up time under several conditions:
%  1) Control Regimen: tumor growth without any therapy implementation
%  2) Rad-223 Therapy: tumor response to Rad-223 therapy
%  3) Cabozantinib Th: tumor responde to Cabozantinib therapy

clear variables

% Experimental In-Vivo Control Data
CONTROL_IN_VIVO_DATA;
     
epoch     = 1; % Starting epoch
tot_epoch = 11; % Final epoch (+1)

tumor_area = [];      % N. Rows    -> number of epochs
                      % N. Columns -> normalized tumor area at each follow-up hour      
resorption_area = []; % N. Rows    -> number of epochs
                      % N. Columns -> resortion area in um^2

tic;
total_elapsed_time = 0;

while (epoch < tot_epoch) 
    
    clearvars -except total_elapsed_time all_tum_bm_ratio all_ocs_number all_obs_number bv_number all_resorption_areas cort_bone resorption_area_2weeks renca_lesion_average renca_lesion_std time_point_renca_lesion control_curve_renca time_points_renca tumor_area resorption_area epoch tot_epoch alpha control_data_interpolated time_points MEAN STD CURVACD31 CONTROL_MEDIA CONTROL_MEDIA_RENCA CONTROL_STD CONTROL_STD_RENCA control_curve_pc3 time_points_pc3
    close all
    clc

    epoch_tic = tic;
   
    %% Select Model Features:
    % 1) Geometry: choose between femur_ls_1, femur_ls_1_lr, tibia_ls_1, tibia_ls_1_lr, tibia_cs_1, tibia_cs_2, tibia_cs_3, calvaria_ls_1,vertebra_cs_1
    % 2) CellLine: choose between c42b, pc3, renca
    global current_geometry cell_line CxT CyT tumor_dim n_cells
    current_geometry = 'tibia_ls_1'; 
    cell_line = 'renca';
    
    % Model Driving Parameters
    ABM_DRIVING_PARAMETERS;
    
    % 2) Tumor Dimension Range from [1, 10] and Tumor Center Coordinates
    tumor_dim = 2;
    % if tumor_dim == 0 then it is possible to chose the number of starting
    % cell of the tumor (ranging from 1 to 3)
    n_cells = 0;

    CxT = 63; CyT = 88; % tibia_cs_2 _ edge

    % 3) Start Therapies
    % start_therapy_cabo = 40* 24; % Cabozantinib Therapy starting hour.
    start_therapy_rad  = 40*24; % Rad223 Therapy starting hour.
    start_therapy_za   = 40 * 24;
    end_therapy_cabo   = 40 * 24; % Withdrawal hour of Cabozantinib.
    flag_resorption    = 0;       % 1 means OCs resorption is activates
    flag_cabozantinib  = 1;
    cabo_percentage = 50;
    cabo_modulation = 'linear';
    start_cabo_administration = 0; % Cabozantinib therapy starting hour 
    end_cabo_subministration = 40*24; % Withdrawal hour of Cabozantinib
    rad_periodic = 0; % flag signaling the periodic subministration of rad 
    periodic_timer = 5 * 24; % after how much time rad should be readministrated
    
    % 4) Display ABM and Area Plot
    show_abm  = 0; % if = 1 -> ABM is shown at the end of each iteration hour
    show_plot = 0; % if = 1 -> Tumor Normalized Area is shown every epoch
    
    % 5) control flag for microtumors
    completed = 0;

    while completed == 0
        try
         %% ABM Initialization
            ABM_Initialization;

            hour = 1;
            follow_up = 15*24;

         %% Main Loop: ABM Dynamics Development
            for hour = 1 : follow_up % Time Step: 1 Hours

                % 1) CONTROL SECTION
                if sum(sum(bone == site.tumor | bone == site.tumor_edge)) == 0
                    completed = 1;
                    break
                end

                % 1a) Loop Control: No PCa Cells -> Break the Loop -> Tumor Eradicated
                if isempty(find(bone, 1))
                    disp(['Tumor Eradicated. Hour: ', num2str(hour), '']);
                    break % We can stop the simulation and exit the cycle
                end

                % 1b) Rad223 Start Therapy Control:
                [Rad] = start_rad(Rad, hour, start_therapy_rad);

                if (hour == start_cabo_administration) && (strcmp(cell_line, 'pc3') == 1)
                    flag_cabozantinib = 1;
                    enter_condition = 1;
                end

                % checking whether the cabo therapy must be withdrawled
                if hour == end_cabo_subministration
                    flag_cabozantinib = 0;
                end

                % 1c) Cabozantinib Therapy Control:
                %         if hour == start_therapy_cabo
                %             [flag_cabo, center_vessels, n_vess_start_cabo, ...
                %                 n_vess_eliminated] = start_cabo(flag_cabo, center_vessels);
                %         end

                % if hour == start_therapy_za
                %     delay_factor = 8;
                %     T_ocs_resorption = ceil(T_ocs_resorption * 8);
                %     for row = 1 : rows
                %         for col = 1 : columns
                %             if bone(row, col) == site.osteoclast
                %                 internal_clock_ocs(row, col) = ceil(internal_clock_ocs(row, col) * 8);
                %             end
                %         end
                %     end
                % end

                % 2) MAIN LOOP BODY

                % 2a) Update PCa Cells Internal Clock Every Hour
                [internal_clock, internal_clock_ocs] = update_internal_clock(internal_clock, internal_clock_ocs, rows, ...
                    columns, bone, site, T_cell_division, T_ocs_resorption);

                % 2b) Define Mitosis/Apoptosis Probabilities for each PCa agent
                ABM_Tumor_Dynamics;

                % 2c) Grid Re-Organization ~ Tissue Plasticity
                ABM_Grid_Reorganization;

                if tumor_dim == 0
                    if hour == 3 && strcmp(cell_line, 'renca')
                        area_tumor_t0 = sum(sum(bone == site.tumor | bone == site.tumor_edge));
                        tumor_mask_t0 = (bone == site.tumor | bone == site.tumor_edge);
                    end
                end

                % angiogenesis performed only if the tumor actually has any blood
                % vessel
                if sum(sum(bone== site.tumor | bone == site.tumor_edge)) > 9
                    % if the tumor has blood vessels then angiogenesis
                    if vessels_number > 0
                        % 2d) Angiogenesis Process
                        ABM_Angiogenesis;

                        % 2e) Vessels Elimination if Cabozantinib is Administrated
                        % ABM_VesselsResponseToCabo;

                        % else if the tumor is big enough vessles are implemented
                    else
                        if tumor_dim == 0
                            switch cell_line
                                case 'pc3'
                                    set_pc3_vessels;
                                case 'renca'
                                    set_renca_vessels;
                            end
                        end
                    end

                end

                % 2f) Osteoblasts Reduction when Rad223 is Administered
                 ABM_ObsResponseToRad223;
                 ABM_OcsResponseToRad223;

                % 2g) Tumor mediates OCs bone resorption activity
                if flag_resorption == 1
                    Ocs_activity;
                    Obs_Ocs_Cleaning;
                end

                % 3) VARIABLE UPDATE AND UTILS
                BONE(:,:,hour) = bone;
                MITOTIC(:, :, hour) = mitotic_cells;
                APOPTOTIC(:, :, hour) = apoptotic_cells;

                % Update the temporal dynamic matrix
                cells_matrix(1, hour) = length(find(bone == site.tumor | bone == site.tumor_edge));

                % 3a) Show ABM at the end of every simulation hour
                [abm_matrix] = display_abm(rows, columns, bone, show_abm);

                % Get cort bone data
                cort_bone(epoch, hour) = sum(sum(bone == site.cortical_bone));

                % Number of blood vessels
                bv_number(epoch, hour) = sum(sum(bone == site.vessel));

                % To avoid resorption creating a hole in the bone
                a = bone == site.cortical_bone | bone == site.osteoblast | bone == site.osteoclast;
                a = imfill(a, 'holes');
                a = bwperim(a);
                [row_bound, col_bound] = find(a);
                bone(sub2ind(size(bone), row_bound, col_bound)) = site.cortical_bone;

                curr_ocs_number = sum(sum(bone == site.osteoclast));
                all_ocs_number(epoch, hour) = curr_ocs_number;
                curr_obs_number = sum(sum(bone == site.osteoblast));
                all_obs_number(epoch, hour) = curr_obs_number;

                % 3b) Display the current epoch and simulation hours
                disp(['Epoch: ' num2str(epoch), '', ' - Hour: ', num2str(hour), '']);

                if (rad_periodic == 1) && ((start_therapy_rad + periodic_timer) == (hour + 1))
                    start_therapy_rad = hour + 1;
                end
            end

            % changing the flag if the simulation has reached the end
            if hour == follow_up
                completed = 1;
            end


        catch ME
            % Handle the error by restarting the simulation
            disp(['Error occurred: ', ME.message]);
            disp('Restarting simulation with initial conditions...');
            % The loop restarts from the beginning automatically due to the while ~completed structure
        end

    end

    
    % Plot in vivo vs in silico data over the follow up time
    PlotABM_Area;  
    
    % Plot cabo data
    PlotCabo_Area;
    
    % Compute and plot resorption area
    PlotResorption_Area;
    
    % We update tumor and bone marrow area development for each epoch 
    tumor_area(epoch, :) = abm_plot_area; 
    
    % Generate txt file with cortical bone coordinates 
    % BoneCoordinatesGenerator;

    tumor_at_1 = display_abm(rows, columns, BONE(:,:,1*24), 1);
    name = ['Results/CaboResults/RENCA/Dose/Linear/50/sim',num2str(epoch),'_tumorat1.png'];
    exportgraphics(gcf, name, 'Resolution',1000);
    tumor_at_5 = display_abm(rows, columns, BONE(:,:,5*24), 1);
    name = ['Results/CaboResults/RENCA/Dose/Linear/50/sim',num2str(epoch),'_tumorat5.png'];
    exportgraphics(gcf, name, 'Resolution',1000);
    tumor_at_9 = display_abm(rows, columns, BONE(:,:,9*24), 1);
    name = ['Results/CaboResults/RENCA/Dose/Linear/50/sim',num2str(epoch),'_tumorat9.png'];
    exportgraphics(gcf, name, 'Resolution',1000);
    tumor_at_12 = display_abm(rows, columns, BONE(:,:,12*24), 1);
    name = ['Results/CaboResults/RENCA/Dose/Linear/50/sim',num2str(epoch),'_tumorat12.png'];
    exportgraphics(gcf, name, 'Resolution',1000);
    tumor_at_15 = display_abm(rows, columns, BONE(:,:,15*24), 1);
    name = ['Results/CaboResults/RENCA/Dose/Linear/50/sim',num2str(epoch),'_tumorat15.png'];
    exportgraphics(gcf, name, 'Resolution',1000);

    epoch = epoch + 1; % update the number of epochs

    epoch_elapsed_time = toc(epoch_tic); % Stop timing the epoch
    disp(['Epoch ', num2str(epoch - 1), ' completed in ', num2str(epoch_elapsed_time), ' seconds.']);

    total_elapsed_time = total_elapsed_time + epoch_elapsed_time;
    
end

disp(['Total simulation time: ', num2str(total_elapsed_time), ' seconds.']);

save('Results/CaboResults/RENCA/Dose/Linear/50/medium tumorarea cabo50.mat', 'tumor_area');
save('Results/CaboResults/RENCA/Dose/Linear/50/medium bv number cabo50.mat', 'bv_number');




