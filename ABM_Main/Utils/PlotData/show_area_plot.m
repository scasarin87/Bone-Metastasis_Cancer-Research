%% Show Area Plot %% 

%  This function creates a plot which compares in vivo experimental tumor 
%  growth under control conditions and our in silico simulation

%  Input  -> MEAN, STD : mean and std values at 5 time points of in vivo
%                        tumor growth
%           control_data_interpolated : interpolated in vivo tumor growth
%           time_point : hours at which the data were retrieved
%           BONE       : contains the bone matrix values for each hour
%           follow_up  : in hour; stands for the in silico follow up time
%           tumor_area_start: tumor area at the first simulation hours
%           site       : struct containing the ABM site identif code
%           show_plot  : if == 1 the area plot is shown by the epoch end
%
%  Output -> abm_plot_area : 

function [abm_plot_area] = show_area_plot(MEAN, STD, control_data_interpolated, time_points, ...
                                              BONE, follow_up, tumor_area_start, site, show_plot, ...
                                              flag_cabozantinib)
                                          
    global cell_line
    
%     tumor_area_start = sum(sum(BONE(:, :, 192) == site.tumor_edge | BONE(:, :, 192) == site.tumor | BONE(:, :, 192) == site.vessel | BONE(:, :, 192) == site.vessel_cabo));
    % Compute normalized area wrt its initial size
    for i = 1 : follow_up 
        abm_plot_area(i) = sum(sum(BONE(:, :, i) == site.tumor_edge | BONE(:, :, i) == site.tumor | BONE(:, :, i) == site.vessel | BONE(:, :, i) == site.vessel_cabo)) / tumor_area_start;
    end
    
    if show_plot

        clear i
        
        % If cabo is administrated we just consider 15 days follow up to match in vivo data
        if strcmp(cell_line, 'pc3') && flag_cabozantinib == 1
            control_data_interpolated = control_data_interpolated(1 : follow_up);
            time_points = time_points(1 : end-2);
            MEAN = MEAN(1 : end-2);
            STD = STD(1 : end-2);
        end
        
        % Plot of the tumor AREA in time
        figure
        errorbar(time_points, MEAN, STD, 'o', 'HandleVisibility', 'off')
        hold on
        plot(1 : 1 : follow_up, control_data_interpolated, 'LineWidth', 1.5)
        hold on
        plot(1 : 1 : follow_up, abm_plot_area, 'LineWidth', 1.5, 'Color', 'blue')
        
        % Set Labels
        xlabel('Time [h]')
        ylabel('Normalized Tumor Area [ ]')
        
        % Set Legend
        legend('Experimental In Vivo Control Data','ABM Simulation Results', '2', 'FontSize',10, 'Location', 'northwest')   
    end   
    
end