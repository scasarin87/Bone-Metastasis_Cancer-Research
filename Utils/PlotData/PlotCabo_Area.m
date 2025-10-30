%% Plot Cabo Area %%

% This function plots the ABM area throughout a follow up time when 
% cabozantinib therapy is administrated against in vivo data
if flag_cabozantinib
    
    % PC3 cell line
    if strcmp(cell_line, 'pc3')
        % Plot in silico and in vivo cabo data
        show_cabo_plot_pc3(abm_plot_area, show_plot, follow_up, control_curve_pc3, MEAN, STD);
    end

    % RENCA cell line
    if strcmp(cell_line, 'renca')
        % Plot in silico and in vivo cabo data
%         show_cabo_plot_renca(abm_plot_area, show_plot, follow_up, control_curve_renca);    
    end
    
end
