% Plot in silico and in vivo tumor growth
if strcmp(cell_line, 'c42b')
    [abm_plot_area] = show_area_plot(MEAN, STD, control_data_interpolated, time_points, ...
                                     BONE, follow_up, tumor_area_start, site, show_plot, flag_cabozantinib);
elseif strcmp(cell_line, 'pc3')
    % Plot in silico and in vivo tumor growth
    [abm_plot_area] = show_area_plot(CONTROL_MEDIA, CONTROL_STD, control_curve_pc3, time_points_pc3, ...
                                     BONE, follow_up, tumor_area_start, site, show_plot, flag_cabozantinib);
elseif strcmp(cell_line, 'renca')
    [abm_plot_area] = show_area_plot(CONTROL_MEDIA_RENCA, CONTROL_STD_RENCA, control_curve_renca, time_points_renca, ...
                                     BONE, follow_up, tumor_area_start, site, show_plot, flag_cabozantinib);

end
