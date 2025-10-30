%% Show cabo plot RENCA %%

% This function shows in vivo vs in silico renca tumor progression trend when
% cabozantinib therapy is administered.

function show_cabo_plot_renca(abm_plot_area, show_plot, follow_up, control_curve_renca)
    
    if show_plot 
        
        % In vivo cabo data (x Luca Cabo.xlsx file, sheet_name = RENCA,tumor growth curve)
        time_points = [0, 8, 11, 15] * 24;
        time_points(1) = 1;
        mean_tumor_area_cabo = [1, 5.422, 9.777, 16.745];
        std_tumor_area_cabo = [0, 2.673, 3.22, 8.765];

        % Get pc3 tumor growth trajectory with cabo
        tumor_growth_cabo = interp1(time_points, mean_tumor_area_cabo, 1:1:follow_up, 'pchip');

        % Plot in vivo cabo data
        figure
        plot(1 : 1 : follow_up, tumor_growth_cabo, 'r', 'LineWidth',1.2)           
        hold on
        errorbar(time_points, mean_tumor_area_cabo, std_tumor_area_cabo, 'o','HandleVisibility','off')

        % Plot in silico cabo data
        hold on
        plot(1 : 1 : follow_up, abm_plot_area, 'b', 'LineWidth', 1.2)

        xlabel('Time [hours]');
        ylabel('Tumor Size Normalized [ ]');
        legend('In vivo','In silico (avegrage trend)', '2', 'FontSize',10, 'Location', 'northwest')
        title('In Vivo vs In Silico Tumor Growth with Cabo');
        
    end

end
