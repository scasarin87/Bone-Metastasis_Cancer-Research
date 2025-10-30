%% Show cabo plot PC3 %%

% This function shows in vivo vs in silico pc3 tumor progression trend when
% cabozantinib therapy is administered.

function show_cabo_plot_pc3(abm_plot_area, show_plot, follow_up, control_curve_pc3, MEAN, STD)
    
    if show_plot 
          
%         % OLD IN VIVO DATA
%         % In vivo cabo data
%         time_points = [0, 4, 8, 11, 17, 21] * 24;
%         time_points(1) = 1;
%         mean_tumor_area_cabo = [1, 4.773, 17.654, 25.458, 44.012, 89.601];
%         std_tumor_area_cabo = [0, 3.17, 14.752, 15.874, 20.144, 31.982];

        % NEW IN VIVO DATA
        time_points = [0, 5, 8, 12, 15] * 24;
        time_points(1) = 1;
        mean_tumor_area_cabo = [1, 4.752, 12.372, 18.609, 25.219];
        std_tumor_area_cabo = [0, 2.17, 5.765, 6.983, 10.876];        

        % Get pc3 tumor growth trajectory with cabo
        tumor_growth_cabo = interp1(time_points, mean_tumor_area_cabo, 1:1:follow_up, 'pchip');
        
        %% Plot In Vivo vs In silico Cabo Data

        % Plot in vivo cabo data
        figure
        plot(1 : 1 : follow_up, tumor_growth_cabo, 'r', 'LineWidth',1.2)           
        hold on
        errorbar(time_points, mean_tumor_area_cabo, std_tumor_area_cabo, 'o','HandleVisibility','off')

        % Plot in silico cabo data
        hold on
        plot(1 : 1 : follow_up, abm_plot_area, 'b', 'LineWidth', 1.2)

        %grid on
        xlabel('Time [hours]');
        ylabel('Tumor Size Normalized [ ]');
        legend('In vivo','In silico (avegrage trend)', '2', 'FontSize',10, 'Location', 'northwest')
        title('In Vivo vs In Silico Tumor Growth with Cabo');
        
        %% Plot In Vivo control data against in vivo and in silico cabo data
        % Plot in vivo cabo data
        figure
        plot(1 : 1 : follow_up, tumor_growth_cabo, 'r', 'LineWidth',1.2)           
        hold on
        errorbar(time_points, mean_tumor_area_cabo, std_tumor_area_cabo, 'o','HandleVisibility','off')

        % Plot in silico cabo data
        hold on
        plot(1 : 1 : follow_up, abm_plot_area, 'b', 'LineWidth', 1.2)
        
        hold on
        plot(1 : 1: follow_up, control_curve_pc3(1 : follow_up), 'g', 'LineWidth', 1.2)
        hold on
        errorbar([1, 96, 192], [1, control_curve_pc3(96), control_curve_pc3(192)], STD(1:end-2), 'o','HandleVisibility','off')

        %grid on
        xlabel('Time [hours]');
        ylabel('Tumor Size Normalized [ ]');
        legend('Cabo In vivo','Cabo In silico','Control in vivo', '2', 'FontSize',10, 'Location', 'northwest')
        title('In Vivo vs In Silico Tumor Growth with Cabo');
        
    end

end
