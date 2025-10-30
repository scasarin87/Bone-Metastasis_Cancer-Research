%% Plot curves %%

% This function is used to plot in silico results against in vivo data

clear all
close all

% Define cell line
cell_line = 'renca';

% Is the loaded variable a cabo result?
flag_cabo = 1;

% Plot independent curves
flag_ind_curves = 0;

% Load current in silico results
var_name = 'tibia_cs_2_osteolysis_50_sim_renca.mat';
in_silico_results = load(var_name);

c = struct('rr', [0.9047, 0.1918, 0.1988], ...  %Your required color
    'bb', [0.2941, 0.5447, 0.7494], ... %Your required color
    'um', [0.0824, 0.1294, 0.4196], ... %ultra marine
    'lb', [0.1059, 0.9020, 0.9216], ... %light blue
    'lg', [0.1255, 0.9020, 0.2353], ... %light green
    'or', [1.0000, 0.5961, 0.0039], ... %ultra marine
    'ly', [0.9412, 0.9020, 0.0863], ... %ultra marine
    'br', [0.6510, 0.5725, 0.3412], ... %bronze
    'grey', [0.8863,0.8863,0.8863], ... %Your required color
    'gl', [0.8314, 0.7020, 0.7843] );   %greyed lavender

% If the loaded variable is a struct
if isstruct(in_silico_results)
    
    % Get the field names dynamically
    fields = fieldnames(in_silico_results);
    
    % Check if the struct has at least one field
    if ~isempty(fields)
        % Access the value of the first field
        in_silico_results = getfield(in_silico_results, fields{1});
    else
        error('The struct has no fields.');
    end
    
end

% Load in vivo data
if strcmp(cell_line, 'pc3')
    
    % Define PC3 control in vivo data
    mean_control = [1 8.0622 20.768 53.94 69.061 132.573]; 
    std_control = [0 4.141 14.23  21.45 41.47 62.536];
    time_points = [0, 4, 8, 15, 17, 21] * 24;    
    
    % Define PC3 cabo in vivo data
    mean_cabo = [1, 4.752, 12.372, 18.609, 25.219];
    std_cabo = [0, 2.17, 5.765, 6.983, 10.876];   
    time_points_cabo = [0, 5, 8, 12, 15] * 24;
    time_points_cabo(1) = 1;
    
elseif strcmp(cell_line, 'renca')
    
    % Define RENCA control in vivo data
    mean_control = [1, 17.612, 32.674, 141.202];
    std_control = [0, 10.125, 22.602, 81.016];
    time_points = [0, 7, 10, 15] * 24;
    
    % Define RENCA cabo in vivo data
    mean_cabo = [1, 5.422, 9.777, 16.745];
    std_cabo = [0, 2.673, 3.22, 8.765];
    time_points_cabo = [0, 8, 11, 15] * 24;
    time_points_cabo(1) = 1;
    
else
    error('Wrong cell line typed');
end

% Get interpolated curves
control_curve = interp1(time_points, mean_control, 1:1:time_points(end), 'pchip'); 
cabo_curve = interp1(time_points_cabo, mean_cabo, 1:1:time_points_cabo(end), 'pchip');

% Process loaded variable
mean_insilico = mean(in_silico_results);
std_insilico = std(in_silico_results);

% Extract values for the plot
if flag_cabo == 1    
    % Get mean and std
    mean_insilico_tp = mean_insilico(time_points_cabo);
    std_insilico_tp = std_insilico(time_points_cabo); 
    time_points_plot = time_points_cabo;
else    
    % Get mean and std
    mean_insilico_tp = mean_insilico(time_points);
    std_insilico_tp = std_insilico(time_points); 
    time_points_plot = time_points;
end

% Plot data
figure

% In vivo control data
plot(1 : 1 : time_points(end), control_curve, 'r', 'LineWidth', 1.2)
hold on
errorbar(time_points, mean_control, std_control, 'o', 'HandleVisibility', 'off');

% In vivo cabo data
hold on
plot(1 : 1 : time_points_cabo(end), cabo_curve, 'g', 'LineWidth', 1.2)
hold on
errorbar(time_points_cabo, mean_cabo, std_cabo, 'o', 'HandleVisibility', 'off');

% In silico data
hold on
plot(1 : 1 : time_points_plot(end), mean_insilico, 'b', 'LineWidth', 1.2)
hold on
errorbar(time_points_plot, mean_insilico_tp, std_insilico_tp, 'o', 'HandleVisibility', 'off');

% Plot labels
xlabel('Time [hours]');
ylabel('Tumor Size Normalized [ ]');
legend('In vivo control data', 'In vivo cabo data', 'In silico data (av)', '2', 'FontSize', 10, 'Location', 'northwest')

% If cabo data, plot just cabo data
if flag_cabo == 1
    
    % If i have pc3 cabo data, I plot them also on the 15 days follow up
    if strcmp(cell_line, 'pc3')
        
        figure

        % In vivo control data
        plot(1 : 1 : 360, control_curve(1:360), 'r', 'LineWidth', 1.2)
        hold on
        errorbar(time_points(1:4), mean_control(1:4), std_control(1:4), 'o', 'HandleVisibility', 'off');

        % In vivo cabo data
        hold on
        plot(1 : 1 : time_points_cabo(end), cabo_curve, 'g', 'LineWidth', 1.2)
        hold on
        errorbar(time_points_cabo, mean_cabo, std_cabo, 'o', 'HandleVisibility', 'off');

        % In silico data
        hold on
        plot(1 : 1 : time_points_plot(end), mean_insilico, 'b', 'LineWidth', 1.2)
        hold on
        errorbar(time_points_plot, mean_insilico_tp, std_insilico_tp, 'o', 'HandleVisibility', 'off');

        % Plot labels
        xlabel('Time [hours]');
        ylabel('Tumor Size Normalized [ ]');
        legend('In vivo control data', 'In vivo cabo data', 'In silico data (av)', '2', 'FontSize', 10, 'Location', 'northwest')

    end
    
    figure
    
    % In vivo cabo data
    plot(1 : 1 : time_points_cabo(end), cabo_curve, 'g', 'LineWidth', 1.2)
    hold on
    errorbar(time_points_cabo, mean_cabo, std_cabo, 'o', 'HandleVisibility', 'off');

    % In silico data
    hold on
    plot(1 : 1 : time_points_plot(end), mean_insilico, 'b', 'LineWidth', 1.2)
    hold on
    errorbar(time_points_plot, mean_insilico_tp, std_insilico_tp, 'o', 'HandleVisibility', 'off');

    % Plot labels
    xlabel('Time [hours]');
    ylabel('Tumor Size Normalized [ ]');
    legend('In vivo cabo data', 'In silico data (av)', '2', 'FontSize', 10, 'Location', 'northwest')

end


% Plot 10 independent curves
if flag_ind_curves == 1
    
    figure 
    
    for i = 21 : 30
        % Plot current area        
        hold on 
        plot(1 : 1 : time_points_cabo(end), in_silico_results(i, :), 'Color', c.grey, 'LineWidth', 0.08)
    end 
    
    % Plot in vivo cabo
    hold on
    plot(1 : 1 : time_points_cabo(end), cabo_curve, 'g', 'LineWidth', 1.2)
    hold on
    errorbar(time_points_cabo, mean_cabo, std_cabo, 'o', 'HandleVisibility', 'off');    
    
end
