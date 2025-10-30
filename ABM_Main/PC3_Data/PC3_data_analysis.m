%% PC3 IN VIVO CONTROL DATA ANALYSIS %%

%  This script analysis the in vivo control data obtained from the PC3
%  cell line tumor growth and extract the tumor growth curve within 15 days

clear all
close all
clc
show_c42b = 0;

%% Read the excel file with the pc3 data

% Define path
filename = 'PC3_data.xlsx';
% Sheet names
sheet_1_name = 'exp1';
sheet_2_name = 'exp2';
% Read excel file
[exp1_data, columns1] = xlsread(filename, sheet_1_name);
[exp2_data, columns2] = xlsread(filename, sheet_2_name);

%% Extract Data

% Extract time points
time_points = exp1_data(:, 1)';
time_points = time_points * 24;

% Extract tumor cell numbers
exp1_data = exp1_data(:, 2 : size(exp1_data, 2));
exp2_data = exp2_data(:, 2 : size(exp2_data, 2));
% Concatenate data into a single matrix
exp_data = [exp1_data, exp2_data];


%% Data analysis

% Get the mean and std of luciferase at every time point
% mean_exp_data = zeros(size(exp_data, 1), 1);
% std_exp_data = zeros(size(exp_data, 1), 1);
% for ind_row = 1 : size(exp_data, 1)    
%     mean_exp_data(ind_row) = mean(exp_data(ind_row, :));    
%     std_exp_data(ind_row) = std(exp_data(ind_row, :));    
% end
% 
% % Normalize wrt the first time point
% for ind_row = 1 : size(exp_data, 2)    
%     exp_data(:, ind_col) = exp_data(:, ind_col) / exp_data(1, ind_col);    
% end


% Normalize data
for ind_col = 1 : size(exp_data, 2)    
    exp_data(:, ind_col) = exp_data(:, ind_col) / exp_data(1, ind_col);    
end

% Get the mean, std, median, and quartiles at the different time points

% Define empty variables
mean_exp_data = zeros(1, size(exp_data, 1));
std_exp_data = zeros(1, size(exp_data, 1));

% Get mean and std
for ind_row = 1 : size(exp_data, 1)
    mean_exp_data(ind_row) = mean(exp_data(ind_row, :));
    std_exp_data(ind_row) = std(exp_data(ind_row, :));
end

% Interpolate data
control_data_interpolated = interp1(time_points, mean_exp_data, 1:1:360, 'pchip'); 

% Plot of the tumor AREA in time
figure
errorbar(time_points, mean_exp_data, std_exp_data, 'o', 'HandleVisibility', 'off')
hold on
plot(1 : 1 : 360, control_data_interpolated, 'LineWidth', 2, 'Color', 'r')
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 1) / exp1_data(1, 1))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 2) / exp1_data(1, 2))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 3) / exp1_data(1, 3))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 4) / exp1_data(1, 4))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 5) / exp1_data(1, 5))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 6) / exp1_data(1, 6))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 7) / exp1_data(1, 7))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 8) / exp1_data(1, 8))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp1_data(:, 9) / exp1_data(1, 9))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 1) / exp1_data(1, 1))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 2) / exp1_data(1, 2))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 3) / exp1_data(1, 3))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 4) / exp1_data(1, 4))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 5) / exp1_data(1, 5))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 6) / exp1_data(1, 6))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 7) / exp1_data(1, 7))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
hold on
plot(1 : 1 : 360, interp1(time_points, (exp2_data(:, 8) / exp1_data(1, 8))', 1:1:360, 'pchip'), 'LineWidth', 0.2)
if show_c42b
    % Average mice tumor growth in c42b cell line 
    MEAN = [1 7.128197054 16.08757084 29.31455403 45.79737765];
    % Standard Deviation at the same time points
    STD = [0 3.719900472 5.890299748 10.73679398 16.58286703];
    time_points = [1, 5 * 24, 8 * 24, 12 * 24, 15 * 24];
    % Experimental Data Interpolated Data
    control_data_interpolated_c42b = interp1(time_points, MEAN, 1:1:360, 'pchip'); 
    hold on
    plot(1 : 1 : 360, control_data_interpolated_c42b, 'LineWidth', 1.5)
    hold on
    errorbar(time_points, MEAN, STD, 'o', 'HandleVisibility', 'off')
end

xlabel('Time [h]')  
ylabel('Normalized Area [ ]')

if show_c42b
    title('In Vivo Control Data Comparison')
    legend('PC3', 'C42B')
else
    title('In Vivo Control Data')
    legend('PC3 Normalized Mean Growth')
end
    






