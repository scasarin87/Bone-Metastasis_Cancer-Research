%% PC3 VESSEL ANALYSIS %%

%  This script analysis the in vivo vessel data obtained from the PC3
%  cell line. The difference between the probability distribution chosen 
%  for the c42b and pc3 is that for both major and minor pc3 axis I chose
%  the InverseGaussian proability density, while for the c42b the major
%  axis are modeled with a KernelDistribution. The reason is that pc3 data
%  are hugely right skewed for both axis, and an InverseGaussian fits
%  perfectly with this data unbalancement.

clear all
close all
clc

show_plot = 1;
compare_with_c42b = 1;
flag_save = 1;

%% Read the excel file with the pc3 data

% Define path
filename = 'PC3_vessels.xlsx';
% Sheet names
sheet_1_name = 'perp';
% sheet_2_name = 'parall';
% Read excel file
[vess_data_perp, columns_perp] = xlsread(filename, sheet_1_name);
% [vess_data_paral, columns_paral] = xlsread(filename, sheet_2_name);

%% Extract Data

% Extract major and minor axis
major_axis = vess_data_perp(:, 2);
minor_axis = vess_data_perp(:, 3);

%% Major Axis 

% Fit a kernel density distribution
pd_major_axis = fitdist(major_axis, 'InverseGaussian');

% Show distribution plots
if show_plot
    % Generate values from the fitted distribution for visualization
    x_values = linspace(0, max(major_axis), 1000);
    y_values = pdf(pd_major_axis, x_values);

    % Plot the fitted kernel density distribution
    figure;
    plot(x_values, y_values, 'LineWidth', 0.5, 'Color', 'r');
    title('Major Axis');
    xlabel('Dimension [um]');
    ylabel('Probability');
end

% Plot the comparison with c42b distribution
if compare_with_c42b
    load('VesselsProperties_C42B.mat')
    y_values_c42b = pdf(pd_a_kernels, x_values);
    hold on
    plot(x_values, y_values_c42b, 'LineWidth', 0.5, 'Color', 'b');
    legend('pc3', 'c42b')
end

%% Minor Axis 

% Fit a kernel density distribution
pd_minor_axis = fitdist(minor_axis, 'InverseGaussian');

% Show distribution plots
if show_plot
    % Generate values from the fitted distribution for visualization
    x_values = linspace(0, max(minor_axis), 1000);
    y_values = pdf(pd_minor_axis, x_values);

    % Plot the fitted kernel density distribution
    figure;
    plot(x_values, y_values, 'LineWidth', 0.5, 'Color', 'r');
    title('Minor Axis');
    xlabel('Dimension [um]');
    ylabel('Probability');
end

% Plot the comparison with c42b distribution
if compare_with_c42b
    y_values_c42b = pdf(pd_b_InverseGaussian, x_values);
    hold on
    plot(x_values, y_values_c42b, 'LineWidth', 0.5, 'Color', 'b');
    legend('pc3', 'c42b')
end

% Save distributions
if flag_save
    % I need the name of the PC3 variables to be equal to the C42B ones
    pd_a_kernels = pd_major_axis;
    pd_b_InverseGaussian = pd_minor_axis;
    % Save in .mat file
    save('VesselsProperties_PC3_new.mat', 'pd_a_kernels', ...
         'pd_b_InverseGaussian', 'Mean_Vessel_Density', ...
         'pd_angles_vessels')    
end


