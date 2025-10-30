%% PC3 BV Analysis %%

% This script generates the major and minor axis BV distributions to
% implement initialization and angiogenesis of PC3 BVs in the model

% Clear previous data, if any
clear all
close all 
clc

flag_save = 1;

renca_ellipse_bv = 0.3;
renca_elongat_bv = 0.7;

% Define Excel filepath
filepath = 'RENCA_vessels.xlsx';

%% Load Excel file
df_vess_size = readtable(filepath, 'Sheet', 'RENCAVesselSize');
df_vess_dist = readtable(filepath, 'Sheet', 'RENCAVesselDistance');


%% Split between elliptical and elongated central vessels
df_vess_ellipse = contains(df_vess_size.Shape, 'ellipse');
df_vess_ellipse = df_vess_size(df_vess_ellipse, :);
df_vess_elongat = contains(df_vess_size.Shape, 'elongated');
df_vess_elongat = df_vess_size(df_vess_elongat, :);

%% Extract Relevant Data

% Extract Major and Minor Axis from central elliptical vess
maj_axis_ellipse = df_vess_ellipse.MajorAxis;
min_axis_ellipse = df_vess_ellipse.MinorAxis;

% Extract Major and Minor Axis from central elongated vess
maj_axis_elongat = df_vess_elongat.MajorAxis;
min_axis_elongat = df_vess_elongat.MinorAxis;

% Extract PC3 BV reciprocal distances
distances = df_vess_dist.ClosestVess;

%% Clear Useless Data
clear df_vess_size df_vess_dist df_vess_ellipse df_vess_elongat

%% Generate Distributions

% Generate kernel normal distributions for central elliptical vess
pd_maj_axis_ellipse = fitdist(maj_axis_ellipse, 'Normal');
pd_min_axis_ellipse = fitdist(min_axis_ellipse, 'Normal');

% Generate kernel inv gaussian distributions for central elongated vess
pd_maj_axis_elongat = fitdist(maj_axis_elongat, 'InverseGaussian');
pd_min_axis_elongat = fitdist(min_axis_elongat, 'Normal');

% Generate kernel inv gaussian distributions for reciprocal vess distance
pd_distances = fitdist(distances, 'InverseGaussian');

%% Plot  Distributions

% Plot major axis for central elliptic vessels
x_values = linspace(min(maj_axis_ellipse), max(maj_axis_ellipse), 1000);
pdf_values = pdf(pd_maj_axis_ellipse, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Major Axis Elliptical Vessels Probability');

% Plot minor axis for central elliptic vessels
x_values = linspace(min(min_axis_ellipse), max(min_axis_ellipse), 1000);
pdf_values = pdf(pd_min_axis_ellipse, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Minor Axis Elliptical Vessels Probability');

% Plot major axis for central elongated vessels
x_values = linspace(min(maj_axis_elongat), max(maj_axis_elongat), 1000);
pdf_values = pdf(pd_maj_axis_elongat, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Major Axis Elongated Vessels Probability');

% Plot minor axis for central elongated vessels
x_values = linspace(min(min_axis_elongat), max(min_axis_elongat), 1000);
pdf_values = pdf(pd_min_axis_elongat, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Minor Axis Elongated Vessels Probability');

% Plot minor axis for boundary vessels
x_values = linspace(min(distances), max(distances), 1000);
pdf_values = pdf(pd_distances, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Distance [um]');
ylabel('Probability [ ]');
title('RENCA Vessels Distance Distribution');

% Save distributions
if flag_save
    
    % Save in .mat file
    save('VesselsProperties_RENCA.mat', 'pd_maj_axis_ellipse', ...
         'pd_min_axis_ellipse', 'pd_maj_axis_elongat', ...
         'pd_min_axis_elongat', 'pd_distances', ...
         'renca_ellipse_bv', 'renca_elongat_bv')    
end
