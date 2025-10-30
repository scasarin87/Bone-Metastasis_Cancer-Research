%% PC3 BV Analysis %%

% This script generates the major and minor axis BV distributions to
% implement initialization and angiogenesis of PC3 BVs in the model

% Clear previous data, if any
clear all
close all 
clc

flag_save = 1;

pc3_ellipse_bv = 0.6;
pc3_elongat_bv = 0.4;

% Define Excel filepath
filepath = 'PC3_vessels_cabo.xlsx';

%% Load Excel file
df_vess_covr = readtable(filepath, 'Sheet', 'PC3VessBoundary');
df_vess_cent = readtable(filepath, 'Sheet', 'PC3VessCenter');
df_vess_boun = readtable(filepath, 'Sheet', 'PC3VessBound');
df_vess_dist = readtable(filepath, 'Sheet', 'PC3VessDistances');

%% Define PC3 Vessel Tumor Boundary Coverage
ratios = df_vess_covr.Ratio;
pc3_vess_boundary_coverage = mean(ratios);

%% Split between elliptical and elongated central vessels
df_vess_cent_ellipse = contains(df_vess_cent.Shape, 'ellipse');
df_vess_cent_ellipse = df_vess_cent(df_vess_cent_ellipse, :);
df_vess_cent_elongat = contains(df_vess_cent.Shape, 'elongated');
df_vess_cent_elongat = df_vess_cent(df_vess_cent_elongat, :);

%% Extract Relevant Data

% Extract Major and Minor Axis from central elliptical vess
maj_axis_cent_ellipse = df_vess_cent_ellipse.MajorAxis;
min_axis_cent_ellipse = df_vess_cent_ellipse.MinorAxis;

% Extract Major and Minor Axis from central elongated vess
maj_axis_cent_elongat = df_vess_cent_elongat.MajorAxis;
min_axis_cent_elongat = df_vess_cent_elongat.MinorAxis;

% Extract Major and Minor Axis from boundary vess
maj_axis_boun = df_vess_boun.MajorAxis;
min_axis_boun = df_vess_boun.MinorAxis;

% Extract PC3 BV reciprocal distances
distances = df_vess_dist.ClosestVess;

%% Clear Useless Data
clear df_vess_cent df_vess_boun df_vess_dist df_vess_cent_ellipse df_vess_cent_elongat

%% Generate Distributions

% Generate kernel normal distributions for central elliptical vess
pd_maj_axis_cent_ellipse = fitdist(maj_axis_cent_ellipse, 'Normal');
pd_min_axis_cent_ellipse = fitdist(min_axis_cent_ellipse, 'Normal');

% Generate kernel inv gaussian distributions for central elongated vess
pd_maj_axis_cent_elongat = fitdist(maj_axis_cent_elongat, 'InverseGaussian');
pd_min_axis_cent_elongat = fitdist(min_axis_cent_elongat, 'Normal');

% Generate kernel inv gaussian distributions for boundary vess
pd_maj_axis_boun = fitdist(maj_axis_boun, 'InverseGaussian');
pd_min_axis_boun = fitdist(min_axis_boun, 'InverseGaussian');

% Generate kernel inv gaussian distributions for reciprocal vess distance
pd_distances = fitdist(distances, 'InverseGaussian');

%% Plot  Distributions

% Plot major axis for central elliptic vessels
x_values = linspace(min(maj_axis_cent_ellipse), max(maj_axis_cent_ellipse), 1000);
pdf_values = pdf(pd_maj_axis_cent_ellipse, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Major Axis Central Elliptical Vessels Probability');

% Plot minor axis for central elliptic vessels
x_values = linspace(min(min_axis_cent_ellipse), max(min_axis_cent_ellipse), 1000);
pdf_values = pdf(pd_min_axis_cent_ellipse, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Minor Axis Central Elliptical Vessels Probability');

% Plot major axis for central elongated vessels
x_values = linspace(min(maj_axis_cent_elongat), max(maj_axis_cent_elongat), 1000);
pdf_values = pdf(pd_maj_axis_cent_elongat, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Major Axis Central Elongated Vessels Probability');

% Plot minor axis for central elongated vessels
x_values = linspace(min(min_axis_cent_elongat), max(min_axis_cent_elongat), 1000);
pdf_values = pdf(pd_min_axis_cent_elongat, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Minor Axis Central Elongated Vessels Probability');

% Plot major axis for boundary vessels
x_values = linspace(min(maj_axis_boun), max(maj_axis_boun), 1000);
pdf_values = pdf(pd_maj_axis_boun, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Major Axis Boundary Vessels Probability');

% Plot minor axis for boundary vessels
x_values = linspace(min(min_axis_boun), max(min_axis_boun), 1000);
pdf_values = pdf(pd_min_axis_boun, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Axis Size [um]');
ylabel('Probability [ ]');
title('Minor Axis Boundary Vessels Probability');

% Plot minor axis for boundary vessels
x_values = linspace(min(distances), max(distances), 1000);
pdf_values = pdf(pd_distances, x_values);
figure
plot(x_values, pdf_values, 'LineWidth', 1.5);
xlabel('Distance [um]');
ylabel('Probability [ ]');
title('PC3 Vessels Distance Distribution');

% Save distributions
if flag_save
    
    % Save in .mat file
    save('VesselsProperties_PC3.mat', 'pd_maj_axis_cent_ellipse', ...
         'pd_min_axis_cent_ellipse', 'pd_maj_axis_cent_elongat', ...
         'pd_min_axis_cent_elongat', 'pd_maj_axis_boun', ...
         'pd_min_axis_boun', 'pd_distances', 'pc3_vess_boundary_coverage',...
         'pc3_ellipse_bv', 'pc3_elongat_bv')    
end
