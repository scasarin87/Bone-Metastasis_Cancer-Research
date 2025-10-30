%% Control In-Vivo Data %%

%  This function contains the std and the average value of mice tibiae 
%  tumor growth under control condition performed in laboratory. This data
%  are our gold standard to replicate tumor growth with our ABM.

%% C42B Cell Line Data
% Experimental Time Points
time_points = [1, 5 * 24, 8 * 24, 12 * 24, 15 * 24];
% Average mice tumor growth 
MEAN = [1 7.128197054 16.08757084 29.31455403 45.79737765];
% Standard Deviation at the same time points
STD = [0 3.719900472 5.890299748 10.73679398 16.58286703];

% Experimental Data Interpolated Data
control_data_interpolated = interp1(time_points, MEAN, 1:1:360, 'pchip'); 

%% PC3 Cell Line Data

% Define PC3 Control Data and STD
CONTROL_MEDIA = [1 8.0622 20.768 69.061 132.573]; 
CONTROL_STD = [0 4.141 14.23 41.47 62.536];
% Time points
time_points_pc3=[0, 4, 8, 17, 21] * 24;  
% Get interpolated curves
control_curve_pc3 = interp1(time_points_pc3, CONTROL_MEDIA, 1:1:21*24, 'pchip'); 
%control_curve_pc3_old = interp1(time_points_pc3_old, CONTROL_MEDIA_OLD, 1:1:21*24, 'pchip'); 
% I need the curve for the first 15 days
% control_curve_pc3 = control_curve_pc3(1:360);
% CONTROL_MEDIA = CONTROL_MEDIA(1:4);
% CONTROL_STD = CONTROL_STD(1:4);
% time_points_pc3 = time_points_pc3(1:4);

%% RENCA CELL LINE DATA

CONTROL_MEDIA_RENCA = [1, 17.612, 32.674, 141.202];
CONTROL_STD_RENCA = [0, 10.125, 22.602, 81.016];
time_points_renca = [0, 7, 10, 15] * 24;
control_curve_renca = interp1(time_points_renca, CONTROL_MEDIA_RENCA, 1:1:15*24, 'pchip');

%% RENCA LESION DEVELOPMENT

renca_lesion_average = [0.070, 0.128, 0.437, 0.564];
renca_lesion_std = [0.005, 0.030, 0.141, 0.106];
time_point_renca_lesion = [14, 21, 28, 35] * 7;

%% Cabo Data from Varkaris Paper
% Cabozantinib Effect Experimental Data (from Varkaris Paper)
cabo_data;


%% Osteoclast in vivo distribution
global top_slice_ocs mid_slice_ocs low_slice_ocs calvaria_ocs_min calvaria_ocs_max vertebra_ocs_min vertebra_ocs_max;

top_slice_ocs = [31.492 47.346 27.122 54.002];
top_slice_ocs = mean(top_slice_ocs);
mid_slice_ocs = [14.838 31.515 21.855 8.448];
mid_slice_ocs = mean(mid_slice_ocs);
low_slice_ocs = [5.876 5.876 0.314 0.668];
low_slice_ocs = mean(low_slice_ocs);

% First (min) and third (max) quartile values of ocs in in vivo distribution
calvaria_ocs_min = 9/2; 
calvaria_ocs_max = 14/2;
vertebra_ocs_min = 13/2;
vertebra_ocs_max = 20/2;

%% Osteoblast in vivo distribution

global obs_coverage obs_rad_decrease_curve;

obs_coverage = [95.31 96.389 92 98.317 98.999 97.59 ...
                99.106 97.989 99.387 77.276 95.49 99.66];
obs_coverage = mean(obs_coverage); 

% PC3 Vessels ratio against tumor
global ratio_pc3_vess_tum ratio_pc3_vess_tum_cabo;
ratio_pc3_vess_tum = 6/100;
ratio_pc3_vess_tum_cabo = 3/100;



