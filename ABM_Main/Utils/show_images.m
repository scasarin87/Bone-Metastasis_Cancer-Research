clear 
close all
clc

% Get screen size to customize visualization
screenSize = get(0, 'ScreenSize');
width = screenSize(3);
heigth = screenSize(4);

% Define filename to load
filename = 'pc3_tibia_cs_2_control';

% Load filename
% load(filename)
%load('BONE_8.mat')

% Split the name components
parts = strsplit(filename, '_');

% Get cell line info
cell_line = parts{1};
% Get bone
bone = parts{2};    
% Get section type
section = parts{3};

% Define the number of temporal screenshots
if strcmp(cell_line, 'renca')
    time = [1, 7, 3*24, 6*24, 9*24, 11*24, 13*24, 15*24];
end

if strcmp(cell_line, 'pc3')
    time = [1, 7, 3*24, 6*24, 9*24, 11*24, 13*24, 15*24, 17*24, 19*24, 21*24];
end

% Define figure size
if strcmp(bone, 'femur')
    if strcmp(section, 'ls')
        fig_size = [0, 0, 220, 850]; % [left, bottom, width, height]
    end
    if strcmp(section, 'cs')
        fig_size = [0, 100, 300, 280]; % [left, bottom, width, height]
    end
end

if strcmp(bone, 'tibia')
    if strcmp(section, 'ls')
        fig_size = [0, 0, 185, 850]; % [left, bottom, width, height]
    end
    if strcmp(section, 'cs')
        fig_size = [0, 100, 400, 300]; % [left, bottom, width, height]
    end
end

if strcmp(bone, 'calvaria')
    fig_size = [0, 0, 185, 850]; % [left, bottom, width, height]
end

if strcmp(bone, 'vertebra')
    fig_size = [0, 100, 400, 300]; % [left, bottom, width, height]
end

% initialize bone
bone = BONE(:, :, 1);
init_tum = find(bone == 2);
init_ves = find(bone == 3);
bone(init_tum) = 0;
bone(init_ves) = 0;
BONE(:, :, 1) = bone;


for i = length(time) : -1 : 1
    figure(i)
    h = gcf; % Get the handle to the current figure
    set(h, 'Position', fig_size); 
    [ABM_matrix] = display_abm(size(BONE, 1), size(BONE, 2), BONE(:, :, time(i)), 1); 
    if fig_size(1) + 2 * fig_size(3) < width
        fig_size(1) = fig_size(1) + fig_size(3);
    else
        fig_size(1) = 0;
    end
end
