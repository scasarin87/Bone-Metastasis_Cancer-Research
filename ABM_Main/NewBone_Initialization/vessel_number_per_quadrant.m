%% Vessel Number Per Quadrant %%

%  This function is used to define the tumor quadrant in which the fewer
%  number of vessels were generated. The idea is to provide a uniform
%  distribution of the vessels in the tumor area for the PC3 cells.

%  Input  -> tumor  : binary mask defining tumor location in the 2D space
%            vessels: binary mask defining vessel locations in the 2D space
%
%  Output -> out_area: tumor area with fewer vessels

function [out_area] = vessel_number_per_quadrant(tumor, vessels)
    
    % Use regionprops to calculate properties of connected regions
    props = regionprops(tumor, 'Area', 'Centroid');

    % Extract the areas of all regions (tumor can split into smaller parts)
    areas = [props.Area];

    % Find the index of the region with the largest area
    [maxArea, maxIndex] = max(areas);

    % Extract the centroid of the largest region
    centroids = props(maxIndex).Centroid;
    
    % Extract row and col centroid 
    c_row = round(centroids(2));
    c_col = round(centroids(1));
    
    % Define the upper right quadrant
    upper_right = zeros(size(tumor));
    upper_right(tumor == 1) = 1; 
    upper_right(c_row + 1 : size(tumor, 1), :) = 0;
    upper_right(:, 1 : c_col - 1) = 0;
    % Define the lower right quadrant
    lower_right = zeros(size(tumor));
    lower_right(tumor == 1) = 1; 
    lower_right(1 : c_row - 1, :) = 0;
    lower_right(:, 1 : c_col - 1) = 0;
    % Define the lower left quadrant
    lower_left = zeros(size(tumor));
    lower_left(tumor == 1) = 1; 
    lower_left(1 : c_row - 1, :) = 0;
    lower_left(:, c_col + 1 : size(tumor, 2)) = 0;
    % Define the upper left quadrant
    upper_left = zeros(size(tumor));
    upper_left(tumor == 1) = 1; 
    upper_left(c_row + 1 : size(tumor, 1), :) = 0;
    upper_left(:, c_col + 1 : size(tumor, 2)) = 0;
    
    % Assess the number of vessels for each quadrant 
    vessels_upper_right = size(find(vessels == 1 & upper_right), 1);
    vessels_lower_right = size(find(vessels == 1 & lower_right), 1);
    vessels_lower_left = size(find(vessels == 1 & lower_left), 1);
    vessels_upper_left = size(find(vessels == 1 & upper_left), 1);
    % Store values in a vector
    vessel_numbers = [vessels_upper_right, vessels_lower_right, vessels_lower_left, vessels_upper_left];
    % Define the max value
    min_vessels = min(vessel_numbers);
    % Find the indices of the maximum values in the vector
    min_indices = find(vessel_numbers == min_vessels);
    
    % If multiple quadrants have the same number of vessels
    if size(min_indices, 2) > 1
        % Randomly select one index from maxIndices
        out_index = min_indices(randi(length(min_indices)));
    else 
        out_index = min_indices;
    end
    
    if out_index == 1
        out_area = upper_right;
    elseif out_index == 2
        out_area = lower_right;
    elseif out_index == 3
        out_area = lower_left;
    elseif out_index == 4
        out_area = upper_left;
    end 

end