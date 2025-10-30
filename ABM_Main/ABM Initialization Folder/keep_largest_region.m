%% Keep Largest Region %%

%  This function outputs a new binary mask that considers only the largest 
%  area among all the areas within the original one. It is useful when i do 
%  not want to take into account the smallest regions but just the largest.

% Input  -> mask: binary mask 
%
% Output -> mask: same binary mask with only the largest region considered

function [mask] = keep_largest_region(mask)
    
    % Find connected components in the binary mask.
    mask = bwlabel(mask);
    % Use regionprops to extract properties of connected components.
    stats = regionprops(mask, 'Area');
    areas = [stats.Area];
    % Find the index of the biggest connected component.
    [~, largest_region_index] = max(areas);
    % Create a new binary mask, setting all pixels except the biggest component to zero.
    mask = mask == largest_region_index;
end