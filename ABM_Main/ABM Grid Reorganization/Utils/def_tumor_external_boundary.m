%% Define Tumor External Boundary %%

%  This function is used to define the ABM site at the interface between
%  tumor and bone marrow. Basically we obtain the bone marrow site adjacent
%  to the tumor cells

%  Input  -> marrow  : mask complementary to tumor

%  Output -> marrow_edge: BM sites adjacent to tumor

function [marrow_edge] = def_tumor_external_boundary(marrow, bone, site)

    % Thanks to https://www.youtube.com/watch?v=5yaj7QWurgM
    method = strel('disk', 1);
    marrow_edge = imerode(marrow, method);
    marrow_edge = marrow - marrow_edge;
    
    % Find Marrow_Edge Mask Coordinates of ones
    [row_marrow, col_marrow] = find(marrow_edge == 1);
    
    % Check for each mask 1 if belongs to BM, otherwise put it to 0
    for cell = 1 : size(row_marrow, 1)
        if bone(row_marrow(cell), col_marrow(cell)) == site.cortical_bone ...
                || bone(row_marrow(cell), col_marrow(cell)) == site.osteoblast ...
                || bone(row_marrow(cell), col_marrow(cell)) == site.osteoclast ...
                || bone(row_marrow(cell), col_marrow(cell)) == site.vessel ...
                || bone(row_marrow(cell), col_marrow(cell)) == site.vessel_cabo
            marrow_edge(row_marrow(cell), col_marrow(cell)) = 0;
        end
    end
    

    % Add the following line to get only the outer boundary
    % marrow_edge = imfill(marrow_edge, 'holes');
    % marrow_edge = bwperim(marrow_edge);
    
end