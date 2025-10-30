%% Compute Distance %%

%  This function is used to compute distance in the hexagonal grid between
%  two points -> Look Utils

%  Input  -> X, Y     : hexagonal grid 
%            x1, y1   : row and column coordinates of point A
%            x2, y2   : row and column coordinates of point B
%
%  Output -> distance : distance between point A and B 

function [distance] = compute_distance(X, Y, x1, y1, x2, y2)

    distance = sqrt((X(x1, y1) - X(x2, y2)) ^ 2 + (Y(x1, y1) - Y(x2, y2)) ^ 2); 

end