%% Grid Features Definition
%  This function is used to set important features of the grid, which will
%  be used in the following steps of the ABM dynamics

% Input  -> X, Y   : hexagonal grid 
%
% Output -> ax, bx : min & max grid values on x-axis
%           ay, by : min & max grid values on y-axis
%           Cx, Cy : center of the grid
%           directionx, directiony : defines path orientation for mitosis and apoptosis

function [ax, ay, bx, by, Cx, Cy, directionx, directiony] = grid_FeaturesDefinition(X, Y)

    % Dimensions of the box
    ax=min(min(X)); bx=max(max(X));
    ay=min(min(Y)); by=max(max(Y)); 

    % Center of the grid
    Cx = (ax+bx)/2; Cy = (ay+by)/2; 

    % Grid navigation vectors
    directionx      = zeros(6,2);
    directiony      = zeros(1,6);

    % if k is odd: 
    directionx(1,1) =  0; directiony(1) =  1; % north-west
    directionx(2,1) =  1; directiony(2) =  1; % north-east
    directionx(3,1) =  1; directiony(3) =  0; % east
    directionx(4,1) =  1; directiony(4) = -1; % south-east
    directionx(5,1) =  0; directiony(5) = -1; % south-west
    directionx(6,1) = -1; directiony(6) =  0; % west

    % if k is even:
    directionx(1,2) = -1; directiony(1) =  1; % north-west
    directionx(2,2) =  0; directiony(2) =  1; % north-east
    directionx(3,2) =  1; directiony(3) =  0; % east
    directionx(4,2) =  0; directiony(4) = -1; % south-east
    directionx(5,2) = -1; directiony(5) = -1; % south-west
    directionx(6,2) = -1; directiony(6) =  0; % west
end
