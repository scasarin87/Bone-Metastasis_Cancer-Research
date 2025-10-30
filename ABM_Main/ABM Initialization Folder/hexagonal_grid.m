%% Hexagonal Grid %%

%  This function creates the hexagonal grid upon which the ABM will be
%  built. It also generates important grid features.

%  Input  -> rows, columns   : size of the ABM grid
%
%  Output -> X, Y   : hexagonal grid 
%            ax, bx : min & max grid values on x-axis
%            ay, by : min & max grid values on y-axis
%            Cx, Cy : center of the grid (central_col, central_row)
%            directionx, directiony : defines path orientation for mitosis and apoptosis

function [X, Y, ax, ay, bx, by, ...
            Cx, Cy, directionx, directiony] = hexagonal_grid(rows, columns)
        
    X = zeros(rows, columns); % X refers to grid columns 
    Y = zeros(rows, columns); % Y refers to grid rows 
    
    for jj = 1:rows
        for kk = 1:columns

            if mod(kk, 2) == 0
                Y(jj, kk) = 2 * jj; 
                X(jj, kk) = 2 * kk;
            else
                Y(jj, kk) = 2 * jj + 1; 
                X(jj, kk) = 2 * kk;
            end

        end
    end
    
    clear jj kk
    
    X = 0.5 .* X; Y = 0.5 .* Y;
    
    % Dimensions of the box
    ax = min(min(X)); bx = max(max(X));
    ay = min(min(Y)); by = max(max(Y)); 

    % Center of the grid
    Cx = (ax + bx)/2; Cy = (ay + by)/2; 

    % Grid navigation vectors
    directionx      = zeros(6, 2);
    directiony      = zeros(1, 6);

    % if k is odd: 
    directionx(1, 1) =  0; directiony(1) =  1; % north-west
    directionx(2, 1) =  1; directiony(2) =  1; % north-east
    directionx(3, 1) =  1; directiony(3) =  0; % east
    directionx(4, 1) =  1; directiony(4) = -1; % south-east
    directionx(5, 1) =  0; directiony(5) = -1; % south-west
    directionx(6, 1) = -1; directiony(6) =  0; % west

    % if k is even:
    directionx(1, 2) = -1; directiony(1) =  1; % north-west
    directionx(2, 2) =  0; directiony(2) =  1; % north-east
    directionx(3, 2) =  1; directiony(3) =  0; % east
    directionx(4, 2) =  0; directiony(4) = -1; % south-east
    directionx(5, 2) = -1; directiony(5) = -1; % south-west
    directionx(6, 2) = -1; directiony(6) =  0; % west
        
end