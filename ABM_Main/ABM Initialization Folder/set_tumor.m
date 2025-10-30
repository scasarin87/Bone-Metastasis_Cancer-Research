%% Set Tumor %%

%  This function is used to create the tumor mask and initialize their
%  internal clock with a random hour between 1 and 24h. Remember that each
%  PCa cell will undergo either mitosis or apoptosis at the end of the 24h.

%  Input  -> rows, columns : Y and X size of the grid
%            X, Y          : hexagonal grid
%            CxT, CyT      : tumor center coordinates
%            T_cell_divis  : After * hours a PCa will undergo mit/apo event
%            tumor_dim     : tumor dimension chosen by the users
%
%  Output -> tumor         : tumor mask -> = 1 if site belongs to tumor
%            internal_clock: matrix that stores each PCa cell "age" [h] 

function [tumor, vess_prob, internal_clock, a_tumor, b_tumor] = set_tumor(rows, columns, X, Y, T_cell_division)
    
    global CxT CyT tumor_dim n_cells cell_line alpha
    
    internal_clock  = zeros(rows, columns); % keep track of the cell time across the division cycle
    tumor           = zeros(rows, columns); % tumor = 1 if site belongs to tumor, 0 otherwise. 
    vess_prob       = zeros(rows, columns);
        
    % Intial tumor dimension is chosen according to experimental plan
    [a_tumor , b_tumor] = tumor_dimension(tumor_dim, n_cells); % a_tumor/b_tumor are the ellipse semiaxis
    
    for jj = 1:rows
        for kk = 1:columns
        
            % Get the distance from the center 
            dist = compute_distance(X, Y, jj, kk, CxT, CyT);        
            % Get the radius
            radius = ellipse_radius(a_tumor, b_tumor, abs(X(jj, kk) - X(CxT, CyT)), abs(Y(jj, kk) - Y(CxT, CyT)));
                    
            % Fill the matrix
            if dist <= radius
                tumor(jj, kk) = 1;
                internal_clock(jj, kk) = round(rand(1) * T_cell_division); % also assign a random clock between 1 and 24
                if strcmp(cell_line, 'pc3')
                    % Keep decreasing radius value to assign to each circular crown a different probability
                    for rad_perc = 99 : -2 : 1
                        curr_rad = ellipse_radius((rad_perc / 100) * a_tumor, (rad_perc / 100) * b_tumor, abs(X(jj, kk) - X(CxT, CyT)), abs(Y(jj, kk) - Y(CxT, CyT)));
                        % Fill the matrix
                        if dist < curr_rad - 0.15 * curr_rad
                            curr_prob = custom_pdf(curr_rad, radius, 'exp', 'false');
                            vess_prob(jj, kk) = curr_prob;
                        end    
                    end
                end 
            end
        end         
    end         
     
    clear jj kk