%% Define Split Tumor Quadrant %%

%  This function is used to split tumor mask into 4 quadrants. Depending on
%  the relationship between tumor centroid and the current PCa cell which
%  is undergoing mitosis or apoptosis, we will consider only the correct
%  portion of bone marrow sites adjacent to the tumor.

%  Input  -> marrow_edge         : BM sites adjacent to tumor
%            row/col_centroid    : tumor row and column centroid
%            row/col_agent       : current PCa cell row and col
%
%  Output -> marrow_edge_quadrant: BM sites adjacent to tumor in the chosen
%                                  quadrant


function [marrow_edge_quadrant] = def_split_tumor_quadrant(marrow_edge, row_centroid, ...
                                                             col_centroid, row_agent, col_agent)
    
    % UPPER RIGHT QUADRANT
    if col_agent >= col_centroid && row_agent < row_centroid
        marrow_edge(row_centroid : end, :)   = 2;
        marrow_edge(:, 1 : col_centroid - 1) = 2;
    end
    
    % UPPER LEFT QUADRANT
    if col_agent < col_centroid && row_agent < row_centroid
        marrow_edge(row_centroid : end, :) = 2;
        marrow_edge(:, col_centroid : end) = 2;
    end
        
    % LOWER LEFT QUADRANT
    if col_agent <= col_centroid && row_agent >= row_centroid
        marrow_edge(1 : row_centroid - 1, :)   = 0;
        marrow_edge(:, col_centroid + 1 : end) = 0;
    end

    % LOWER RIGHT QUADRANT
    if col_agent > col_centroid && row_agent >= row_centroid
        marrow_edge(1 : row_centroid - 1, :) = 2;
        marrow_edge(:, 1 : col_centroid)     = 2;
    end
    
    marrow_edge_quadrant = marrow_edge(:, :);                                               
                                                         
end                                                         