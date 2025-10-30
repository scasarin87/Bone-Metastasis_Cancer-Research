%% Tumor Dimension %% 

% the code defines the initial tumor size according to experimental plan

%  Input  -> Dim_Flag  : tumor dimension chosen by the user
%                 N_Cells : the number of cells a microtumor has initially
%                 (may go from 0 to 4)
%
%  Output -> a_tumor   : tumor semi-axis on the X-axis
%            b_tumor   : tumor semi-axis on the Y-axis


function [a_tumor,b_tumor] = tumor_dimension(Dim_Flag, N_Cells)
    
    % Initialization of tumor semi-axis dimensions
    a_tumor_init = [2, 2, 5, 7, 9, 12, 15, 18, 20, 23]; 
    b_tumor_init = [1, 3, 4, 6, 7, 10, 12, 15, 17, 20];
        
    % Select the corresponding dimensions (according to Dim_Flag value)
    % if Dim_Flag > 0 then normal tumor dynamics
    if Dim_Flag > 0
    for i_size = 1:10
        if i_size == Dim_Flag
           a_tumor = a_tumor_init(i_size);
           b_tumor = b_tumor_init(i_size);
        end
    end  
    % otherwise miroctumor
    else
        switch N_Cells
            case 1
                a_tumor = 0;
                b_tumor = 1;
            case 2
                if rand > 0
                    [a_tumor, b_tumor] = deal(0, 2);
                else
                    [a_tumor, b_tumor] = deal(2, 0);
                end
            case 3
                if rand > 0
                    [a_tumor, b_tumor] = deal(1, 2);
                else
                    [a_tumor, b_tumor] = deal(2, 1);
                end
        end
    end

end












