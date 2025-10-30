%% Monte Carlo Algorithm %%

% This function performs the montecarlo algorithm: given the mitosis and
% apoptosis probabilities of the current active tumor cell (at coordinates
% row_agent, col_agent), it defines whether it undergoes mitosis, 
% apoptosis or stays quiescient by comparing the computed probabilities 
% with test probabilities. It guarantees stochastic character to the virtual
% experiment.
%
%  Input  -> p_mit, p_apo    : current tumor cell computed probabilities
%            row/col_agent   : current tumor site coordinates (row/column)
%            change_cell_old : is the current change_cell matrix (see
%                              FixedMatrix function)
%
% Output  -> change_cell     : updated matrix depending on the algorithm
%                              results


function [change_cell, mitotic_cells, apoptotic_cells] = monte_carlo_algorithm(p_mit, p_apo, row_agent, ...
                                               col_agent, change_cell_old, ...
                                               mitotic_cells, apoptotic_cells, site)
     
    % Define Indexes
    site_in_mitosis = 1; site_in_apoptosis = -1;
    
    % Define the test mitosis/apoptosis number for the M.Carlo Method
    p_mit_test  = rand(1); p_apo_test = rand(1);   
    
    % CASE 1: Apoptosis with no Conflict
    if p_apo_test < p_apo && p_mit_test > p_mit
       change_cell_old(row_agent, col_agent) = site_in_apoptosis;
       apoptotic_cells(row_agent, col_agent) = site.apoptotic_cell;
    end  
                
    % CASE 2: Mitosis with no Conflict
    if p_apo_test > p_apo && p_mit_test < p_mit
       change_cell_old(row_agent, col_agent) = site_in_mitosis;
       mitotic_cells(row_agent, col_agent) = site.mitotic_cell;
    end
    
    % CASE 3: When Both Events can Potentially Happen
    if (p_apo_test < p_apo && p_mit_test < p_mit) 
                    
       % The event that will occur is the one with most prevalence 
       if (p_apo - p_apo_test) / p_apo > (p_mit - p_mit_test) / p_mit
           change_cell_old(row_agent, col_agent) = site_in_apoptosis;
           apoptotic_cells(row_agent, col_agent) = site.apoptotic_cell;
       else
           change_cell_old(row_agent, col_agent) =  site_in_mitosis;
           mitotic_cells(row_agent, col_agent) = site.mitotic_cell;
       end 
                    
    end    
    
    change_cell = change_cell_old;
                                              
end
