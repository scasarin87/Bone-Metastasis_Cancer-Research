%% Define ABM Driving Parameters %%

%  This function creates the struct containing the 6 model driving parameters
%  You can import the values by writing alpha.nameparameter (e.g.
%  alpha.p_mit_max)

global alpha

% C42B
if strcmp(cell_line, 'c42b')

    alpha = struct(...
            'p_mit_max'     , 8.67779667e-01, ...  % Maximum probability of mitosis
            'p_mit_min'     , 2.85883871e-01, ...  % Minimum probability of mitosis
            'slope_prob'    , 1.33384558e-01, ...  % Tumor Rate Growth
            'p_apo_max'     , 4.20590184e-01, ...  % Maximum probability of apoptosis
            'p_apo_min'     , 2.27033839e-01, ...  % Minimum probability of mitosis
            'vess_influence', 1.41262732e+02);     % Variable in µm describing the maximum 
                                                   % cell-vessel distance according to which
                                                   % a vessel is still able to provide a Mitosis 
                                                   % probability increase.
%PC3    
elseif strcmp(cell_line, 'pc3')
    
    alpha = struct(...
            'p_mit_max'     , 0.55667,  ... % Maximum probability of mitosis
            'p_mit_min'     , 0.43449,  ... % Minimum probability of mitosis
            'slope_prob'    , 0.062574, ... % Tumor Rate Growth
            'p_apo_max'     , 0.39985,  ... % Maximum probability of apoptosis
            'p_apo_min'     , 0.087551, ... % Minimum probability of mitosis
            'vess_influence', 157.4313);    % Variable in µm describing the maximum 
                                            % cell-vessel distance according to which
                                            % a vessel is still able to provide a Mitosis 
                                            % probability increase.           
                                                   
elseif strcmp(cell_line, 'renca')
    
    alpha = struct(...
            'cabo_max',        5.19, ...
            'p_mit_max'     , 6.56599849e-01, ...%MAPE6.55844045e-01,  ... % Maximum probability of mitosis
            'p_mit_min'     , 3.34848582e-01, ...%3.35406850e-01,  ... % Minimum probability of mitosis
            'slope_prob'    , 1.79787160e-01, ...%1.74113444e-01, ... % Tumor Rate Growth
            'p_apo_max'     , 1.81329122e-01, ...%1.38261550e-01,  ... % Maximum probability of apoptosis
            'p_apo_min'     , 1.23429699e-01, ...%1.25070755e-01, ... % Minimum probability of mitosis
            'vess_influence', 1.49908143e+02); ...%1.53999838e+02);    % Variable in µm describing the maximum 
                                            % cell-vessel distance according to which
                                            % a vessel is still able to provide a Mitosis 
                                            % probability increase.      
    
end