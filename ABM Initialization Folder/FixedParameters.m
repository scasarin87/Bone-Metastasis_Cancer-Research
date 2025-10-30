%% Fixed Parameters %%

%  This function defines some fixed parameters which I won't need to modify
%  during the ABM simulations

% 1) Standard Time and Site Constants

% hour            = 0;  % Initialization of the follow-up hours
days_hours      = 24; % Set the day hour
T_cell_division = 24; % Time after which PCa cell will undergo mit/apo event
if strcmp(cell_line, 'pc3')
    T_ocs_resorption = 36; % Time after which Ocs resorption is triggered
else
    T_ocs_resorption = 72;
end


if strcmp(cell_line, 'c42b')
    follow_up       = 15 * days_hours;
elseif strcmp(cell_line, 'pc3')
    follow_up       = 21 * days_hours;
elseif strcmp(cell_line, 'renca')
    follow_up       = 28 * days_hours;    
end    

site_dim   = 500/24;        % ABM scaling factor (1 pix = 20.833 um). 
                            % Represents the diameter of a single PCa cell                        

% 2) Rad223 Dynamics Parameters

% 2a - Rad Coefficients
Rad = struct(...
        'mitosis'   , (alpha.p_mit_max + alpha.p_mit_min)/3, ...  
        'apoptosis' , (alpha.p_apo_max + alpha.p_apo_min)/3, ... 
        'quiescence', 1 - ((alpha.p_mit_max + alpha.p_mit_min)/2) + ((alpha.p_apo_max + alpha.p_apo_min)/2), ... % 1 - rad.mitosis - rad.apoptosis 
        'activity'  ,   0);    % Turn Rad-223 activity ON by switching the
                               % activity marker

% 2b - Rad Temporal Parameter: Rad decays exponentially from 1(max) toward 0,
%      by halving itself every 11 days.
rad_max = 1;  % Max activity is intended as unitary
tau     = 16; % decay coefficients fitted from experimental data
half_life_time_rad223 = 11 * days_hours; % rad223 half life time

% 2c - Distance-dependent: The intrinsic impact of Rad on mitosis and
%      apoptosis densities of probability (intrinsic = independent from the 
%      very lesion) is scaled on the current lesion which sees a 50% chanches
%      of mitosis and a 10% of apoptosis (in the innermost portion)
[B_mit, B_apo] = get_radium_coefficients(Rad);

% 3) Cabonzantinib Parameters

flag_cabo = 0; % if flag_cabo = 1 -> cabozantinib therapy is currently 
               % active
               
flag_min_distance = 0; % check ABM Tumor Dynamics -> events_probability_control/Cabo   

vessels_time = zeros(1, 3); % 1st Col -> vessel center row
                            % 2nd Col -> vessel center column
                            % 3rd Col -> time passed after cabo targeting 

vessel_retard = 4.5 * days_hours; % After this time, vessels targeted by cabozantinib
                                  % will stop giving contribution to
                                  % mitosis and apoptosis probabilities
                                  
%res_area_3_days = 1; % defined fixed parameters, it will change when 
                     % zoledronic acid therapy is initiated                                  


