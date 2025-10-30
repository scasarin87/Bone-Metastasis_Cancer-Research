%% Start Cabo %%

%  This function is used to start the administration of Cabozantinib
%  therapy by activating the flag cabo ( = 1)

%  Input  -> old_flag_cabo     : current flag_cabo value
%            old_center_vessels: old Nvess x 4 matrix containing vess info   
%            current_hour      : current simulation hour
%            start_therapy_cabo: hour at which user want cabo therapy to
%                                start
%
%  Output -> new_flag_cabo     : updated flag_cabo value
%            old_center_vessels: updated center_vessels variable
%            n_vess            : vessels number when cabo starts
%            vess_eliminated   : set to 0 the vessels already elim by cabo


function [new_flag_cabo, new_center_vessels, n_vess, vess_eliminated] = ...
                            start_cabo(old_flag_cabo, old_center_vessels)
    
    % Cabozantinib Effect Started
    old_flag_cabo = 1; 
    % Add Column for Vessels Targeted by Cabo
    old_center_vessels(:, 5) = 0; 
    % Initialize Variables to be used in the vessel elimination section
    n_vess = size(old_center_vessels, 1); % Current vessel number
    % Initialize the number of eliminated vessels
    vess_eliminated = 0; 
    
    % Update Variables
    new_flag_cabo = old_flag_cabo;
    new_center_vessels = old_center_vessels;    
                                                
end