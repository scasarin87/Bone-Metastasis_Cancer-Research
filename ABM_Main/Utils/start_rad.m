%% Start Rad %%

%  This function is used to start the administration of Rad223 based
%  therapy by activating the Rad.activity flag ( = 1)

%  Input  -> Old_Rad          : current struct containg Rad properties
%            current_hour     : current simulation hour
%            start_therapy_rad: hour at which user want Rad223 therapy to
%                               start
%
%  Output -> New_Rad          : updated struct containg Rad properties

function [New_Rad] = start_rad(Old_Rad, current_hour, start_therapy_rad)

    if current_hour == start_therapy_rad
        Old_Rad.activity = 1; % Rad223 therapy is activated
    end
    
    New_Rad = Old_Rad;

end