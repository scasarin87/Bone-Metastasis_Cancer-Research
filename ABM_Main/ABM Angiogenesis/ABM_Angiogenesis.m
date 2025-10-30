%% ABM Angiogenesis %%

%  This function performs the angiogenesis: under control condition we set
%  3 different angiogenesis steps:
%  A - Every 6  hours
%  B - Every 20 hours
%  C - Every 30 hours


if strcmp(cell_line, 'c42b') && hour == 1
    load VesselsProperties_C42B
    Prob_ott = zeros(8,1);
elseif strcmp(cell_line, 'pc3') && hour == 1
    load 'ABM Main'/VesselData/VesselsProperties_PC3.mat
elseif strcmp(cell_line, 'renca') && hour == 1
    load 'ABM Main'/VesselData/VesselsProperties_RENCA.mat
end

if strcmp(cell_line, 'pc3') 
    if mod(hour, 3) == 0

        pc3_boundary_vessel_dynamics;

    end
end

if strcmp(cell_line, 'renca') 
    if mod(hour, 3) == 0

        renca_vessel_dynamics;

    end
end
%     
% C42B Angiogenesis every 6 hours
if strcmp(cell_line, 'c42b')
    % Six Hours Angiogenesis Process
    if (mod(hour, 6) == 0 && flag_cabo == 0) 
        SixHour_Angiogenesis_C42B;     
    end  
end
% % PC3 Angiogenesis every 6 hours
% else    
% %     % Vessel push towards the boundaries
% %     if (mod(hour, 6) == 0) 
% %         ReDistribute_Vessels;
% %     end
% %     % Every 6 hours add a vessel if vess/tum area ratio for pc3 is below threshold
% %     if ((mod(hour, 6) == 0 && flag_cabo == 0)) % || mod(hour, 25) == 0) && flag_cabo == 0)  
% %         Angio_PC3;
% %     end     
% %     %if (mod(hour, 24) == 0 && flag_cabo == 0)  
% %     %    Angio_PC3_24Hours;
% %     %end  
% end

% C42B Angiogenesis every 20 hours
if strcmp(cell_line, 'c42b')
    % Twenty Hours Angiogenesis Process
    if (mod(hour, 20) == 0 && flag_cabo == 0)     
        TwentyHour_Angiogenesis_C42B;    
    end
end
% % PC3 Angiogenesis every 20 hours
% else    
%     % Twenty Hour Angiogenesis Process
%     if (mod(hour, 20) == 0 && flag_cabo == 0) 
%         % TwentyHour_Angiogenesis_PC3;     
%     end  
% end
%     
% 
% C42B Angiogenesis every 30 hours
if strcmp(cell_line, 'c42b')
    % Thirty Hour Angiogenesis Process
    if (mod(hour, 30) == 0 && flag_cabo == 0)     
        ThirtyHour_Angiogenesis_C42B;    
    end
end
% % PC3 Angiogenesis every 30 hours
% else    
%     % Thirty Hour Angiogenesis Process
%     if (mod(hour, 30) == 0 && flag_cabo == 0) 
%         % ThirtyHour_Angiogenesis_PC3;     
%     end  
% end
% 
% % if (mod(hour, 30) == 0 && flag_cabo == 1)
% %     CaboAngio;
% % end
% 
% % Update Vessels Number After Angiogenesis
% vessels_number = size(center_vessels, 1);
% 
