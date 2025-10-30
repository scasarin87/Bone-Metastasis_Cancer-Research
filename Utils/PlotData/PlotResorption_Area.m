%% Plot Resorption Area

if flag_resorption
    
    % Compute the resorption area at each simulation area
    cortical_bone_sites = abs(cortical_bone_sites - cortical_bone_sites(1));
    
    % Transform the area to mm^2
    resorption_area = cortical_bone_sites * (site_dim / 1000) * (site_dim / 1000);
    
    % Get the value at two weeks
    %resorption_area_2weeks(epoch, 1) = resorption_area(14 * days_hours);
    
    all_resorption_areas(epoch, :) = resorption_area;
    
end