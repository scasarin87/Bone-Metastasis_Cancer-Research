%% Site Dynamic %%

%  The function performs the cellular activity described by change_cell matrix
%  (i.e. mitosis and apoptosis). First, we look for a random site belonging
%  to the bone marrow, the outward remodeling of the tumor will
%  be directed towards that very site, while the inward remodeling will start 
%  from that very site. The path draw is valid both for mitosis and apoptosis. 
%  First site is always the closest site to the active one (row_agent, 
%  col_agent), and the last site is the first bone site that the path finds
%  on its way.


% Mitosis Event
if change_cell(row_agent, col_agent) == 1 
    
    % Define Mitotic Site
    [row_mitosis, col_mitosis] = define_mitotic_site(bone, site, row_agent, col_agent, X, Y);
    
    % Perform Mitotic Event
    [bone, change_cell, internal_clock] = mitosis_event(bone, change_cell, ...
        internal_clock, site, row_agent, col_agent, row_mitosis, col_mitosis);
               
end 

% Apoptosis Event
if change_cell(row_agent, col_agent) == -1 
    
    % Define Apoptotic Site    
    [jj3, kk3] = define_apoptotic_site(bone, row_target, col_target, ...
                                             directionx, directiony, X, Y, site);
    
    % Perform Apoptotic Event
    [bone, change_cell, internal_clock] = apoptosis_event(bone, change_cell, ...
        internal_clock, site, row_agent, col_agent, jj3, kk3);

end
       