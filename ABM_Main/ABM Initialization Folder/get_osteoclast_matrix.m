%% Get Osteoclast Matrix %%

%  This function generates a matrix of the same shape of the ABM that 
%  identifies the spatial position of a new OCs. The OC shape can range 
%  from 60 to 100um, accordingly it will occupy from 3 to 5 ABM sites

% Input  -> bone   : ABM main matrix. Has labels for all model agents
%           cortical_bone_edge: matrix that identifies the interface 
%                               between cortical bone and bone marrow    
%           row_cortbone_edge: curr row of the selected cortbone_edge site
%           col_cortbone_edge: curr col of the selected cortbone_edge site
%           cb_site: current cortical bone site
%           method: method of the erosion function
%
% Output -> osteoclasts: matrix = 1 where site is the new oc to be added

function [osteoclasts] = get_osteoclast_matrix(bone, cortical_bone_edge, ...
                                               row_cortbone_edge, col_cortbone_edge, ...
                                               cb_site, method);
    
    % Initialize empty matrix
    osteoclasts = ones(size(bone, 1), size(bone, 2));
    % Identify the new OC center
    osteoclasts(row_cortbone_edge(cb_site), col_cortbone_edge(cb_site)) = 0;
    % Erosion to expand OC center
    osteoclasts = imerode(osteoclasts, method);
    osteoclasts = imerode(osteoclasts, method);
    % I want OC site = 1
    osteoclasts = not(osteoclasts);
    % I keep only the sites belonging to the cortical bone edge
    osteoclasts = (osteoclasts & cortical_bone_edge);
    
end