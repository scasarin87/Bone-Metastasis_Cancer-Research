%% Get Osteoclasts Probability %%

%  This function is used to compute the probability of finding OCs at
%  different heights of the cortical bone depending on the distance form
%  the epiphysis. Values are based on in house experimental data.

%  NB: At the moment we assume an exponential decrease of probability from.
%  epiphysis to diaphysis. The exponential function is fitted on 3 
%  experimental points where tha actual OCs coverage percentage has been 
%  computed

% Input  -> rows   : overall number of rows for current bone geometry
%           min_row: heighest row of the bone marrow mask = 1
%           avg_row: average row of the bone marrow mask = 1
%           max_row: lowest point of the bone marrow mask = 1
%
% Output -> ocs_probability: vector containing the probability between 0
%                            and 1 of the

function [ocs_probability, bone_section, bone_type] = get_ocs_probability(rows, min_row, avg_row, max_row)
    
    % Call global variables for ocs distribution in femur
    global top_slice_ocs mid_slice_ocs low_slice_ocs;
    % Call global variable for ocs distribution in vertebra/calvaria
    global calvaria_ocs_min calvaria_ocs_max vertebra_ocs_min vertebra_ocs_max;
    % Call global variable that defines current bone geometry
    global current_geometry;
    global cell_line;
    
    % Define whether geometry is tibia or femur
    bone_geometry = split(current_geometry, '_');
    bone_type = bone_geometry{1};
    bone_section = bone_geometry{2};
    
    % If i have a bone longitudinal section
    if strncmp(bone_section, 'ls', numel('ls'))
        % Build vector that will contain a probability for each row
        ocs_probability = zeros(rows, 1);
        % Zero-based row count
        zero_min_row = min_row - min_row; 
        zero_avg_row = avg_row - min_row;
        row_values = linspace(zero_min_row, zero_avg_row, zero_avg_row + 1);
        zero_max_row = max_row - min_row;
        % Define parameters to fit exponential curve to in vivo data
        get_final_parameters;
        % Define probability value for each row of the ABM
        if strncmp(bone_type, 'tibia', numel('tibia'))
            % Only higher epiphysis has lots of OCs
            row_values = linspace(zero_min_row, zero_max_row, zero_max_row + 1);
            ocs_probability(min_row : max_row) = fun(final_params, row_values);
            ocs_probability = ocs_probability ./ 2; % Because OCs are larger than one site
        elseif strncmp(bone_type, 'femur', numel('femur'))
            % Both epiphysis have lots of OCs
            row_values = linspace(zero_min_row, zero_avg_row, zero_avg_row + 1);
            high_ocs_prob = fun(final_params, row_values);
            low_ocs_prob = flip(high_ocs_prob);
            ocs_probability(min_row : avg_row) = high_ocs_prob;
            ocs_probability(avg_row : avg_row + size(low_ocs_prob, 2) - 1) = low_ocs_prob;
            ocs_probability = ocs_probability ./ 2; % Because OCs are larger than one site
        elseif  strncmp(bone_type, 'calvaria', numel('calvaria'))
            ocs_probability = calvaria_ocs_min + (calvaria_ocs_max - calvaria_ocs_min) * rand;
        else
            error('Bone type not supported')        
        end 
    % If i have a bone cross section   
    elseif strncmp(bone_section, 'cs', numel('cs'))
        if strncmp(bone_type, 'vertebra', numel('vertebra'))
           ocs_probability = vertebra_ocs_min + (vertebra_ocs_max - vertebra_ocs_min) * rand;
        else % Fix coverage to a value
            if strcmp(cell_line, 'pc3')
                ocs_probability = 20;
            else 
                ocs_probability = 6;
            end
        end
    else 
        error('Bone type not supported')
    end
end