%% ABM OCs Response to Rad223 %%

%  This function is used to alter the osteobcasts number when Rad223
%  theraphy starts (and then makes them increase again following the
%  percentages computed from in vivo observations)

% This algorithm implement a specific curve of growth for ocs found trhough
% observational studies 
if hour == 1
    oc_number_start = sum(sum(bone == site.osteoclast));
end 

% Checking whether the period of time is within the observation window
if (hour >= start_therapy_rad) && (hour <= start_therapy_rad + half_life_time_rad223 - 1)

    % Time passed since Rad223 administration (in days)
    time_rad = (hour - start_therapy_rad) / 24;

    % Current osteoclast number
    oc_number = sum(sum(bone == site.osteoclast));

    % Define experimental time points (in days) and normalize OC percentages
    time_points_OC = [1, 2, 7, 11, 18, 25, 32];
    oc_percentages = [100, 164.7, 55.6, 51.3, 28.4, 14.2, 22.5] / 100;

    % Perform interpolation (PCHIP)
    percentage = interp1(time_points_OC, oc_percentages, time_rad, 'pchip');

    % Ensure full recovery after day 32
    if time_rad > 32
        percentage = 1.0;
    end

    % Compute target osteoclast number
    target_oc_number = fix(percentage * oc_number_start);

    if target_oc_number < oc_number  
        % Osteoclasts need to be removed
        oc_to_eliminate = oc_number - target_oc_number;

        for oc = 1:oc_to_eliminate
            % Find OC coordinates
            [row_oc, col_oc] = find(bone == site.osteoclast);

            % If OCs are present, remove them
            if ~isempty(row_oc)
                random_oc = randi(length(row_oc));
                bone(row_oc(random_oc), col_oc(random_oc)) = site.cortical_bone_induced;
            end
        end

    elseif target_oc_number > oc_number
        % Osteoclasts need to be added
        oc_to_add = target_oc_number - oc_number;

        for oc = 1:oc_to_add
            % Find empty spaces where OCs can be added
            [row_empty, col_empty] = find(bone == site.cortical_bone_induced);

            if ~isempty(row_empty)
                random_empty = randi(length(row_empty));
                bone(row_empty(random_empty), col_empty(random_empty)) = site.osteoclast;
            end
        end
    end
end