%% ABM OBs Response to Rad223 %%

%  This function is used to reduce the osteoblasts number when Rad223
%  theraphy starts (and then makes them increase again following the
%  percentages computed from in vivo observations)

% This algorithm performs an obs decrease between the theraphy start and the
% half life time of rad223 (this parameters can be modified if needed) and
% then it makes them increase from day 12 to day 30 (assuming that from
% that moment on the obs number is back to the original one)

% checking whether the period of time is < 30 days 
if (hour >= start_therapy_rad) && (hour <= start_therapy_rad + 30 * 24)  

    % time passed since Rad223 administration (considering it in days)
    time_rad = (hour - start_therapy_rad) / 24;  

    % current osteoblast number
    obs_number = sum(sum(bone == site.osteoblast));

    % interpolation of survival percentage (PCHIP interpolation) to do so the
    % normlaization of the percentages in fraction is performed 
    time_points = [0, 7, 11, 18, 25, 30];  % Time points (in days)
    percentages = [100, 81.87, 70.42, 94.16, 94.06, 100] / 100;

    percentage = interp1(time_points, percentages, time_rad, 'pchip');

    % forcing 100% recovery after 30 days --> the in vivo observation show a
    % recovery after day 25 (the percentages are all close to each other
    % from that moment on)
    if time_rad > 30
        percentage = 1.0;
    end

    % compute target osteoblast number
    target_obs_number = fix(percentage * obs_number_start);
    
    if target_obs_number < obs_number  
        % if number of obs > number of obs to remove than we are at t < day
        % 11 and the obs have to be removed
        obs_to_eliminate = obs_number - target_obs_number;

        for obs = 1:obs_to_eliminate
            % finding osteoblast coordinates
            [row_obs, col_obs] = find(bone == site.osteoblast);

            % the ob is replaced by a special agent that make possible to
            % cell to stay empty instead that to be occupied by a oc
            if ~isempty(row_obs)
                random_obs = randi(length(row_obs));
                bone(row_obs(random_obs), col_obs(random_obs)) = site.cortical_bone_induced;
            end
        end

    elseif target_obs_number > obs_number
        % if the number of obs < number of obs that should be present than
        % we are at t > day 11 and this means that Radium effects have
        % anded and obs need to be reinserted within the model
        obs_to_add = target_obs_number - obs_number;

        for obs = 1:obs_to_add
            % the cell occupied is one of the ones with the special agent
            % that keeps it empty from ocs
            [row_cb, col_cb] = find(bone == site.cortical_bone_induced);

            if ~isempty(row_cb)
                random_cb = randi(length(row_cb));
                bone(row_cb(random_cb), col_cb(random_cb)) = site.osteoblast;
            end
        end
    end
end

