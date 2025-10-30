%% Compute Bone Marrow Area %%

%  This function computes the difference between the initial and final
%  bone marrow area when the Ocs tumor mediated resorption occurs.

%  Input  -> bone           : matrix defining the abm spatial sites
%            site           : struct defining the abm sites values
%            bm_start_area  : bone marrow area before OCs resorption
%            flag_resorption: if = 1 means the resorption activity occured
%  Output -> diff_area      : (final_bm_area - starting_bm_area) scaled to um2

function [bm_area] = show_bm_area(BONE, site, bm_start_area, show_plot, flag_resorption)
    
    if show_plot
        if flag_resorption
            % Compute normalized bone marrow area wrt its initial size
            for hour = 1 : follow_up 
                bm_area(hour) = sum(sum(BONE(:, :, hour) == site.tumor_edge | BONE(:, :, hour) == site.tumor | BONE(:, :, hour) == site.vessel | BONE(:, :, hour) == site.vessel_cabo | BONE(:, :, hour) == site.bone_marrow)) / bm_start_area;
            end
            clear hour
            
            % Plot of the tumor AREA in time
            figure(1)
            hold on
            plot(1 : 1 : follow_up, bm_area, 'LineWidth', 1.5)
            % Set Labels
            xlabel('Time [h]')
            ylabel('Normalized Reabsorbed Bone Marrow Area [ ]')

            % Set Legend
            legend('Reabsorbed BM Area')        

        else
            bm_area = 0;   
        end 
    end 

end 