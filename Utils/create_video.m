%% Create Video %%

%  This function creates the video of the last simulation

% Color Struct Definition
DefineColors;

% ABM Sites Struct Definition
DefineABM_sites;

% RGB Channels
channels   = 3;

% Set Video Features
v = VideoWriter('nome_video.mp4','MPEG-4');
v.Quality = 96;
open(v);

for hours = 1 : follow_up
    
    figure('Position', [150 150 700 500])
    VideoMatrix = zeros(rows, columns, 3);

    for i = 1 : rows
        for j = 1 : columns
            switch BONE(i, j, hours)
                case site.bone_marrow
                    VideoMatrix(i, j, :) = color.light_grey;
                case site.cortical_bone
                    VideoMatrix(i, j, :) = color.dark_grey;
                case site.cortical_bone_induced
                    VideoMatrix(i, j, :) = color.dark_grey;
                case site.outer 
                    VideoMatrix(i, j, :) = color.white;
                case site.tumor_edge
                    VideoMatrix(i, j, :) = color.light_pink;
                case site.tumor
                    VideoMatrix(i, j, :) = color.mill_pink;
                case site.vessel
                    VideoMatrix(i, j, :) = color.red;
                case site.vessel_cabo
                    VideoMatrix(i, j, :) = color.purple;
                case site.osteoblast
                    VideoMatrix(i, j, :) = color.yellow;
                case site.osteoclast
                    VideoMatrix(i, j, :) = color.ocs_blue;
                case site.mitotic_cell
                    VideoMatrix(i, j, :) = color.green;
                case site.apoptotic_cell
                    VideoMatrix(i, j, :) = color.orange;
                otherwise
                    error('Unexpected ABM site type! Check bone matrix');  
            end  
        end  
    end  

    fig = gcf;
    if (strcmp(current_geometry, 'femur_ls_1') == 1) || (strcmp(current_geometry, 'calvaria_ls_1') == 1)
        set(fig, 'Units', 'pixels', 'Position', [450, 50, 400, 710]);
    elseif (strcmp(current_geometry, 'tibia_ls_1') == 1) || (strcmp(current_geometry, 'tibia_ls_1_lr') == 1)
        set(fig, 'Units', 'pixels', 'Position', [450, 50, 400, 710])
    elseif (strcmp(current_geometry, 'tibia_cs_1') == 1) || (strcmp(current_geometry, 'tibia_cs_2') == 1) || (strcmp(current_geometry, 'tibia_cs_3') == 1) || (strcmp(current_geometry, 'vertebra_cs_1' == 1))
        set(fig, 'Units', 'pixels', 'Position', [450, 80, 600, 500]);
    end

    % Convert hours to days
    days = floor(hours / 24) + 1;
    
    image(VideoMatrix);
    axis off

    % Overlay time text in the top-left corner
    text(2, 10, sprintf('Day: %d', days), 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'FontName', 'Comic Sans Ms');

    frame = getframe;
    axis off
    writeVideo(v, frame);
    close
    pause(0.01)
   
end

close(v);
