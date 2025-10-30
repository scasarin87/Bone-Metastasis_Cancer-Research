bone_marrow_mask = bone(:, :) == 0;
cortical_bone_mask = bone(:, :) == -1;
outer_mask = bone(:, :) == -2;