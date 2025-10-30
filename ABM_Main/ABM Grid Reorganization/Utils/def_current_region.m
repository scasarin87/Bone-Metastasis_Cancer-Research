function [] = def_current_region(curr_tumor, row_tum, col_tum)

    centroids = [];
    bound_box = [];
    tum_masks = {};

    centroids = [centroids; regionprops(curr_tumor, 'all').Centroid];
    bound_box = [bound_box; regionprops(curr_tumor, 'all').BoundingBox];
    tum_masks = {tum_masks; regionprops(curr_tumor, 'all').Image};    
    
    
end