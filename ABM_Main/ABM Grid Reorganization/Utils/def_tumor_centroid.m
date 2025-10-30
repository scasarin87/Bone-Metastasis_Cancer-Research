%% Define Tumor Centroid %%

%  This function is used to compute the tumor centroid

%  Input  -> tumor      : mask of PCa cells and blood vessels

%  Output -> row/col_centroid : tumor row and column centroid

function [row_centroid, col_centroid] = def_tumor_centroid(tumor)
    %row_centroid = round(regionprops(tumor, 'centroid')(1).Centroid(2));
    %col_centroid = round(regionprops(tumor, 'centroid').Centroid(1));
    props = regionprops(tumor, 'centroid');
    row_centroid = round(props(1).Centroid(2));
    col_centroid = round(props(1).Centroid(1));
end
