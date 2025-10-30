%% Define Tumor Masks %%

%  This function is used to define the current tumor mask as the enseble of
%  Pca cells and blood vessels and its complementary mask, which contains
%  all the other ABM sites (CB, BM, OBs, OCs, Outer)

%  Input  -> bone               : current bone matrix
%            site               : struct defining the ABM sites

%  Output -> tumor           : tumor mask
%            no_tumor        : complementary to tumor mask

function [tumor, no_tumor] = def_tumor_masks(bone, site)
    
    % Define Tumor Mask Matrix
    tumor = (bone == site.tumor_edge | bone == site.tumor | ...
                        bone == site.vessel | bone == site.vessel_cabo);
    % Define the Tumor Complementary Mask Matrix
    no_tumor = (bone == site.outer | bone == site.cortical_bone | ...
                bone == site.bone_marrow | bone == site.osteoblast | ...
                bone == site.osteoclast);

end