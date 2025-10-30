%% Display Mask
%  This function is used to plot the masks, cortical bone, bone marrow,
%  tumor and osteoblasts, to show also the individual masks

% cortical bone mask
figure('Position',[20 40 450 250])
imagesc(cortical_bone_mask)
colorbar

% bone marrow mask
figure('Position',[20 400 450 250])
imagesc(bone_marrow_mask)
colorbar

% tumor mask
figure('Position',[750 40 450 250])
imagesc(outer_mask)
colorbar