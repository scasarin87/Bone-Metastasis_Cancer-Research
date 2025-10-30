%% Get Final Parameters %%

%  This function provides for each geometry the parameters that allows to 
%  to fit the experimental data with an exponentially decreasing curve

%% This code is used to tune the parameters fitting
% x_data = [zero_min_row, zero_min_row + 0.1 * zero_avg_row, zero_avg_row];
% y_data = [top_slice_ocs, mid_slice_ocs, low_slice_ocs];
% % Perform the curve fitting using lsqcurvefit
fun = @(x,x_data)x(1)*exp(x(2)*x_data)+x(3);
% start_params = [100, -100, 1]; % change initial values
% final_params = lsqcurvefit(fun, start_params, x_data, y_data);

%% Now the correct final parameters for each geometry can be added here

global current_geometry;

if strncmp(current_geometry, 'tibia_ls_1', numel('tibia_ls_1')) || strncmp(current_geometry, 'calvaria_ls_1', numel('calvaria_ls_1'))
    final_params = [36.8158 -0.0247 3.1747];
elseif strncmp(current_geometry, 'tibia_ls_1_lr', numel('tibia_ls_1_lr'))
    final_params = [36.8158 -0.0426 3.1747];
elseif strncmp(current_geometry, 'femur_ls_1', numel('femur_ls_1'))
    final_params = [36.8158 -0.0286 3.1747];
elseif strncmp(current_geometry, 'femur_ls_1_lr', numel('femur_ls_1_lr'))
    final_params = [36.8158 -0.0446 3.1747];    
end