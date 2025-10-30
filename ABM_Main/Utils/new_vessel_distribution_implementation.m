%% New Vessel Distribution Implementation

%  This script is used to implement in the model a new distribution of
%  blood veesels obtained from experimental data

% Define path and sheet name
path = strcat(pwd, "\ABM Main", "\PC3_data.xlsx");
sheetName = 'DatiLuca'; % Replace with the name of your sheet

% Read data from the specified sheet into a table
df = readtable(path, 'Sheet', sheetName);

% Cross Section vessels major axis
perp_max_length = df.('MaxPerpend');
nanMask = isnan(perp_max_length);
perp_max_length = perp_max_length(~nanMask);
% Get the semi-axis
perp_max_length = perp_max_length / 2;
% Cross Section vessels minor axis
perp_min_length = df.('MinPerpend');
nanMask = isnan(perp_min_length);
perp_min_length = perp_min_length(~nanMask);
% Get the semi-axis
perp_min_length = perp_min_length / 2;
% Longitudinal vessels major axis
par_max_length = df.('MaxParall');
nanMask = isnan(par_max_length);
par_max_length = par_max_length(~nanMask);
% Get the semi-axis
par_max_length = par_max_length / 2;
% Longitudinl vessels minor axis
par_min_length = df.('MinParall');
nanMask = isnan(par_min_length);
par_min_length = par_min_length(~nanMask);
% Get the semi-axis
par_min_length = par_min_length / 2;
% Distance between vessels
distance = df.('Distance');
nanMask = isnan(distance);
distance = distance(~nanMask);

% Create a Kernel Distribution object
pd_a_kernels = fitdist(perp_max_length, 'Kernel', 'Kernel', 'normal');  
pd_a_kernels_par = fitdist(par_max_length, 'Kernel', 'Kernel', 'normal');  

% Plot distributions
x_values = 0:0.5:300;
pdf_values_old = pdf(pd_a_kernels, x_values);
pdf_values_new = pdf(pd_a_kernels_par, x_values);
plot(x_values, pdf_values_old);
hold on
plot(x_values, pdf_values_new);
xlabel('Dimension [um]');
ylabel('Probability');
title('Kernel Distribution of Major Semi-Axis');
legend('C4-2B', 'PC3', 'FontSize',10)

         
% Estimate the parameters (mu and lambda) of the Inverse Gaussian distribution
params_perp = fitdist(perp_min_length, 'InverseGaussian');
% Create the Inverse Gaussian distribution object
pd_b_InverseGaussian = makedist('InverseGaussian', 'mu', params_perp.mu, 'lambda', params_perp.lambda);

% Estimate the parameters (mu and lambda) of the Inverse Gaussian distribution
params_par = fitdist(par_min_length, 'InverseGaussian');
% Create the Inverse Gaussian distribution object
pd_b_InverseGaussian_par = makedist('InverseGaussian', 'mu', params_par.mu, 'lambda', params_par.lambda);
          
% Plot the PDF
x_values = 0:0.25:100;  % Adjust the range as needed
pdf_values_old = pdf(pd_b_InverseGaussian, x_values);
pdf_values_new = pdf(pd_b_InverseGaussian_par, x_values);

figure
plot(x_values, pdf_values_old);
hold on
plot(x_values, pdf_values_new);
xlabel('Distance [um]');
ylabel('Probability');
title('Inverse Gaussian Distribution of Minor Semi-Axis');
legend('C4-2B', 'PC3', 'FontSize',10)

pd_Distance = fitdist(distance, 'Normal');

