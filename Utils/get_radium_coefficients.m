%% Get Rad Coefficients %%

%  This function returns the coefficients of the 3rd order polynomial
%  function fitting the distance-dependent activity of Rad on the base of
%  the intrinsic effect of Rad (intended at the max of its temporal
%  activity) and the standard probabilities associated with our current
%  tumor

function [B_mit, B_apop] = get_radium_coefficients(Rad)

    % 3rd order polynomial function
    pol_degree = 3;

    % Matrix A_temp will contain:
    % 1st row --> distance discretization
    % 2nd row --> prbability of mitosis
    % 3rd row --> probability of apoptosis

    %A_temp = xlsread('Rad_Activity.xls');
    A_temp = [0 100 200 300 400;
              0.0100 0.0200 0.0400 0.0600 0.0600
              0.1300 0.0900 0.0500 0.0400 0.0400]; % data found in Rad_Activity excel table

    % Scaling Factors: we need to scale the effect of Rad such as the inside of
    % the tumor (beyond 300 microns) has the same probailities of the Control
    % both in terms of mitosis and apoptosis
    epsilon = zeros(1, 2); % scaling factors vector
    epsilon(1) = Rad.mitosis / A_temp(2, end);
    epsilon(2) = Rad.apoptosis / A_temp(3, end);

    % Apply the scaling factors to the intrinsic distribution of probabilities
    % associate to Rad and get the right distance-dependent distribution of
    % probabilities of Rad-223 effect
    newA_temp(1,:) = A_temp(2, :) * epsilon(1);
    newA_temp(2,:) = A_temp(3, :) * epsilon(2);

    % Get the coefficients for the polynomical functioncs that will fit the
    % distance-dependent Rad activity into the model
    B_mit  = polyfit(A_temp(1, :), newA_temp(1, :), pol_degree);
    B_apop = polyfit(A_temp(1, :), newA_temp(2, :), pol_degree);

end