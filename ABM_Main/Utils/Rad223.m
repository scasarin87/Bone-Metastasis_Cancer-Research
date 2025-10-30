function R = Rad223(start_t, t, Tau, Rad_max)

    % Time dependent function of Rad activity
    t = (t - start_t + 1)/24;
    R = Rad_max * exp(-t / Tau);
    
end