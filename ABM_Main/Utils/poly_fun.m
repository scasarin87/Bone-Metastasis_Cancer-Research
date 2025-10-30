function output = poly_fun(d, A)

    % 3rd order polynomial function for the distance dependet probability of
    % mitosis/apoptosis
    if d > 300
        d = 300; % beyond 400 microns we have a plateau
    end

    output = A(1) * d.^3 + A(2) * d.^2 + A(3) * d + A(4);
    
end