function [pdf_value] = custom_pdf(x, limit, type, init)

    % Define your custom PDF function here
    
    if strcmp(type, 'exp')
        if strcmp(init, 'true')
            y1 = 0.001;
        elseif strcmp(init, 'false')
            y1 = 0.08;
        end
        % Exponential decrease between x=0 and x=15, bounded between 0.1 and 0.9
        % Define the probabilities and x-values for the two points
        x1 = 0; x2 = limit;
        y2 = 0.9;
        % Calculate the rate parameter (lambda) for the exponential distribution
        lambda = -log(y2 / y1) / (x2 - x1);
        % Calculate the PDF using the exponential distribution formula
        pdf_value = y1 * exp(-lambda * (x - x1));
    
    elseif strcmp(type, 'sig')
        % Define the sigmoid function
        sigmoid = @(params, x) params(1) + (params(2) - params(1)) ./ (1 + exp(-(x - params(3)) / params(4)));
        % Define the initial guess for the parameters
        initialGuess = [0.01, 0.9, 7.5, 1];
        % Define the x and y data points
        xData = [0, limit];
        yData = [0.001, 0.8];
        % Perform the curve fitting using lsqcurvefit
        paramsFit = lsqcurvefit(sigmoid, initialGuess, xData, yData);
        % Calculate the corresponding y value for the fitted curve
        pdf_value = sigmoid(paramsFit, x);
    else
        error('wrong pdf type')
    end

end






