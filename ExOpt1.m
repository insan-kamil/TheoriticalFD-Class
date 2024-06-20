function fft_derivative_analysis()
    % Number of discretization points
    N = 256;
    
    % Discretize the interval [0, 2*pi)
    x = 2 * pi * (0:N-1) / N; % x_i
    
    % Define the function
    f = exp(cos(x));
    
    % Create the vector of frequency domain multipliers and compute the FFT of the function
    ik = 1i * [0:N/2-1 0 -N/2+1:-1];
    f_ = fft(f); % hat f_k
    f_(abs(f_) / N < eps) = 0; % Krasny's filter
    
    % Compute the numerical derivatives using inverse FFT
    derivatives = compute_numerical_derivatives(f_, ik);
    
    % Compute the exact derivatives
    exact_derivatives = compute_exact_derivatives(x);
    
    % Calculate the errors
    errors = calculate_errors(derivatives, exact_derivatives);
    
    % Display the errors
    display_errors(errors);
end

function derivatives = compute_numerical_derivatives(f_, ik)
    % Compute the numerical derivatives using inverse FFT
    derivatives.fx = ifft(ik .* f_);
    derivatives.fxx = ifft(ik .^ 2 .* f_);
    derivatives.fxxx = ifft(ik .^ 3 .* f_);
    derivatives.fxxxx = ifft(ik .^ 4 .* f_);
    derivatives.fxxxxx = ifft(ik .^ 5 .* f_);
    derivatives.fxxxxxx = ifft(ik .^ 6 .* f_);
end

function exact_derivatives = compute_exact_derivatives(x)
    % Compute the exact derivatives
    exact_derivatives.fx = -sin(x) .* exp(cos(x));
    exact_derivatives.fxx = (sin(x).^2 - cos(x)) .* exp(cos(x));
    exact_derivatives.fxxx = 0.5 * (1 + 6 * cos(x) + cos(2 * x)) .* sin(x) .* exp(cos(x));
    exact_derivatives.fxxxx = 1/8 * (-1 - 4 * cos(x) + 24 * cos(2 * x) + 12 * cos(3 * x) + cos(4 * x)) .* exp(cos(x));
    exact_derivatives.fxxxxx = -1/8 * (31 + 100 * cos(x) + 96 * cos(2 * x) + 20 * cos(3 * x) + cos(4 * x)) .* sin(x) .* exp(cos(x));
    exact_derivatives.fxxxxxx = -1/32 * (34 - 148 * cos(x) + 191 * cos(2 * x) + 630 * cos(3 * x) + 254 * cos(4 * x) + 30 * cos(5 * x) + cos(6 * x)) .* exp(cos(x));
end

function errors = calculate_errors(numerical, exact)
    % Calculate the errors between numerical and exact derivatives
    errors.fx = norm(numerical.fx - exact.fx, inf);
    errors.fxx = norm(numerical.fxx - exact.fxx, inf);
    errors.fxxx = norm(numerical.fxxx - exact.fxxx, inf);
    errors.fxxxx = norm(numerical.fxxxx - exact.fxxxx, inf);
    errors.fxxxxx = norm(numerical.fxxxxx - exact.fxxxxx, inf);
    errors.fxxxxxx = norm(numerical.fxxxxxx - exact.fxxxxxx, inf);
end

function display_errors(errors)
    % Display the errors for each derivative order
    fprintf(['Error for the first derivative is %e\n' ...
             'Error for the second derivative is %e\n' ...
             'Error for the third derivative is %e\n' ...
             'Error for the fourth derivative is %e\n' ...
             'Error for the fifth derivative is %e\n' ...
             'Error for the sixth derivative is %e\n'], ...
             errors.fx, ...
             errors.fxx, ...
             errors.fxxx, ...
             errors.fxxxx, ...
             errors.fxxxxx, ...
             errors.fxxxxxx);
end

% Run the analysis
fft_derivative_analysis();
