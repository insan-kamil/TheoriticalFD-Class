function hypotrochoid_analysis()
    % Define the parameters
    a = 5;
    b = 7;
    c = 2.3;

    % Determine the period
    T = compute_period(a, b);

    % Discretize the time variable
    N = 1000;
    t = linspace(0, T, N);

    % Define the parametric equations of the hypotrochoid
    [x, y] = compute_hypotrochoid(a, b, c, t);

    % Visualize the hypotrochoid
    plot_hypotrochoid(x, y);

    % Compute the length of the hypotrochoid using FFT-based derivatives
    L = compute_hypotrochoid_length(x, y, t, T);

    % Display the period and the length of the hypotrochoid
    fprintf('The period of the hypotrochoid is T = %f\n', T);
    fprintf('The length of the hypotrochoid curve is L = %f\n', L);
end

function T = compute_period(a, b)
    % Compute the period of the hypotrochoid
    T = 2 * pi * b / gcd(a, b);
end

function [x, y] = compute_hypotrochoid(a, b, c, t)
    % Compute the parametric equations of the hypotrochoid
    x = (a - b) * cos(t) + c * cos((a / b - 1) * t);
    y = (a - b) * sin(t) - c * sin((a / b - 1) * t);
end

function plot_hypotrochoid(x, y)
    % Plot the hypotrochoid
    figure;
    plot(x, y, 'LineWidth', 1.5);
    xlabel('x(t)', 'FontSize', 12);
    ylabel('y(t)', 'FontSize', 12);
    title('Hypotrochoid Curve', 'FontSize', 14);
    grid on;
    axis equal;
    set(gca, 'FontSize', 12);
end

function L = compute_hypotrochoid_length(x, y, t, T)
    % Compute the length of the hypotrochoid curve using FFT-based derivatives
    N = length(t);
    dt = t(2) - t(1);

    % FFT of x and y coordinates
    f_x = fft(x);
    f_y = fft(y);

    % Frequency components
    frequencies = 2 * pi * [0:N/2-1 -N/2:-1] / T;

    % Compute derivatives using inverse FFT
    dx_dt = ifft(1i * frequencies .* f_x);
    dy_dt = ifft(1i * frequencies .* f_y);

    % Compute the integrand for the curve length
    integrand = sqrt(real(dx_dt).^2 + real(dy_dt).^2);

    % Numerical integration to find the length
    L = sum(integrand) * dt;
end

% Run the analysis
hypotrochoid_analysis();
