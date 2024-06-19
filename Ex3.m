clear
close all

M = 4; % number of sides of the polygon
N = 64 * M; % Total number of points
s = 2 * pi * (0:N-1) / N;
T = zeros(N, 3); % 3 columns, 3 coordinates of the tension vector, 2D because the last column is 0
tic

v = VideoWriter('VFE_Ex3_Polygon_RationalMultiples.mp4', 'MPEG-4');
v.Quality = 100;
open(v);

for m = 0:M-1
    T(m * (N/M) + (1:N/M), 1) = cos(2 * pi * m / M);
    T(m * (N/M) + (1:N/M), 2) = sin(2 * pi * m / M);
end

X = (2 * pi / N) * cumsum([0 0 0; T(1:end-1, :)]);

tend = 2 * pi / M^2;
mmax = 2400; % Increased mmax for higher resolution
Xcorner = zeros(mmax + 1, 3); % Evolution of a corner
Xcorner(1, :) = X(1, :);
nplots = 120;
X1data = zeros(N, nplots + 1);
X2data = zeros(N, nplots + 1);
X3data = zeros(N, nplots + 1);
T1data = zeros(N, nplots + 1);
T2data = zeros(N, nplots + 1);
T3data = zeros(N, nplots + 1);
idata = 1;
X1data(:, idata) = X(:, 1);
X2data(:, idata) = X(:, 2);
X3data(:, idata) = X(:, 3);
T1data(:, idata) = T(:, 1);
T2data(:, idata) = T(:, 2);
T3data(:, idata) = T(:, 3);

dt = tend / mmax;
ik = 1i * [0:N/2-1 -N/2:-1]';

figure('Name', 'Polygon Evolution', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 100, 800, 600]);

for m = 1:mmax
    T_ = fft(T);
    T_(abs(T_)/N < eps) = 0; % Krasny's filter
    Ts = ifft(ik .* T_);
    Tss = ifft(ik .^2 .* T_);
    K1X = cross(T, Ts); % cross product
    K1T = cross(T, Tss);

    Taux = T + (dt/2) * K1T;
    T_ = fft(Taux);
    T_(abs(T_)/N < eps) = 0; % Krasny's filter
    Ts = ifft(ik .* T_);
    Tss = ifft(ik .^2 .* T_);
    K2X = cross(T, Ts); % cross product
    K2T = cross(T, Tss);

    Taux = T + (dt/2) * K2T;
    T_ = fft(Taux);
    T_(abs(T_)/N < eps) = 0; % Krasny's filter
    Ts = ifft(ik .* T_);
    Tss = ifft(ik .^2 .* T_);
    K3X = cross(T, Ts); % cross product
    K3T = cross(T, Tss);

    Taux = T + dt * K3T;
    T_ = fft(Taux);
    T_(abs(T_)/N < eps) = 0; % Krasny's filter
    Ts = ifft(ik .* T_);
    Tss = ifft(ik .^2 .* T_);
    K4X = cross(T, Ts); % cross product
    K4T = cross(T, Tss);

    X = X + (dt/6) * (K1X + 2 * K2X + 2 * K3X + K4X);
    T = T + (dt/6) * (K1T + 2 * K2T + 2 * K3T + K4T);
    T = T ./ sqrt(T(:, 1) .^ 2 + T(:, 2) .^ 2 + T(:, 3) .^ 2);

    X = real(X);
    T = real(T);
    if isnan(X(1, 1)), error('Stability error!'); end
    Xcorner(m + 1, :) = X(1, :);
    if mod(m, mmax/nplots) == 0
        idata = idata + 1;
        X1data(:, idata) = X(:, 1);
        X2data(:, idata) = X(:, 2);
        X3data(:, idata) = X(:, 3);
        T1data(:, idata) = T(:, 1);
        T2data(:, idata) = T(:, 2);
        T3data(:, idata) = T(:, 3);
    end
    
    % Plot and update the video
    plot3(X(:,1), X(:,2), X(:,3), 'LineWidth', 1.5)
    grid on
    xlabel('X', 'FontSize', 12)
    ylabel('Y', 'FontSize', 12)
    zlabel('Z', 'FontSize', 12)
    title(sprintf('Polygon Evolution (VFE)\nSolution at t = %.2f', m*dt), 'FontSize', 14)
    legend('Polygon Path')
    axis([-2, 2, -2, 2, -1, 1])
    drawnow
    frame = getframe(gcf);
    writeVideo(v, frame);
end
toc
close(v)

% Determine rational multiples of 2*pi/M^2
ratios = [0.0164, 1/4, 1/3, 1/2, 2/3, 3/4, 7/9, 8/9, 1]; % Define the desired rational multiples
timesteps = round(ratios * nplots); % Corresponding timesteps

% Plot results at these timesteps
figure('Name', 'Simulation Results', 'NumberTitle', 'off');
set(gcf, 'Position', [100, 100, 800, 600]);

for i = 1:length(timesteps)
    idx = timesteps(i);
    if idx > nplots + 1
        continue; % Skip indices that exceed the bounds
    end
    subplot(3, 3, i);
    plot3(X1data(:, idx), X2data(:, idx), X3data(:, idx), 'LineWidth', 1.5);
    grid on
    [N, D] = rat(ratios(i)); % Get the numerator and denominator of the ratio
    if D == 1
        titleStr = sprintf('Time = %d * 2\\pi / M^2', N);
    else
        titleStr = sprintf('Time = %d/%d * 2\\pi / M^2', N, D);
    end
    title(titleStr, 'FontSize', 10);
    xlabel('X', 'FontSize', 8);
    ylabel('Y', 'FontSize', 8);
    zlabel('Z', 'FontSize', 8);
    axis equal;
end

sgtitle('Simulation Results at Rational Multiples of 2\\pi/M^2', 'FontSize', 16);



toc
