clear
close all

% Parameters
M = 4; % Number of sides of the polygon
N = 64 * M; % Total number of points
s = 2 * pi * (0:N-1) / N;
T = zeros(N, 3); % 3D tension vector, initialized to zero

% Setup initial tension vectors
for m = 0:M-1
    idx = m * (N/M) + (1:N/M);
    T(idx, 1) = cos(2 * pi * m / M);
    T(idx, 2) = sin(2 * pi * m / M);
end

% Initial position vector X
X = (2 * pi / N) * cumsum([0 0 0; T(1:end-1, :)]);
tend = 2 * pi / M^2;
mmax = 2400;
dt = tend / mmax;

% Initialize data storage
Xcorner = zeros(mmax + 1, 3); % Evolution of a corner
Xcorner(1, :) = X(1, :);
nplots = 120;
Xdata = zeros(N, 3, nplots + 1);
Tdata = zeros(N, 3, nplots + 1);
idata = 1;

% Store initial data
Xdata(:, :, idata) = X;
Tdata(:, :, idata) = T;

% Fourier transform wavenumbers
ik = 1i * [0:N/2-1 -N/2:-1]';

tic
for m = 1:mmax
    % Runge-Kutta integration
    [K1X, K1T] = computeRKStep(T, ik);
    [K2X, K2T] = computeRKStep(T + (dt/2) * K1T, ik);
    [K3X, K3T] = computeRKStep(T + (dt/2) * K2T, ik);
    [K4X, K4T] = computeRKStep(T + dt * K3T, ik);
    
    % Update position and tension vectors
    X = X + (dt/6) * (K1X + 2 * K2X + 2 * K3X + K4X);
    T = T + (dt/6) * (K1T + 2 * K2T + 2 * K3T + K4T);
    
    % Normalize tension vectors
    T = T ./ sqrt(sum(T .^ 2, 2));
    
    % Store corner data
    Xcorner(m + 1, :) = X(1, :);
    
    % Store data for plotting
    if mod(m, mmax/nplots) == 0
        idata = idata + 1;
        Xdata(:, :, idata) = X;
        Tdata(:, :, idata) = T;
    end
    
    % Plot evolution of the corner
    plot(sqrt(Xcorner(1:m,1) .^2 + Xcorner(1:m,2) .^2) + 1i * Xcorner(1:m,3))
    title('Evolution of X(0, t)');
    axis([0, 0.4, 0, 0.4])
    drawnow
end
toc

% Function to compute Runge-Kutta step
function [KX, KT] = computeRKStep(T, ik)
    T_ = fft(T);
    T_(abs(T_)/numel(T) < eps) = 0; % Krasny's filter
    Ts = ifft(ik .* T_);
    Tss = ifft(ik .^ 2 .* T_);
    KX = cross(T, Ts); % cross product
    KT = cross(T, Tss);
end
