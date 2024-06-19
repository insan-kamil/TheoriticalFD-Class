clear; tic;

L = 20; c0 = 1; ds = 0.05;
t_values = 1:-0.01:0.01;

% Video Writer Settings
v = VideoWriter('Ex1_RungeKutta_insan.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 20;  % Adjust for smoother animation
open(v);

% Figure Setup
figure;  
view(30, 40); % Set a 3D view angle
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

for t = t_values
    % ... (rest of the code for solving the Runge-Kutta equations remains the same) ...
    x0 = 2 * c0 * sqrt(t) * [0, 0, 1];
    t0 = [1 0 0]; % T is Tangent Vector
    n0 = [0 1 0]; % n is Normal Vector
    b0 = [0 0 1]; % b is Binormal Vector
    mmax = round(L / ds);
    XXpos = zeros(mmax + 1, 3);
    TTpos = zeros(mmax + 1, 3);
    nnpos = zeros(mmax + 1, 3);
    bbpos = zeros(mmax + 1, 3);
    XXpos(1,:) = x0;
    TTpos(1,:) = t0;
    nnpos(1,:) = n0;
    bbpos(1,:) = b0;

    X = x0; T = t0; n = n0; b = b0; s = 0;
    for m = 1:mmax
        k1x = T;
        k1t = (c0 / sqrt(t)) * n;
        k1n = -(c0 / sqrt(t)) * T + (s / (2 * t)) * b;
        k1b = -(s / (2 * t)) * n;
        Saux = s + ds / 2;
        Xaux = X + (ds / 2) * k1x;
        Taux = T + (ds / 2) * k1t;
        naux = n + (ds / 2) * k1n;
        baux = b + (ds / 2) * k1b;

        k2x = Taux; 
        k2t = (c0 / sqrt(t)) * naux; 
        k2n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k2b = -(Saux / (2 * t)) * naux; 

        Saux = s + ds / 2;
        Xaux = X + (ds / 2) * k2x;
        Taux = T + (ds / 2) * k2t;
        naux = n + (ds / 2) * k2n;
        baux = b + (ds / 2) * k2b;

        k3x = Taux; 
        k3t = (c0 / sqrt(t)) * naux; 
        k3n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k3b = -(Saux / (2 * t)) * naux; 
        
        Saux = s + ds;
        Xaux = X + ds * k3x;
        Taux = T + ds * k3t;
        naux = n + ds * k3n;
        baux = b + ds * k3b;

        k4x = Taux; 
        k4t = (c0 / sqrt(t)) * naux; 
        k4n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k4b = -(Saux / (2 * t)) * naux; 

        s = Saux; % s = m * ds;
        X = X + (ds / 6) * (k1x + 2 * k2x + 2 * k3x + k4x);
        T = T + (ds / 6) * (k1t + 2 * k2t + 2 * k3t + k4t);
        n = n + (ds / 6) * (k1n + 2 * k2n + 2 * k3n + k4n);
        b = b + (ds / 6) * (k1b + 2 * k2b + 2 * k3b + k4b);

        T = T / norm(T);
        n = n / norm(n);
        b = b / norm(b);
        
        % Storing the variables 
        XXpos(m + 1,:) = X;
        TTpos(m + 1,:) = T;
        nnpos(m + 1,:) = n;
        bbpos(m + 1,:) = b;
    end

    % For the negative values
    XXneg = zeros(mmax + 1, 3);
    TTneg = zeros(mmax + 1, 3);
    nnneg = zeros(mmax + 1, 3);
    bbneg = zeros(mmax + 1, 3);
    XXneg(1,:) = x0;
    TTneg(1,:) = t0;
    nnneg(1,:) = n0;
    bbneg(1,:) = b0;

    X = x0; T = t0; n = n0; b = b0; s = 0;
    ds = -ds; % To calculate the negative values
    for m = 1:mmax
        k1x = T;
        k1t = (c0 / sqrt(t)) * n;
        k1n = -(c0 / sqrt(t)) * T + (s / (2 * t)) * b;
        k1b = -(s / (2 * t)) * n;
        Saux = s + ds / 2;
        Xaux = X + (ds / 2) * k1x;
        Taux = T + (ds / 2) * k1t;
        naux = n + (ds / 2) * k1n;
        baux = b + (ds / 2) * k1b;

        k2x = Taux; 
        k2t = (c0 / sqrt(t)) * naux; 
        k2n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k2b = -(Saux / (2 * t)) * naux; 

        Saux = s + ds / 2;
        Xaux = X + (ds / 2) * k2x;
        Taux = T + (ds / 2) * k2t;
        naux = n + (ds / 2) * k2n;
        baux = b + (ds / 2) * k2b;

        k3x = Taux; 
        k3t = (c0 / sqrt(t)) * naux; 
        k3n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k3b = -(Saux / (2 * t)) * naux; 
        
        Saux = s + ds;
        Xaux = X + ds * k3x;
        Taux = T + ds * k3t;
        naux = n + ds * k3n;
        baux = b + ds * k3b;

        k4x = Taux; 
        k4t = (c0 / sqrt(t)) * naux; 
        k4n = -(c0 / sqrt(t)) * Taux + (Saux / (2 * t)) * baux; 
        k4b = -(Saux / (2 * t)) * naux; 

        s = Saux; % s = m * ds;
        X = X + (ds / 6) * (k1x + 2 * k2x + 2 * k3x + k4x);
        T = T + (ds / 6) * (k1t + 2 * k2t + 2 * k3t + k4t);
        n = n + (ds / 6) * (k1n + 2 * k2n + 2 * k3n + k4n);
        b = b + (ds / 6) * (k1b + 2 * k2b + 2 * k3b + k4b);

        T = T / norm(T);
        n = n / norm(n);
        b = b / norm(b);
        
        % Storing the variables 
        XXneg(m + 1,:) = X;
        TTneg(m + 1,:) = T;
        nnneg(m + 1,:) = n;
        bbneg(m + 1,:) = b;
    end
    
    % Reset ds to positive
    ds = -ds;
    
    % Combining all data together
    XX = [XXneg(end: -1:2,:);XXpos];
    TT = [TTneg(end: -1:2,:);TTpos];
    nn = [nnneg(end: -1:2,:);nnpos];
    bb = [bbneg(end: -1:2,:);bbpos];
    
    % Plotting
    clf; % Clear the figure    
    % Plotting with Enhancements
    plot3(XX(:, 1), XX(:, 2), XX(:, 3), 'LineWidth', 1.5, 'Color', 'red'); % Thicker red line
    
    % Dynamically adjust axes limits
    xlim([min(XX(:,1))-1, max(XX(:,1))+1]); 
    ylim([min(XX(:,2))-1, max(XX(:,2))+1]);
    zlim([min(XX(:,3))-1, max(XX(:,3))+1]);

    title(sprintf('Runge-Kutta (VFE) \n Solution at t = %.2f', t));

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
toc;
