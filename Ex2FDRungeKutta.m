clear; close all; tic;

v = VideoWriter('VFE_Ex2_FiniteDifferences_improved.mp4', 'MPEG-4');
v.Quality = 100;
v.FrameRate = 20;  % Adjust for smoother animation
open(v);

L = 20;
c0 = 1;
t = 1;
ds = 0.01;
X0 = 2*c0*sqrt(t)*[0 0 1];
T0 = [1 0 0];
n0 = [0 1 0];
b0 = [0 0 1];
options = odeset('RelTol',5e-14,'AbsTol',5e-14);

% ... (rest of the code for calculating initial conditions is the same) ...
XTnb0=[X0 T0 n0 b0];

[~,XTnb]=ode45(@(s,XTnb)rhsVFE(s,XTnb,c0,t), 0:ds:L, XTnb0, options);
XXpos=XTnb(:,1:3);
TTpos=XTnb(:,4:6);
nnpos=XTnb(:,7:9);
bbpos=XTnb(:,10:12);

[~,XTnb]=ode45(@(s,XTnb)rhsVFE(s,XTnb,c0,t), 0:-ds:-L, XTnb0, options);
XXneg=XTnb(:,1:3);
TTneg=XTnb(:,4:6);
nnneg=XTnb(:,7:9);
bbneg=XTnb(:,10:12);

XX0=[XXneg(end:-1:2,:);XXpos];
TT0=[TTneg(end:-1:2,:);TTpos];
nn0=[nnneg(end:-1:2,:);nnpos];
bb0=[bbneg(end:-1:2,:);bbpos];

% Figure Setup
figure;  
view(30, 40); % Set a 3D view angle
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;

mmax = 50000;
dt = -t/mmax;
XX = XX0;
TT = TT0;
bb = bb0;
nn = nn0;

for m = 1:mmax
    % ... (Runge-Kutta calculations remain the same) ...
    XXs = [zeros(1,3);(XX(3:end,:)-XX(1:end-2,:)) / (2*ds); zeros(1,3)]; %first derivative
    XXss = [zeros(1,3);(XX(3:end,:)-2*XX(2:end-1,:) + XX(1:end-2,:)) / ds^2; zeros(1,3)]; %second derivative
    K1=cross(XXs,XXss); %Runge-Kutta

    XXaux=XX+.5*dt*K1;
    XXs = [zeros(1,3);(XXaux(3:end,:)-XXaux(1:end-2,:)) / (2*ds); zeros(1,3)]; %first derivative
    XXss = [zeros(1,3);(XXaux(3:end,:)-2*XXaux(2:end-1,:) + XXaux(1:end-2,:)) / ds^2; zeros(1,3)]; %second derivative
    K2=cross(XXs,XXss); %Runge-Kutta

    XXaux=XX+.5*dt*K2;
    XXs = [zeros(1,3);(XXaux(3:end,:)-XXaux(1:end-2,:)) / (2*ds); zeros(1,3)]; %first derivative
    XXss = [zeros(1,3);(XXaux(3:end,:)-2*XXaux(2:end-1,:) + XXaux(1:end-2,:)) / ds^2; zeros(1,3)]; %second derivative
    K3=cross(XXs,XXss); %Runge-Kutta

    XXaux=XX+dt*K3;
    XXs = [zeros(1,3);(XXaux(3:end,:)-XXaux(1:end-2,:)) / (2*ds); zeros(1,3)]; %first derivative
    XXss = [zeros(1,3);(XXaux(3:end,:)-2*XXaux(2:end-1,:) + XXaux(1:end-2,:)) / ds^2; zeros(1,3)]; %second derivative
    K4=cross(XXs,XXss); %Runge-Kutta
    
    XX=XX+(dt/6)*(K1+2*K2+2*K3+K4);
    if isnan(XX(2,1)),error('Stability error!'),end
    
    % Plotting and Visualization Improvements
    if mod(m, 500) == 0
        ts = t + m * dt;
        
        plot3(XX(:,1), XX(:,2), XX(:,3), 'LineWidth', 1.5, 'Color', 'green'); % Thicker green line

        % Dynamically adjust axes limits
        xlim([min(XX(:,1))-1, max(XX(:,1))+1]);
        ylim([min(XX(:,2))-1, max(XX(:,2))+1]);
        zlim([min(XX(:,3))-1, max(XX(:,3))+1]); 
        
        title(sprintf('Finite Differences \n Solution at t = %.2f', ts));

        drawnow;
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
end

close(v);
toc;
