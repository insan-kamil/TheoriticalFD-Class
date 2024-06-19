clear; close all; tic;

v = VideoWriter('Ex1_ode45_insan.mp4','MPEG-4');
v.Quality = 100;
v.FrameRate = 20;  % Adjust for smoother animation
open(v);

L = 20;
c0 = 1;
ds = 0.01;
T0 = [1 0 0];
n0 = [0 1 0];
b0 = [0 0 1];
options = odeset('RelTol',5e-14,'AbsTol',5e-14);

figure;  % Create figure outside the loop
view(30, 40); % Set a 3D view angle
xlabel('X'); ylabel('Y'); zlabel('Z'); 
grid on;

for t = 1:-0.01:0.01
    X0 = 2*c0*sqrt(t)*[0 0 1];
    XTnb0 = [X0 T0 n0 b0];

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
    
    XX=[XXneg(end:-1:2,:);XXpos];
    TT=[TTneg(end:-1:2,:);TTpos];
    nn=[nnneg(end:-1:2,:);nnpos];
    bb=[bbneg(end:-1:2,:);bbpos];

    plot3(XX(:,1), XX(:,2), XX(:,3), 'LineWidth', 1.5, 'Color', 'blue');  % Thicker blue line

    title(sprintf('ode45 (VFE) \n Solution at t = %.2f', t));
    
    % Dynamically adjust axes limits
    xlim([min(XX(:,1)), max(XX(:,1))]);
    ylim([min(XX(:,2)), max(XX(:,2))]);
    zlim([min(XX(:,3)), max(XX(:,3))]);

    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);
toc;
