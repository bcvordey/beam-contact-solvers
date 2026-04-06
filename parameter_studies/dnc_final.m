clc; close all; clear;

%% --- Description ---
% This code solves the free vibration of a 1D beam using a Finite Difference approach.
% It demonstrates how an obstacle boundary condition is enforced at the beam's free end.

%% --- User Data ---
L    = 1;            % Total length of the beam (m)
dx   = 0.025;        % Step size for the spatial grid (m)
x    = 0:dx:L;       % Vector of node locations
M    = length(x);    % Number of spatial nodes

T    = 2;            % Final computation time (s)
dt   = 0.000004;     % Time step (s)
t    = 0:dt:T;       % Vector of time points
N    = length(t);    % Number of time steps

% Special Parameters
beta      = 1;
kappa     = 0.01;
obstacle  = -0.05;
Tol       = 1e-10;
force     = -0.01;

% Aluminum beam with rectangular cross section (example dimensions)
A         = 1;        % Cross-sectional area (m^2)
rho       = 1;        % Density (kg/m)
E         = 1;        % Young's modulus (Pa)
I         = 1;        % Second moment of area (m^4)

alpha     = sqrt((E*I)/(rho*A*L^4));  % Wave speed-like parameter

disp('lambda must be less than or equal 1 for stability.');
lambda    = dt / ( (dx^2) * (alpha) );  % Courant-like number

% Additional coefficients
tau1      = alpha^2 / dx^3;
tau2      = beta / dt;

phi1      = 2 * tau1 + kappa - tau2;
phi2      = -tau1 - 2*kappa + 2*tau2;

%% --- Boundary & Initial Conditions ---
% Displacement array: u(i,j) = displacement at node i and time step j
u = zeros(M, N);

% Boundary conditions at x = 0 (first two nodes)
u(1,:) = 0;   
u(2,:) = u(1,:);

% Initial condition for displacement (a cubic polynomial shape)
P        = 3 * 0.1 / (L^3);       % For a 10-cm upward displacement at the free end
u(:,1)   = P * (x.^2) .* (3*L - x) / 6;
u(:,2)   = u(:,1);

% For plotting
ymax     = max(abs(u(:,1)));

%% --- Main Time-Stepping Loop ---
for j = 3:N   
    for i = 3:M-2
        % Finite difference update for interior nodes
        u(i,j) = -(lambda^2) * ( u(i-2,j-1) - 4*u(i-1,j-1) - 4*u(i+1,j-1) + ...
                     u(i+2,j-1) + 6*u(i,j-1) ) ...
                 - u(i,j-2) + 2*u(i,j-1) + (dt^2) * force;
    end
    
    % Boundary conditions near the obstacle (end of beam, i.e., node M)
    if (abs(u(M,j-1) - obstacle) < Tol) || (u(M,j-1) < obstacle)
        % Contact / obstacle condition
        u(M-1,j) = ( phi1 * u(M-2,j) - tau1 * u(M-3,j) - ...
                     tau2 * u(M,j-1) + kappa * obstacle ) / phi2;
        
        u(M,j)   = 2 * u(M-1,j) - u(M-2,j);  % Equivalent to reflection
    else
        % Regular free-end condition
        u(M-1,j) = 2 * u(M-2,j) - u(M-3,j);
        u(M,j)   = 3 * u(M-2,j) - 2 * u(M-3,j);
    end
end

%% --- Visualization ---
% % Uncomment the following lines if you want to create a video:
% video_filename = 'beam_vibration.mp4'; 
% v = VideoWriter(video_filename);
% v.FrameRate = 10;
% open(v);

figure;
h = plot(x, u(:,1), '-', 'LineWidth', 5);
xlabel('x');
ylabel('Displacement, u(x,t)');
title('Vibration of a Beam with Obstacle');
axis([0 1.1*L -1.1*ymax 1.1*ymax]);
hold on;
% Add obstacle visualization at x = L
line([L, L], [-1.1*ymax, obstacle], 'Color', 'k', 'LineWidth', 10);
hold off;

% Update plot in time
for kk = 2:1000:N
    set(h, 'YData', u(:,kk));
    drawnow;
    
    % % To record a video, uncomment:
    % frame = getframe(gcf);
    % writeVideo(v, frame);
end

% % If recording video, close the file:
% close(v);

%% --- Post-Processing Measurements --

% 2) Distance between the obstacle and the end of the beam, |obstacle - u(M,t)|
dist_obsM   = abs(obstacle - u(M,:));
dObsM_sub   = dist_obsM(1:500:end);
t_sub  = t(1:500:end);               % sub-sample time

figure;
plot(t_sub, dObsM_sub, 'b-', 'LineWidth', 1, 'MarkerIndices', 1:length(t_sub));
xlabel('Time (s)');
ylabel('|obstacle - u(M,t)|');
title('Distance between the Obstacle and the Beam End (Every 500 Steps)');
grid on;

% 3) Distance between 0 and u(M,t), simply |u(M,t)|
dist0M      = u(M,:);
dist0M_sub  = dist0M(1:500:end);


figure;
plot(t_sub, dist0M_sub, '-', 'LineWidth', 1, 'MarkerIndices', 1:length(dist0M_sub));
xlabel('Time (s)');
ylabel('u(M,t)');
title('Beam End Displacement Over Time (Every 500 Steps)');
grid on;
