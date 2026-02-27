% Finite-Difference Simulation with Right-End Contact (MATLAB)
clear; clc;

%% Parameters (match Python)
L = 1;            % domain length
T = 20;           % total time
M = 10;           % number of spatial subintervals (M+1 nodes)
N = 10000;        % number of time steps (will be adjusted only if needed)
c = 2;            % wave speed
f = 1;            % constant forcing

% Contact/obstacle parameters
beta = 1;
obstaclePosition = 1;
kappa = 10;

% Initial velocity and t0 (same formula as Python)
initialVelocity = 1;
t0 = (1/f)*(initialVelocity - sqrt(initialVelocity^2 + 2*f*obstaclePosition));

% Tolerance
Tol = 1e-5;

%% Discretization
dx = L / M;
dt = (T - t0) / N;

% Enforce CFL for convergence (explicit 1D wave scheme: r <= 1)
r = c*dt/dx;
if r > 1
    dt = 0.95*dx/c;                 % shrink dt safely
    N  = ceil((T - t0)/dt);         % adjust number of steps
    dt = (T - t0)/N;                % recompute dt to end exactly at T
    r  = c*dt/dx;                   % update CFL
end
r2 = r^2;

% Contact coefficient (same as Python)
phi = 1 / (1 + kappa*dx + beta*(dx/dt));

%% Allocate and set initial conditions
u = zeros(N+1, M+1);                % (time, space)
u(1,:) = 0;                         % u^0
u(2,:) = initialVelocity * dt;      % u^1 from initial velocity

%% Time stepping
for j = 2:N
    % Interior nodes: central in time & space (vectorized)
    u(j+1,2:M) = 2*(1 - r2)*u(j,2:M) + r2*(u(j,3:M+1) + u(j,1:M-1)) ...
                 - u(j-1,2:M) + f*dt^2;

    % Left boundary (zero-slope Neumann: mirror)
    u(j+1,1) = u(j+1,2);

    % Right boundary with unilateral contact
    if abs(u(j,M+1) - obstaclePosition) < Tol || u(j,M+1) > obstaclePosition
        u(j+1,M+1) = phi * ( u(j+1,M) + kappa*obstaclePosition*dx + beta*(dx/dt)*u(j,M+1) );
    else
        u(j+1,M+1) = u(j+1,M);
    end
end

%% Time vector and auxiliary quantities
time_vec = linspace(t0, T, N+1).';  % column vector for plotting
ult = u(:,1) - L;                   % left end minus L (to match Python)

%% Plot 1: End-point displacements vs time (like the Python orientation)
figure;
plot(u(:,end), time_vec, 'b', 'DisplayName','Right end','LineWidth',1.2); hold on;
plot(ult,       time_vec, 'r', 'DisplayName','Left end','LineWidth',1.2);
xline(obstaclePosition, '--k', 'Obstacle Position','LineWidth',1.2);
xlabel('Position'); ylabel('Time');
title('End-Point Displacements vs Time (Finite-Difference Simulation)');
legend('Location','best'); grid on;

%% Plot 2: Length of the rod over time
rod_len = u(:,end) - u(:,1) + L;    % right − left + initial length
figure;
plot(time_vec, rod_len, 'LineWidth', 1.6);
xlabel('Time'); ylabel('Length');
title('Length of the rod over time');
grid on;

% Optional: show CFL used
fprintf('CFL number r = c*dt/dx = %.4f (<= 1 for stability)\n', r);
