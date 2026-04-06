% fdm_contact.m   — explicit finite‑difference scheme with unilateral contact
% -------------------------------------------------------------------------
clear;  clc;

%% Parameters (identical to the Python version)
L  = 1;          % domain length
T  = 20;         % total time
M  = 10;         % spatial divisions   (=> M+1 nodes)
N  = 10000;      % time steps          (=> N+1 time levels)
c  = 2;
f  = 1;

beta            = 0.1;
obstaclePosition = 1;   % “h” in the paper
kappa           = 0.1;

initialVelocity = 1;
t0 = (1/f) * ( initialVelocity ...
      - sqrt(initialVelocity^2 + 2*f*obstaclePosition) );

Tol = 1.0e-5;

%% Derived grid sizes
dx = L / M;
dt = (T - t0) / N;

phi = 1 / (1 + kappa*dx + beta*dx/dt);   % Eq. (ψ)
lambda = (c*dt)/dx;


%% Allocate and set initial data
u = zeros(N+1, M+1);         % u(j,i)  with j=1..N+1,  i=1..M+1
u(1,:) = 0;                   % u^0_i  = 0
u(2,:) = initialVelocity * dt; % u^1_i  = v0 * Δt   (row‑2 in MATLAB)

%% Time‑marching loop
for j = 2:N                     % j = 2 … N   (corresponds to Python j = 1 … N-1)
    % Internal nodes: i = 2 … M   (Python 1 … M-1)
    for i = 2:M
        u(j+1,i) = 2*(1 - c^2*(dt/dx)^2) * u(j,i) ...
                  + c^2*(dt/dx)^2 * (u(j,i+1) + u(j,i-1)) ...
                  - u(j-1,i) + f * dt^2;
    end
    
    % Left boundary (Neumann, Eq. 62)
    u(j+1,1) = u(j+1,2);
    
    % Right boundary (contact, Eq. 63)
    tip_now = u(j, M+1);               % u_M^j   (M+1 due to MATLAB indexing)
    if abs(tip_now - obstaclePosition) < Tol || tip_now > obstaclePosition
        u(j+1, M+1) = phi * ( ...
              u(j+1, M) ...                    % u_{M-1}^{j+1}
            + kappa * obstaclePosition * dx ...
            + beta * (dx/dt) * tip_now ...
        );
    else
        u(j+1, M+1) = u(j+1, M);
    end
end

%% Plotting — identical visual to the Python version
time_vec = linspace(t0, T, N+1);

ult = u(:,1) - L;                        % left end, shifted like Python

figure;  hold on;
plot( u(:, M+1), time_vec, 'r', 'DisplayName','right end','LineWidth',1.2 );
plot( ult,        time_vec, 'b', 'DisplayName','left end','LineWidth',1.2 );
xline( obstaclePosition, '--k', 'Obstacle', 'DisplayName','Obstacle','LineWidth',1.2 );

xlabel('Position');
ylabel('Time');
title('End‑Point Displacements vs Time (Finite‑Difference Simulation)');
legend('Location','best');
set(gca,'YLim',[0, T - t0]); 
grid on;


rod_len = u(:, M+1) - u(:,1) + L;        % u_right + L − u_left

figure;                                  % new window
plot(time_vec, rod_len,'k', 'LineWidth',1.2);
xlabel('time');
ylabel('Length');
title('Length of the rod over time');
set(gca,'XLim',[0, T - t0]); 
grid on;
