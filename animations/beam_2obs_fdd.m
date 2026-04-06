function beam_dnc_two_obstacles_FDD()
% Compare right-end displacement for multiple beta values on one plot.
% Implicit FD scheme with TWO-SIDED contact at x=1.
% Bottom and top obstacles can have DIFFERENT (kappa, beta).

%% Fixed parameters
L = 1.0; T = 10.0; dx = 0.01; dt = 0.01;
alpha2 = 1.296;

% --- Obstacles and DNC law ---
y_minus   = -0.02;         % lower obstacle
y_plus    =  0.02;         % upper obstacle
kappa_bot = 0.001;             % κ_- (bottom)
kappa_top = 1000;             % κ_+ (top)   <-- choose independently

% Sweep bottom damping β_-; use a fixed β_+ for the top
betas      = 0.001;   % β_- set
beta_top   = 0;                 % β_+ (top)  <-- choose independently

% --- Time-varying load f(x,t) = F0*sin(omega*t) ---
F0    = 0.6;
omega = 12.0;

fig = figure('Visible','on'); hold on; grid on;
for b = betas
    % b is the bottom β_-
    [t, u_right] = simulate_right_end_twoobs( ...
        b, kappa_bot, ...         % bottom β_-, κ_-
        beta_top, kappa_top, ...  % top β_+,   κ_+
        y_minus, y_plus, L, T, dx, dt, alpha2, F0, omega);

    plot(t, u_right, 'LineWidth', 1.6, 'DisplayName', sprintf('\\beta_- = %g', b));
end

% Obstacle lines
yline(y_minus, '--','LineWidth',1.6, 'DisplayName', sprintf('y_- = %.3g', y_minus));
yline(y_plus,  '--','LineWidth',1.6, 'DisplayName', sprintf('y_+ = %.3g', y_plus));

xlim([0 T]);
xlabel('time');
ylabel('u(1,t)  (right end position)', 'Interpreter','latex');
title(sprintf(['Right-end displacement vs time ($\\kappa_- = %.3g$, $\\kappa_+ = %.3g$, ', ...
               '$\\beta_+ = %.3g$), $f(t)=%.3g\\,\\sin(%.3g\\,t)$'], ...
               kappa_bot, kappa_top, beta_top, F0, omega), 'Interpreter','latex');

legend('Location','northeast');
hold off;
end

% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end_twoobs( ...
    beta_bot, kappa_bot, beta_top, kappa_top, ...
    y_minus, y_plus, L, T, dx, dt, alpha2, F0, omega)
% One simulation; two-sided DNC with distinct (beta,kappa) at bottom/top.

M  = round(L/dx);              % nodes: 0..M
N  = round(T/dt);              % times: 0..N
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

phi = @(tt) 0.0;               % u(0,t) = 0
u0  = @(xx) 0.0;               % initial displacement
v0  = @(xx) 0.0;               % initial velocity

gamma = alpha2*(dt^2)/(dx^4);

% Storage (u(i+1,j+1) ~ u_i^j)
u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0;                    % Dirichlet at x=0

for j = 1:(N-1)                % compute u^{j+1} from u^{j}, u^{j-1}

    % Implicit load at t_{j+1}
    fjp1 = F0 * sin(omega * t(j+2));

    % Unknowns: w = [u_1^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*fjp1 + 2*u(2, j+1) - u(2, j);

    % Interior rows: i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1, A(r, i-2) = A(r, i-2) + gamma; end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;
        b(r) = b(r) + dt^2*fjp1 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i=M-1 row
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)     = dt^2*fjp1 + 2*u(M, j+1) - u(M, j);

    % i=M row: decide regime using u_M^j (stored as u(M+1, j+1))
    u_tip_prev = u(M+1, j+1);

    if (u_tip_prev > y_minus) && (u_tip_prev < y_plus)
        % -------- No contact --------
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*fjp1 + 2*u(M+1, j+1) - u(M+1, j);

    else
        % -------- Contact (bottom OR top) with separate (kappa,beta) --------
        if u_tip_prev <= y_minus
            % Bottom contact: use (kappa_bot, beta_bot)
            y_wall = y_minus;
            beta_c = beta_bot;
            kappa_c = kappa_bot;
        else
            % Top contact: use (kappa_top, beta_top)
            y_wall = y_plus;
            beta_c = beta_top;
            kappa_c = kappa_top;
        end

        c = (dx^3/alpha2) * (kappa_c + beta_c/dt);
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);

        b(M) = dt^2*fjp1 ...
             + 2*(1 + gamma*(beta_c*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             -   u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa_c*y_wall;
    end

    % Solve and advance
    w = A \ b;
    u(1,   j+2) = 0;       % phi = 0
    u(2:end,j+2) = w;
end

u_right = u(end,:);         % trace at x=1
end
