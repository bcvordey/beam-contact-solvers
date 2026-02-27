function beam_dnc_two_obstacles_FDt()
% Compare right-end displacement for multiple beta values on one plot.
% Implicit FD scheme with TWO-SIDED contact at x=1 (lower y_- and upper y_+).
% Time-varying distributed load: f(x,t) = F0*sin(omega*t) (uniform in x).

%% Fixed parameters (mirroring the 1-obstacle driver)
L = 1.0; T = 50.0; dx = 0.01; dt = 0.01;
alpha2 = 1.296;

% --- Obstacles and DNC law ---
y_minus = -0.02;                % lower obstacle
y_plus  =  0.02;                % upper obstacle
kappa   = 0.1;                 % contact stiffness
betas   = [0.001,0.1,1,5];% beta sweep

% --- Time-varying load parameters f(x,t) = F0*sin(omega*t) ---
F0    = 0.6;
omega = 12.0;

fig = figure('Visible','on'); hold on; grid off;
for b = betas
    [t, u_right] = simulate_right_end_twoobs(b, kappa, y_minus, y_plus, ...
        L, T, dx, dt, alpha2, F0, omega);

    % one figure per beta (no stacking)
    

    p = plot(t, u_right, 'LineWidth', 1.6, 'DisplayName', sprintf('\\beta = %g', b));
   
    % % Save as its own PDF (filename includes kappa and beta)
    % kap_str = regexprep(sprintf('%.3g', kappa), '\.', 'p');
    % bet_str = regexprep(sprintf('%g', b),      '\.', 'p');
    % fname   = sprintf('right_end_kappa_%s_beta_%s.pdf', kap_str, bet_str);
    % exportgraphics(fig, fname, 'ContentType','vector');  % or: print(fig, fname, '-dpdf','-painters')
    % 
    % close(fig);  % close figure so the next beta gets a fresh one
end

yline(y_minus, '--','LineWidth',1.6, 'HandleVisibility','off');
yline(y_plus,  '--','LineWidth',1.6, 'HandleVisibility','off');

xlim([0 T]);
xlabel('time');
ylabel('u(1,t)  (right end position)', 'Interpreter','latex');
title(sprintf('Right-end displacement vs time ($\\kappa = %.3g$), $f(t)=%.3g\\sin(%.3gt)$', ...
          kappa, F0, omega), 'Interpreter','latex');

legend('Location','northeast');box on;
hold off;
end

% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end_twoobs(beta, kappa, y_minus, y_plus, ...
    L, T, dx, dt, alpha2, F0, omega)
% Runs one simulation for a given beta; returns time vector and right-end trace.
% This is the 1-obstacle routine with ONLY the top-contact branch added,
% and with a time-dependent load f(t)=F0*sin(omega*t) evaluated at t_{j+1}.

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
u(1,:) = 0;                    % enforce Dirichlet at x=0

for j = 1:(N-1)                % compute u^{j+1} from u^{j}, u^{j-1}

    % -------- time-dependent load evaluated at t_{j+1} (implicit) --------
    fjp1 = F0 * sin(omega * t(j+2));   % j runs 1..N-1 -> j+2 in [3..N+1]

    % Unknowns are w = [u_1^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*fjp1 + 2*u(2, j+1) - u(2, j);  % phi=0 so +4γ*phi term vanishes

    % Interior rows: i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            % would add gamma*phi to RHS; phi=0 -> no effect
        end
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

    % i=M row: three regimes (no contact / bottom / top), decided by u_M^j
    u_tip_prev = u(M+1, j+1);   % this is u_M^j in our storage convention

    if (u_tip_prev > y_minus) && (u_tip_prev < y_plus)
        % ------------------ No contact ------------------
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*fjp1 + 2*u(M+1, j+1) - u(M+1, j);

    else
        % ------------------ Contact (bottom or top) ------------------
        % Same algebra as the 1-obstacle contact row, with y_wall = y_- or y_+.
        if u_tip_prev <= y_minus
            y_wall = y_minus;   % bottom
        else
            y_wall = y_plus;    % top
        end

        c = (dx^3/alpha2) * (kappa + beta/dt);
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);

        b(M) = dt^2*fjp1 ...
             + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             -   u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa*y_wall;
    end

    % Solve and write back
    w = A \ b;
    u(1,   j+2) = 0;       % phi = 0
    u(2:end,j+2) = w;
end

u_right = u(end,:);         % trace at x=1
end
