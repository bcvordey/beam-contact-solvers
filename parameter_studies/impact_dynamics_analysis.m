function restitution()
%% Fixed parameters
L = 1.0;  T = 10.0;  dx = 0.01;  dt = 0.01;
alpha2 = 1;
f0 = -0.2;
y_minus = -0.02;

% ------------------ Beta range (use linspace) ------------------
beta_min = 0.001;
beta_max = 1;
n_beta   = 10;                     % number of beta values (controls "step")
betas    = linspace(beta_min, beta_max, n_beta);

% ------------------ Kappa values ------------------
kappas = [1, 5, 10, 15, 20];

for kk = 1:numel(kappas)
    kappa = kappas(kk);

    % One figure per kappa
    figure('Color','w');
    ax = axes; hold(ax,'on');

    % Use MATLAB's default color order automatically (different color each line)
    colororder(ax, lines(numel(betas)));

    for ib = 1:numel(betas)
        beta = betas(ib);

        [t_event, e_event] = simulate_tip_restitution(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

        % Plot with color only (solid line, no markers)
        if isempty(t_event)
            plot(ax, nan, nan, '-', ...
                'LineWidth', 1.6, ...
                'DisplayName', sprintf('\\beta = %.3g', beta));
        else
            plot(ax, t_event, e_event, '-', ...
                'LineWidth', 1.6, ...
                'DisplayName', sprintf('\\beta = %.3g', beta));
        end
    end

    % Cosmetics
    box(ax,'on');
    grid(ax,'off');
    xlim(ax,[0 T]);

    xlabel(ax,'$t$','Interpreter','latex');
    ylabel(ax,'$e(t)$','Interpreter','latex');

    title(ax, sprintf('Restitution Coefficient at $x=1$ ($\\kappa=%.3g$)', kappa), ...
        'Interpreter','latex');

    legend(ax,'Location','best');
end

end

% -------------------------------------------------------------------------
function [t_exit, e_exit] = simulate_tip_restitution(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Runs one simulation and returns (event time, restitution) at each EXIT event.

[t, U] = simulate_right_end_trace(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

Vhalf = diff(U) / dt;

tol = 1e-10;
isAbove_n   = (U(1:end-1) >  y_minus + tol);
isBelow_np1 = (U(2:end)   <= y_minus + tol);
entryIdx = find(isAbove_n & isBelow_np1);

isBelow_n   = (U(1:end-1) <= y_minus + tol);
isAbove_np1 = (U(2:end)   >  y_minus + tol);
exitIdx  = find(isBelow_n & isAbove_np1);

t_exit = [];
e_exit = [];

p = 1;
for k = 1:numel(entryIdx)
    n_in = entryIdx(k);

    while p <= numel(exitIdx) && exitIdx(p) <= n_in
        p = p + 1;
    end
    if p > numel(exitIdx), break; end

    n_out = exitIdx(p);

    denom_out = (U(n_out+1) - U(n_out));
    if abs(denom_out) < eps
        theta_out = 0.5;
    else
        theta_out = (y_minus - U(n_out)) / denom_out;
        theta_out = min(max(theta_out,0),1);
    end
    t_out = t(n_out) + theta_out*dt;

    v_in  = Vhalf(n_in);
    v_out = Vhalf(n_out);

    if abs(v_in) < 1e-14
        e = NaN;
    else
        e = -v_out / v_in;
    end

    t_exit(end+1,1) = t_out;
    e_exit(end+1,1) = e;

    p = p + 1;
end

end

% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end_trace(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Implicit FD solver for one obstacle; returns tip trace u(1,t).

M  = round(L/dx);
N  = round(T/dt);
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

u0  = @(xx) 0.0;
v0  = @(xx) 0.0;

gamma = alpha2*(dt^2)/(dx^4);

u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0;

for j = 1:(N-1)
    in_contact = (u(M+1, j+1) <= y_minus);

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);

    % interior i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i=M-1
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)      = dt^2*f0 + 2*u(M, j+1) - u(M, j);

    % i=M
    if ~in_contact
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        c = (dx^3/alpha2) * (kappa + beta/dt);
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);
        b(M)      =  dt^2*f0 ...
                   + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
                   -   u(M+1, j) ...
                   + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    w = A \ b;
    u(1,   j+2)   = 0;
    u(2:end,j+2)  = w;
end

u_right = u(end,:).';
end
