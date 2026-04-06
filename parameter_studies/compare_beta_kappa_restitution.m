function beam_dnc1obsFD_compare_beta_kappa_restitution()
% Restitution coefficient vs time:
%   - One figure per kappa
%   - All betas plotted on the same axes
%   - Distinct line/marker styles per beta
%   - No grid, box ON

%% Fixed parameters
L = 1.0;  T = 10.0;  dx = 0.01;  dt = 0.01;
alpha2 = 1;
f0 = -0.2;
y_minus = -0.02;

betas  = [1e-3, 1e-2, 1e-1, 0.2,0.5,0.7,0.9, 1];
kappas = [1, 5, 10, 15, 20, 40, 50];   % edit as needed

% ---- Style map for betas (distinct markers/lines) ----
% You can customize these freely.
style = cell(numel(betas),1);
style{1} = struct('ls','-','mk','*','lw',1.5,'ms',6);  % beta=1e-3 : star
style{2} = struct('ls','-','mk','o','lw',1.5,'ms',5);  % beta=1e-2 : circle
style{3} = struct('ls','-','mk','d','lw',1.5,'ms',5);  % beta=1e-1 : diamond
style{4} = struct('ls','-','mk','s','lw',1.5,'ms',5);  % beta=1    : square
style{5} = struct('ls','-','mk','^','lw',1.5,'ms',5);  % beta=5    : triangle
style{6} = struct('ls','-','mk','d','lw',1.5,'ms',5);  % beta=1e-1 : diamond
style{7} = struct('ls','-','mk','s','lw',1.5,'ms',5);  % beta=1    : square
style{8} = struct('ls','-','mk','^','lw',1.5,'ms',5);  % beta=5    : triangle

for kk = 1:numel(kappas)
    kappa = kappas(kk);

    % ---- One figure per kappa ----
    figure('Color','w');
    ax = axes; hold(ax,'on');

    for ib = 1:numel(betas)
        beta = betas(ib);

        [t_event, e_event] = simulate_tip_restitution(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

        if isempty(t_event)
            % show in legend even if empty
            plot(nan, nan, style{ib}.ls, ...
                'Marker', style{ib}.mk, ...
                'LineWidth', style{ib}.lw, ...
                'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g', beta));
        else
            plot(t_event, e_event, style{ib}.ls, ...
                'Marker', style{ib}.mk, ...
                'LineWidth', style{ib}.lw, ...
                'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g', beta));
        end
    end

    % ---- Cosmetics (as requested) ----
    box(ax,'on');
    grid(ax,'off');

    xlim(ax,[0 T]);
    % ylim([0 1]);
    % ylim(ax, [0 1]); % Set y-axis limits for restitution coefficient
    xlabel(ax,'$t$','Interpreter','latex');
    ylabel(ax,'$e(t)$','Interpreter','latex');

    % Clear, descriptive title
    title(ax, sprintf(['Restitution Coefficient at $x=1$  ', ...
        '($\\kappa=%.3g$)'], ...
        kappa), 'Interpreter','latex');

    legend(ax,'Location','best');
end

end


% -------------------------------------------------------------------------
function [t_exit, e_exit] = simulate_tip_restitution(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Runs one simulation and returns (event time, restitution) at each EXIT event.

% ---- Run your beam solver, but only return the tip trace U(t) ----
[t, U] = simulate_right_end_trace(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

% ---- Compute piecewise-linear velocities V^{n+1/2} ----
Vhalf = diff(U) / dt;        % length N, corresponds to interval [t_n, t_{n+1}]
t_mid = t(1:end-1) + 0.5*dt; % midpoint times (optional)

% ---- Detect entry and exit indices (based on crossing y_minus) ----
tol = 1e-10; % you can set tol = 1e-10 for noisy signals
isAbove_n   = (U(1:end-1) >  y_minus + tol);
isBelow_np1 = (U(2:end)   <= y_minus + tol);
entryIdx = find(isAbove_n & isBelow_np1);  % indices n where [n -> n+1] enters

isBelow_n   = (U(1:end-1) <= y_minus + tol);
isAbove_np1 = (U(2:end)   >  y_minus + tol);
exitIdx  = find(isBelow_n & isAbove_np1);  % indices n where [n -> n+1] exits

% ---- Pair each entry with the first subsequent exit ----
t_exit = [];
e_exit = [];

p = 1; % pointer into exitIdx
for k = 1:numel(entryIdx)
    n_in = entryIdx(k);

    % Find next exit after this entry
    while p <= numel(exitIdx) && exitIdx(p) <= n_in
        p = p + 1;
    end
    if p > numel(exitIdx)
        break; % entered but never exited before T
    end

    n_out = exitIdx(p);

    % --- Optional: compute interpolated crossing times within each step ---
    % Entry crossing time:
    denom_in = (U(n_in+1) - U(n_in));
    if abs(denom_in) < eps
        theta_in = 0.5;
    else
        theta_in = (y_minus - U(n_in)) / denom_in;  % in [0,1] ideally
        theta_in = min(max(theta_in,0),1);
    end
    t_in = t(n_in) + theta_in*dt;

    % Exit crossing time:
    denom_out = (U(n_out+1) - U(n_out));
    if abs(denom_out) < eps
        theta_out = 0.5;
    else
        theta_out = (y_minus - U(n_out)) / denom_out;
        theta_out = min(max(theta_out,0),1);
    end
    t_out = t(n_out) + theta_out*dt;

    % --- Velocities on the entry/exit intervals ---
    v_in  = Vhalf(n_in);   % approx u_t on [t_n_in, t_{n_in+1}]
    v_out = Vhalf(n_out);  % approx u_t on [t_n_out, t_{n_out+1}]

    % --- Restitution coefficient (standard positive definition) ---
    % e = - v_out / v_in   (assuming v_in < 0 and v_out > 0)
    if abs(v_in) < 1e-14
        e = NaN; % avoid blow-up
    else
        e = -v_out / v_in;
    end

    % If you *insist* on "velocity_in divided by velocity_out", replace above by:
    % e = v_in / v_out;

    % Store event value at exit time (common choice)
    t_exit(end+1,1) = t_out;
    e_exit(end+1,1) = e;

    p = p + 1; % move to next exit for next entry
end

end

% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end_trace(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Your original solver, slightly modified to return only tip trace u(1,t).

M  = round(L/dx);              % nodes: 0..M
N  = round(T/dt);              % times: 0..N
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
    % contact state at current time level j (index j+1 in storage)
    in_contact = (u(M+1, j+1) <= y_minus);

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);

    % Interior rows: i=2..M-2
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

    % i=M-1 row
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)      = dt^2*f0 + 2*u(M, j+1) - u(M, j);

    % i=M row: contact vs no-contact
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
