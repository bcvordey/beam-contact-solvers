function beam_dnc_two_obs_restitution()
% =========================================================================
% TWO-OBSTACLE EULER-BERNOULLI BEAM (FD + DNC) + RESTITUTION COEFFICIENTS
% =========================================================================
% This single file:
%   1) Solves the two-obstacle beam problem with DNC at x=1 (lower y_- and upper y_+).
%   2) Computes event-based restitution coefficients for:
%        - Lower wall impacts  e_-(t)
%        - Upper wall impacts  e_+(t)
%   3) Plots restitution vs time for multiple beta values.
%
% Restitution per wall is computed from the right-tip trace U^n = u(1,t_n).
% Discrete slab velocity:
%       V^{n+1/2} = (U^{n+1}-U^n)/dt
% Entry/exit events (lower wall y_-):
%       entry: U^n > y_-  and U^{n+1} <= y_-
%       exit : U^n <= y_- and U^{n+1} >  y_-
% Upper wall y_+:
%       entry: U^n < y_+  and U^{n+1} >= y_+
%       exit : U^n >= y_+ and U^{n+1} <  y_+
% Restitution (sign-consistent, positive on rebound):
%       e = - v_out / v_in
% where v_in  = V^{n_in+1/2},  v_out = V^{n_out+1/2}.
%
% Notes:
% - With forcing and compliant contact, e can exceed 1 (energy injected by forcing).
% - If "chattering" near the stops occurs, set epsH > 0 below.
% =========================================================================

%% ------------------------ USER PARAMETERS -------------------------------
L = 1.0;     T = 10.0;     dx = 0.01;     dt = 0.01;
alpha2 = 1.296;

% Stops
y_minus = -0.02;
y_plus  =  0.02;

% Distinct stiffnesses
kappa_minus = 1000;
kappa_plus  = 1000;

% Damping sweep (dashpot parameter)
betas = [0.001, 0.01, 0.1, 1, 5];

% Distributed load: f(x,t) = F0*sin(omega*t) (uniform in x)
F0    = 0.6;
omega = 12.0;

% Event hysteresis/tolerance (use small positive value if needed)
epsH = 1e-6;          % e.g. 1e-6 to avoid toggling when U ~ y_wall

% Plot style per beta (distinct markers)
style = cell(numel(betas),1);
style{1} = struct('ls','-','mk','*','lw',1.5,'ms',6); % 0.001
style{2} = struct('ls','-','mk','^','lw',1.5,'ms',5); % 0.01
style{3} = struct('ls','-','mk','o','lw',1.5,'ms',5); % 0.1
style{4} = struct('ls','-','mk','s','lw',1.5,'ms',5); % 1
style{5} = struct('ls','-','mk','d','lw',1.5,'ms',5); % 5
%% ------------------------------------------------------------------------

% ======================== FIGURE: LOWER WALL ============================
fig1 = figure('Color','w'); ax1 = axes(fig1); hold(ax1,'on');
box(ax1,'on'); grid(ax1,'off');

% ======================== FIGURE: UPPER WALL ============================
fig2 = figure('Color','w'); ax2 = axes(fig2); hold(ax2,'on');
box(ax2,'on'); grid(ax2,'off');

for ib = 1:numel(betas)
    beta = betas(ib);

    % ---- Solve beam and get tip trace U(t) ----
    [t, U] = simulate_right_end_twoobs( ...
        beta, kappa_minus, kappa_plus, y_minus, y_plus, ...
        L, T, dx, dt, alpha2, F0, omega);

    % ---- Compute restitution events for both walls ----
    [tL, eL, tU, eU] = restitution_two_walls(t, U, y_minus, y_plus, dt, epsH);

    % ---- Plot lower restitution (e_-) ----
    if isempty(tL)
        plot(ax1, nan, nan, style{ib}.ls, ...
            'Marker', style{ib}.mk, 'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
            'DisplayName', sprintf('\\beta=%g ', beta));
    else
        plot(ax1, tL, eL, style{ib}.ls, ...
            'Marker', style{ib}.mk, 'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
            'DisplayName', sprintf('\\beta=%g', beta));
    end

    % ---- Plot upper restitution (e_+) ----
    if isempty(tU)
        plot(ax2, nan, nan, style{ib}.ls, ...
            'Marker', style{ib}.mk, 'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
            'DisplayName', sprintf('\\beta=%g ', beta));
    else
        plot(ax2, tU, eU, style{ib}.ls, ...
            'Marker', style{ib}.mk, 'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
            'DisplayName', sprintf('\\beta=%g', beta));
    end
end

% Finalize lower-wall plot
xlim(ax1,[0 T]);
xlabel(ax1,'$t$','Interpreter','latex');
ylabel(ax1,'$e_-(t)$','Interpreter','latex');
title(ax1, sprintf(['Lower-wall restitution at exit events  ' ...
    '($y_-=%.3g$, $\\kappa_-=%.3g$, $\\kappa_+=%.3g$, $F_0=%.3g$, $\\omega=%.3g$)'], ...
    y_minus, kappa_minus, kappa_plus, F0, omega), 'Interpreter','latex');
legend(ax1,'Location','best');

% Finalize upper-wall plot
xlim(ax2,[0 T]);
xlabel(ax2,'$t$','Interpreter','latex');
ylabel(ax2,'$e_+(t)$','Interpreter','latex');
title(ax2, sprintf(['Upper-wall restitution at exit events  ' ...
    '($y_+=%.3g$, $\\kappa_-=%.3g$, $\\kappa_+=%.3g$, $F_0=%.3g$, $\\omega=%.3g$)'], ...
    y_plus, kappa_minus, kappa_plus, F0, omega), 'Interpreter','latex');
legend(ax2,'Location','best');

end

% =========================================================================
%                       RESTITUTION POST-PROCESSING
% =========================================================================
function [tL, eL, tU, eU] = restitution_two_walls(t, U, y_minus, y_plus, dt, epsH)
% Returns event-based restitution coefficients for both walls:
%   (tL,eL) at lower-wall exit times, (tU,eU) at upper-wall exit times.

U = U(:);
Vhalf = diff(U)/dt; % V^{n+1/2} on interval [t_n, t_{n+1}]

% ---- LOWER wall (y_minus): enter from above, exit upward ----
entryL = find( (U(1:end-1) >  y_minus + epsH) & (U(2:end) <= y_minus - epsH) );
exitL  = find( (U(1:end-1) <= y_minus - epsH) & (U(2:end) >  y_minus + epsH) );
[tL, eL] = pair_and_compute(t, U, Vhalf, dt, y_minus, entryL, exitL);

% ---- UPPER wall (y_plus): enter from below, exit downward ----
entryU = find( (U(1:end-1) <  y_plus - epsH) & (U(2:end) >= y_plus + epsH) );
exitU  = find( (U(1:end-1) >= y_plus + epsH) & (U(2:end) <  y_plus - epsH) );
[tU, eU] = pair_and_compute(t, U, Vhalf, dt, y_plus, entryU, exitU);

end

% -------------------------------------------------------------------------
function [tExit, eExit] = pair_and_compute(t, U, Vhalf, dt, yWall, entryIdx, exitIdx)
% Pair each entry with the next exit on the same wall; compute e=-v_out/v_in at exit time.

tExit = [];
eExit = [];

p = 1;
for k = 1:numel(entryIdx)
    n_in = entryIdx(k);

    while p <= numel(exitIdx) && exitIdx(p) <= n_in
        p = p + 1;
    end
    if p > numel(exitIdx)
        break; % entered but did not exit before final time
    end

    n_out = exitIdx(p);

    % Linear interpolation for crossing time inside the exit interval
    denom = (U(n_out+1) - U(n_out));
    if abs(denom) < eps
        theta = 0.5;
    else
        theta = (yWall - U(n_out)) / denom;
        theta = min(max(theta,0),1);
    end
    t_out = t(n_out) + theta*dt;

    v_in  = Vhalf(n_in);
    v_out = Vhalf(n_out);

    % Sign-consistent restitution; positive for rebound
    if abs(v_in) < 1e-14
        e = NaN;
    else
        e = -v_out / v_in;
    end

    tExit(end+1,1) = t_out;
    eExit(end+1,1) = e;

    p = p + 1;
end

end

% =========================================================================
%                         TWO-OBSTACLE BEAM SOLVER
% =========================================================================
function [t, u_right] = simulate_right_end_twoobs(beta, kappa_minus, kappa_plus, ...
    y_minus, y_plus, L, T, dx, dt, alpha2, F0, omega)
% Implicit FD scheme with two-sided DNC at right tip x=1.
% Returns:
%   t        : time grid (N+1)
%   u_right  : tip trace U^n = u(1,t_n), n=0..N

M  = round(L/dx);              % nodes: 0..M
N  = round(T/dt);              % times: 0..N
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

u0  = @(xx) 0.0;
v0  = @(xx) 0.0;

gamma = alpha2*(dt^2)/(dx^4);

% Storage u(i+1,j+1) ~ u_i^j
u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0; % Dirichlet at x=0

for j = 1:(N-1)

    % load at t_{j+1} (implicit)
    fjp1 = F0 * sin(omega * t(j+2));

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*fjp1 + 2*u(2, j+1) - u(2, j);

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

        b(r) = b(r) + dt^2*fjp1 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i=M-1 row
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)      = dt^2*fjp1 + 2*u(M, j+1) - u(M, j);

    % i=M row: two-sided contact decision using previous tip value
    u_tip_prev = u(M+1, j+1);

    if (u_tip_prev > y_minus) && (u_tip_prev < y_plus)
        % No contact
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*fjp1 + 2*u(M+1, j+1) - u(M+1, j);

    else
        % Contact with one wall
        if u_tip_prev <= y_minus
            y_wall = y_minus;
            kappa_local = kappa_minus;
        else
            y_wall = y_plus;
            kappa_local = kappa_plus;
        end

        c = (dx^3/alpha2) * (kappa_local + beta/dt);

        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);

        b(M) = dt^2*fjp1 ...
             + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             -   u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa_local*y_wall;
    end

    w = A \ b;
    u(1,   j+2)   = 0;
    u(2:end,j+2)  = w;
end

u_right = u(end,:).'; % column vector
end
