function beam_dnc1obs_restitution_vs_damping()
% =========================================================================
% ONE-OBSTACLE (lower stop) TIP-CONTACT BEAM (FD + DNC) :
% Compare restitution coefficient vs DNC damping term magnitude.
%
% Outputs per beta:
%   - Tip trace U(t) = u(1,t)
%   - Event-based restitution e_k at rebound (exit) times
%   - Episode damping metrics during contact:
%         peakD_k  = max_{contact} |beta * u_t|
%         meanD_k  = mean_{contact} |beta * u_t|
%         intD_k   = sum_{contact}  |beta * u_t| dt   (integral-like)
%
% Plots:
%   (1) e(t) at exit times for all betas
%   (2) peak |beta*u_t| per episode vs time (same exit times)
%   (3) scatter: e vs peak |beta*u_t| (correlation picture)
% =========================================================================

%% ------------------------ USER PARAMETERS -------------------------------
L = 1.0;     T = 10.0;     dx = 0.01;     dt = 0.01;
alpha2 = 1.0;

% Load (constant or time-varying). Here constant matches your earlier code.
f0 = -0.2;

% Obstacle level and stiffness
y_minus = -0.02;

% Allow kappa vector (we plot one figure per kappa, betas on same axes)
kappas = [1, 10, 100, 1000];

% beta sweep
betas = [1e-2, 1e-1, 1, 5];

% Event hysteresis (use small positive value if chattering occurs)
epsH = 0;    % e.g. 1e-6

% Style per beta (distinct markers)
style = cell(numel(betas),1);
style{1} = struct('ls','-','mk','*','lw',1.5,'ms',6); % 0.01
style{2} = struct('ls','-','mk','o','lw',1.5,'ms',5); % 0.1
style{3} = struct('ls','-','mk','s','lw',1.5,'ms',5); % 1
style{4} = struct('ls','-','mk','d','lw',1.5,'ms',5); % 5
%% ------------------------------------------------------------------------

for kk = 1:numel(kappas)
    kappa = kappas(kk);

    % ==============================================================
    % FIGURE A: restitution e(t) for this kappa (all betas)
    % ==============================================================
    % figE = figure('Color','w');
    % axE = axes(figE); hold(axE,'on'); box(axE,'on'); grid(axE,'off');

    % ==============================================================
    % FIGURE B: episode peak damping metric vs time (all betas)
    % ==============================================================
    figD = figure('Color','w');
    axD = axes(figD); hold(axD,'on'); box(axD,'on'); grid(axD,'off');

    % ==============================================================
    % FIGURE C: scatter e vs peak damping metric (all betas)
    % ==============================================================
    figS = figure('Color','w');
    axS = axes(figS); hold(axS,'on'); box(axS,'on'); grid(axS,'off');

    for ib = 1:numel(betas)
        beta = betas(ib);

        % ---- Solve beam and get tip trace ----
        [t, U] = simulate_right_end_oneobs(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

        % ---- Postprocess: restitution + damping metrics per episode ----
        out = restitution_and_damping_metrics(t, U, y_minus, dt, beta, epsH);

        % out fields:
        %   out.t_exit, out.e_exit, out.peakD, out.meanD, out.intD

        % ---- Plot e(t) at exit times ----
        % figure(figE);
        % if isempty(out.t_exit)
        %     plot(axE, nan, nan, style{ib}.ls, 'Marker', style{ib}.mk, ...
        %         'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
        %         'DisplayName', sprintf('\\beta=%g (no contact)', beta));
        % else
        %     plot(axE, out.t_exit, out.e_exit, style{ib}.ls, 'Marker', style{ib}.mk, ...
        %         'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
        %         'DisplayName', sprintf('\\beta=%g', beta));
        % end

        % ---- Plot peak |beta*u_t| per episode at same exit times ----
        figure(figD);
        if isempty(out.t_exit)
            plot(axD, nan, nan, style{ib}.ls, 'Marker', style{ib}.mk, ...
                'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g (no contact)', beta));
        else
            plot(axD, out.t_exit, out.peakD, style{ib}.ls, 'Marker', style{ib}.mk, ...
                'LineWidth', style{ib}.lw, 'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g', beta));
        end

        % ---- Scatter: e vs peakD (correlation picture) ----
        figure(figS);
        if ~isempty(out.e_exit)
            plot(axS, out.peakD, out.e_exit, 'LineStyle','none', ...
                'Marker', style{ib}.mk, 'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g', beta));
        else
            plot(axS, nan, nan, 'LineStyle','none', ...
                'Marker', style{ib}.mk, 'MarkerSize', style{ib}.ms, ...
                'DisplayName', sprintf('\\beta=%g (no contact)', beta));
        end
    end

    % Finalize FIGURE A
    % figure(figE);
    % xlim(axE, [0 T]);
    % xlabel(axE,'$t$','Interpreter','latex');
    % ylabel(axE,'$e(t)$','Interpreter','latex');
    % title(axE, sprintf(['Restitution coefficient at exit events  ' ...
    %     '($\\kappa=%.3g$, $y_-=%.3g$, $\\Delta x=%.3g$, $\\Delta t=%.3g$)'], ...
    %     kappa, y_minus, dx, dt), 'Interpreter','latex');
    % legend(axE,'Location','best');

    % Finalize FIGURE B
    figure(figD);
    xlim(axD, [0 T]);
    xlabel(axD,'$t$','Interpreter','latex');
    ylabel(axD,'$\max_{\mathrm{contact}}|\beta\,u_t(1,t)|$','Interpreter','latex');
    title(axD, sprintf(['Peak DNC dashpot magnitude per contact episode  ' ...
        '($\\kappa=%.3g$, $y_-=%.3g$)'], kappa, y_minus), 'Interpreter','latex');
    legend(axD,'Location','best');

    % Finalize FIGURE C
    figure(figS);
    xlabel(axS,'$\max_{\mathrm{contact}}|\beta\,u_t(1,t)|$','Interpreter','latex');
    ylabel(axS,'$e$','Interpreter','latex');
    title(axS, sprintf(['Episode-wise comparison (correlation view)  ' ...
        '($\\kappa=%.3g$, $y_-=%.3g$)'], kappa, y_minus), 'Interpreter','latex');
    legend(axS,'Location','best');
end

end

% =========================================================================
%            POSTPROCESS: RESTITUTION + DAMPING METRICS PER EPISODE
% =========================================================================
function out = restitution_and_damping_metrics(t, U, y_minus, dt, beta, epsH)
% Compute:
%   - exit times t_exit
%   - restitution e_exit at exits
%   - damping metrics during each contact episode:
%       peakD = max |beta*ut|
%       meanD = mean |beta*ut|
%       intD  = sum  |beta*ut| dt
%
% We use slab velocities V^{n+1/2} = (U^{n+1}-U^n)/dt.

U = U(:);
Vhalf = diff(U)/dt;                 % length N
absDash = abs(beta * Vhalf);        % |beta*u_t| on each slab

% Entry (above -> below) and exit (below -> above)
entryIdx = find( (U(1:end-1) >  y_minus + epsH) & (U(2:end) <= y_minus - epsH) );
exitIdx  = find( (U(1:end-1) <= y_minus - epsH) & (U(2:end) >  y_minus + epsH) );

t_exit = [];
e_exit = [];
peakD  = [];
meanD  = [];
intD   = [];

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

    % Exit crossing time (linear interpolation)
    denom = (U(n_out+1) - U(n_out));
    if abs(denom) < eps
        theta = 0.5;
    else
        theta = (y_minus - U(n_out)) / denom;
        theta = min(max(theta,0),1);
    end
    t_out = t(n_out) + theta*dt;

    v_in  = Vhalf(n_in);
    v_out = Vhalf(n_out);

    % Restitution (positive on rebound)
    if abs(v_in) < 1e-14
        e = NaN;
    else
        e = -v_out / v_in;
    end

    % Damping metrics across the contact interval in slab indices:
    % slabs correspond to Vhalf(n) on [t_n, t_{n+1}]
    slabs = n_in:n_out; % covers from entry slab through exit slab
    dvals = absDash(slabs);

    peak_val = max(dvals);
    mean_val = mean(dvals);
    int_val  = sum(dvals) * dt;

    t_exit(end+1,1) = t_out;
    e_exit(end+1,1) = e;
    peakD(end+1,1)  = peak_val;
    meanD(end+1,1)  = mean_val;
    intD(end+1,1)   = int_val;

    p = p + 1;
end

out.t_exit = t_exit;
out.e_exit = e_exit;
out.peakD  = peakD;
out.meanD  = meanD;
out.intD   = intD;

end

% =========================================================================
%                   ONE-OBSTACLE BEAM SOLVER (TIP CONTACT)
% =========================================================================
function [t, u_right] = simulate_right_end_oneobs(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Implicit FD scheme with unilateral DNC contact at right tip x=1 (lower stop y_-).
% Returns:
%   t       : time grid (N+1)
%   u_right : tip trace U^n = u(1,t_n), n=0..N

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
    % contact state decided from previous tip value at time level j
    in_contact = (u(M+1, j+1) <= y_minus);

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);

    % interior rows i=2..M-2
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

    % i=M row: no-contact vs contact
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

        b(M) = dt^2*f0 ...
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
