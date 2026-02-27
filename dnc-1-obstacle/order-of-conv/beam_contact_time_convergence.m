function beam_contact_time_convergence()
% beam_contact_time_convergence.m
%
% Goal: Measure the **temporal** order of convergence of the implicit
% Euler–Bernoulli contact solver by refining Δt and using the **finest**
% time step as a reference ("approximate exact").
%
% What this program does
% ----------------------
% 1) Fix a sufficiently fine spatial grid (Δx small) and a final time T.
% 2) Choose a list of time steps Δt_k = Δt_0 / 2^(k-1). For each Δt_k:
%    - Solve the time-dependent contact problem (implicit FD).
%    - Extract the right-tip trace u_tip^k(t) = u(1,t).
% 3) Use the finest run (smallest Δt_ref) as the reference solution.
%    Since Δt_k are powers of 2, coarse time nodes are a subset of the
%    fine grid: we directly subsample the reference on coarse times.
% 4) For each Δt_k, compute temporal errors over [0,T] using the tip trace:
%       L_inf_t = max_j |u_tip^k(t_j) - u_tip^ref(t_j)|
%       L2_t    = ( sum_j |...|^2 * Δt_k )^{1/2}
% 5) Print a rate table with two-grid observed orders:
%       p_k = log(E_k / E_{k+1}) / log(Δt_k / Δt_{k+1})
% 6) Plot ONLY the log–log error curves vs Δt.
%
% Notes:
% - We measure **time** convergence (not space). Keep Δx fixed and small.
% - Contact makes solutions nonsmooth in time near impact; orders may drop.
% - No least-squares slopes; no solution-visualization plots.

%% ------------------------ User settings ------------------------
L      = 1.0;        % beam length
T      = 10.0;        % final time horizon for the test
dx     = 0.001;       % fixed spatial step (fine enough so spatial error is negligible)
alpha2 = 1.0;        % stiffness parameter
f0     = -0.2;       % distributed forcing (constant)
y_minus= -0.02;      % obstacle height
kappa  = 1e-2;       % contact stiffness
beta   = 1.0;        % contact damping
phi    = @(tt) 0.0;  % left displacement BC
u0     = @(xx) 0.0;  % initial displacement
v0     = @(xx) 0.0;  % initial velocity

% Time-step list (powers of 2 so coarse times subset fine times)
dt_list = [0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125];   % refine in time
numK    = numel(dt_list);

%% ------------------------ Storage -----------------------------
L2t_err   = zeros(numK,1);   % L2-in-time error (tip trace)
Linf_terr = zeros(numK,1);   % Linf-in-time error (tip trace)
dt_vals   = dt_list(:);

tip_traces = cell(numK,1);
time_grids = cell(numK,1);

%% ------------------------ Runs at each Δt ----------------------
for k = 1:numK
    dt = dt_list(k);
    [~, t, u] = solve_beam_contact_implicit_FD(L, T, dx, dt, alpha2, f0, y_minus, kappa, beta, phi, u0, v0);
    time_grids{k} = t;              % row vector 1x(N+1)
    tip_traces{k} = u(end, :).';    % column vector (N+1)x1, tip at x=1
end

%% ------------------------ Reference & errors -------------------
% Finest run is the last (smallest dt)
t_ref   = time_grids{end};
tip_ref = tip_traces{end};

for k = 1:numK
    dt  = dt_list(k);
    t_k = time_grids{k};
    uk  = tip_traces{k};

    % Since dt_ref divides dt_k, coarse times align with fine indices:
    m = round(dt / dt_list(end));     % integer ratio
    idx = 1:m:numel(t_ref);           % matching indices on the fine grid
    tref_sub = t_ref(idx).';
    uref_sub = tip_ref(idx);

    % Sanity: sizes should match uk
    if numel(uk) ~= numel(uref_sub)
        error('Time grids do not align as expected. Ensure dt_list are powers of 2.');
    end

    e = uk - uref_sub;
    Linf_terr(k) = max(abs(e));
    L2t_err(k)   = sqrt(sum(e.^2) * dt);   % discrete L2 in time
end

%% ------------------------ Two-grid rates -----------------------
p_L2t   = nan(numK,1);
p_Linft = nan(numK,1);
for k = 1:numK-1
    p_L2t(k)   = log(L2t_err(k)   / L2t_err(k+1))   / log(dt_vals(k) / dt_vals(k+1));
    p_Linft(k) = log(Linf_terr(k) / Linf_terr(k+1)) / log(dt_vals(k) / dt_vals(k+1));
end

%% ------------------------ Print table -------------------------
fprintf('\n%-8s %-10s %-16s %-16s %-10s %-10s\n', ...
    'Index','dt','L2_time_error','Linf_time_error','p_L2(t)','p_Linf(t)');
fprintf('%s\n', repmat('-',1,78));
for k = 1:numK
    p1 = p_L2t(k);   p2 = p_Linft(k);
    s1 = rate_str(p1); s2 = rate_str(p2);
    fprintf('%-8d %-10.6f %-16.6e %-16.6e %-10s %-10s\n', ...
        k, dt_vals(k), L2t_err(k), Linf_terr(k), s1, s2);
end

%% ------------------------ Log–log plot only -------------------
figure('Name','Temporal convergence (tip trace)'); hold on; grid off;
loglog(dt_vals, L2t_err,   'o--', 'LineWidth', 1.6, 'DisplayName','L^2 in time (tip)');
loglog(dt_vals, Linf_terr, 's--', 'LineWidth', 1.6, 'DisplayName','L^\infty in time (tip)');
set(gca,'XDir','reverse'); % optional: finest on the right or left (remove if undesired)
xlabel('$\Delta t$','Interpreter','latex');
ylabel('Error (tip trace)','Interpreter','latex');
title('Log--log plot of temporal error at the right tip','Interpreter','latex');
legend('Location','southwest');

end % ===== end of main =====

%% ===============================================================
%% Helper: format rate strings
function s = rate_str(p)
    if isnan(p)
        s = '--';
    else
        s = sprintf('%8.3f', p);
    end
end

%% ===============================================================
%% Parameterized implicit FD solver (based on your routine)
function [x,t,u] = solve_beam_contact_implicit_FD(L,T,dx,dt,alpha2,f0,y_minus,kappa,beta,phi,u0,v0)
% Implicit second-order FD for Euler–Bernoulli beam with unilateral contact.
% Unknowns at each step: w = [u_1^{j+1},...,u_M^{j+1}]^T.
% BC: u(0,t)=phi(t), u_x(0,t)=0 (ghost), u_xx(1,t)=0, and u_xxx(1,t)
%      from contact law when in contact (else 0).
    M = round(L/dx);
    N = round(T/dt);
    x = linspace(0,L,M+1).';
    t = linspace(0,T,N+1);
    u = zeros(M+1, N+1);

    gamma = alpha2*(dt^2)/(dx^4);

    % Initial data
    u(:,1) = u0(x);
    u(:,2) = u0(x) + dt*v0(x);

    % Enforce left Dirichlet on known time rows we touch
    u(1,:) = phi(t);

    for j = 1:(N-1)
        % Contact state at right end based on u_M^j (current known time t_j)
        in_contact = (u(M+1, j+1) <= y_minus);

        A = spalloc(M, M, 5*M);  % banded
        b = zeros(M,1);

        % Row i=1 (uses u_{-1}=u_1 via u_x(0,t)=0)
        A(1,1) = 1 + 7*gamma;
        if M>=2, A(1,2) = -4*gamma; end
        if M>=3, A(1,3) =  gamma;   end
        b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1)); % phi=0 here

        % Interior rows i=2..M-2
        for i = 2:(M-2)
            r = i;
            if i-2 >= 1
                A(r, i-2) = A(r, i-2) + gamma;
            else
                b(r) = b(r) - gamma*phi(t(j+1));
            end
            A(r, i-1) = A(r, i-1) - 4*gamma;
            A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
            A(r, i+1) = A(r, i+1) - 4*gamma;
            A(r, i+2) = A(r, i+2) + gamma;

            b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
        end

        % Row i=M-1
        if M>=2
            A(M-1, M-3) = gamma;
            A(M-1, M-2) = -4*gamma;
            A(M-1, M-1) = 1 + 5*gamma;
            A(M-1, M  ) = -2*gamma;
            b(M-1) = dt^2*f0 + 2*u(M, j+1) - u(M, j);
        end

        % Row i=M: contact vs no-contact at right end
        if ~in_contact
            A(M, M-2) =  2*gamma;
            A(M, M-1) = -4*gamma;
            A(M, M  ) =  1 + 2*gamma;
            b(M)      =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
        else
            % -alpha^2 u_xxx(1) = kappa (u(1)-y_-) + beta * (u(1)^{j+1}-u(1)^j)/dt
            c = (dx^3/alpha2)*(kappa + beta/dt);
            A(M, M-2) =  2*gamma;
            A(M, M-1) = -4*gamma;
            A(M, M  ) =  1 + 2*gamma*(c + 1);
            b(M) = dt^2*f0 ...
                 + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt)) * u(M+1, j+1) ...
                 - u(M+1, j) ...
                 + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
        end

        % Solve and write back
        w = A \ b;
        u(1,   j+2) = phi(t(j+2));
        u(2:end,j+2) = w;
    end
end
