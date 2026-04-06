function beam_contact_coupled_convergence()
% beam_contact_coupled_convergence.m
%
% Purpose
% --------
% Observe the **combined** order of convergence when refining space and time
% together for the Euler–Bernoulli beam with unilateral contact.
% We refine with Δt_k = c * Δx_k^2 so that the scheme parameter
%     gamma = alpha^2 * (Δt^2)/(Δx^4)
% remains approximately constant across grids.
%
% What it does
% ------------
% 1) Choose nested meshes M_k (powers of 2): Δx_k = L/M_k.
% 2) For each k, set Δt_k = c * Δx_k^2, then snap to N_k steps so that the
%    final time is exactly T (dt_k_eff = T/N_k).
% 3) Solve the time-dependent contact problem (implicit FD).
% 4) Use the **finest** mesh as the reference; compare u(·,T) at coincident nodes:
%       L∞(Ω): max_i |u_k(x_i,T) - u_ref(x_i,T)|
%       L2(Ω): ( Σ_i |...|^2 * Δx_k )^{1/2}
% 5) Report two-grid observed orders:
%       p = log(E_k/E_{k+1}) / log(Δx_k/Δx_{k+1})
% 6) Plot ONLY log–log error vs Δx (no other plots, no LS slope).
%
% Notes
% -----
% - The presence of contact (impact) can reduce observed smoothness and thus
%   the effective order (often to ~1 in time-dominated regimes).
% - Here we keep gamma ~ const, so the observed rate reflects the **coupled**
%   space+time refinement with Δt ∝ Δx^2.

%% -------------------- User settings --------------------
L      = 1.0;        % beam length
T      = 5.0;        % final time
alpha2 = 1.0;        % stiffness parameter
f0     = -0.2;       % distributed forcing (constant)
y_minus= -0.02;      % obstacle level
kappa  = 1e-2;       % contact stiffness
beta   = 1.0;        % contact damping

phi = @(tt) 0.0;     % left displacement BC
u0  = @(xx) 0.0;     % initial displacement
v0  = @(xx) 0.0;     % initial velocity

% Meshes (powers of 2 so nodes nest): coarse -> fine
M_list = [10, 20, 40, 80, 160, 320];   % extend to 320, 640 if runtime allows
numK   = numel(M_list);

% Choose Δt = c * Δx^2 (pick c so N isn't too big)
c_dt = 0.5;  % you can adjust; smaller makes time steps smaller

%% -------------------- Storage --------------------------
dx_vals   = zeros(numK,1);
dt_vals   = zeros(numK,1);
L2_err    = zeros(numK,1);   % L2 over x at final time
Linf_err  = zeros(numK,1);   % Linf over x at final time
profiles  = cell(numK,1);    % u(:,end) for each mesh
x_grids   = cell(numK,1);

%% -------------------- Solve for each (Δx_k, Δt_k) ------
for k = 1:numK
    M  = M_list(k);
    dx = L / M;
    % target dt via c*dx^2, then snap to N steps so that t_N = T exactly
    dt_target = c_dt * dx^2;
    N = max(1, round(T / dt_target));
    dt = T / N;  % effective dt so the final time aligns exactly

    dx_vals(k) = dx;
    dt_vals(k) = dt;

    [x, t, u] = solve_beam_contact_implicit_FD(L, T, dx, dt, ...
        alpha2, f0, y_minus, kappa, beta, phi, u0, v0);

    x_grids{k}  = x;           % (M+1)x1
    profiles{k} = u(:, end);   % spatial profile at final time T
end

%% -------------------- Reference & spatial errors --------
% Finest grid (last one) as reference
x_ref  = x_grids{end};
u_refT = profiles{end};

for k = 1:numK
    xk  = x_grids{k};
    ukT = profiles{k};
    % Because M_list are powers of 2, coarse nodes are a subset of fine nodes:
    stride = (M_list(end) / M_list(k));   % integer
    idx    = 1:stride:numel(x_ref);       % fine indices matching coarse nodes
    uref_sub = u_refT(idx);

    ek = ukT - uref_sub;
    Linf_err(k) = max(abs(ek));
    L2_err(k)   = sqrt(sum(ek.^2) * dx_vals(k));   % discrete L2 over [0,L]
end

%% -------------------- Two-grid rates (w.r.t. Δx) -------
p_L2   = nan(numK,1);
p_Linf = nan(numK,1);
for k = 1:numK-1
    p_L2(k)   = log(L2_err(k)  / L2_err(k+1))  / log(dx_vals(k) / dx_vals(k+1));
    p_Linf(k) = log(Linf_err(k)/ Linf_err(k+1))/ log(dx_vals(k) / dx_vals(k+1));
end

%% -------------------- Print table -----------------------
fprintf('\n%-8s %-10s %-10s %-14s %-14s %-10s %-10s\n', ...
    'M','dx','dt','L2_error','L_inf_error','p_L2','p_Linf');
fprintf('%s\n', repmat('-',1,92));
for k = 1:numK
    s1 = rate_str(p_L2(k)); s2 = rate_str(p_Linf(k));
    fprintf('%-8d %-10.6f %-10.6f %-14.6e %-14.6e %-10s %-10s\n', ...
        M_list(k), dx_vals(k), dt_vals(k), L2_err(k), Linf_err(k), s1, s2);
end

%% -------------------- Log–log plot ONLY -----------------
figure('Name','Coupled space–time convergence @ final time'); hold on; grid off;
loglog(dx_vals, L2_err,  'o--','LineWidth',1.6,'DisplayName','L^2 (space) @ T');
loglog(dx_vals, Linf_err,'s--','LineWidth',1.6,'DisplayName','L^\infty (space) @ T');
xlabel('$\Delta x$','Interpreter','latex');
ylabel('Error at $t=T$','Interpreter','latex');
title('Log--log plot: error vs. $\Delta x$ with $\Delta t \propto \Delta x^2$','Interpreter','latex');
legend('Location','southwest');

end % ================== end main =========================

%% =======================================================
%% Helper: format rate strings
function s = rate_str(p)
    if isnan(p), s = '--'; else, s = sprintf('%8.3f', p); end
end

%% =======================================================
%% Implicit FD contact solver (parameterized)
function [x,t,u] = solve_beam_contact_implicit_FD(L,T,dx,dt,alpha2,f0,y_minus,kappa,beta,phi,u0,v0)
% Implicit second-order FD for Euler–Bernoulli beam with unilateral contact.
% Unknowns at each step: w = [u_1^{j+1},...,u_M^{j+1}]^T.
% BC: u(0,t)=phi(t), u_x(0,t)=0 (ghost), u_xx(1,t)=0, and u_xxx(1,t)
%     from contact law when in contact (else 0).

    M = round(L/dx);
    N = round(T/dt);
    x = linspace(0,L,M+1).';
    t = linspace(0,T,N+1);
    u = zeros(M+1, N+1);

    gamma = alpha2*(dt^2)/(dx^4);

    % Initial data
    u(:,1) = u0(x);
    u(:,2) = u0(x) + dt*v0(x);

    % Enforce left Dirichlet on all touched rows
    u(1,:) = phi(t);

    for j = 1:(N-1)
        in_contact = (u(M+1, j+1) <= y_minus);

        A = spalloc(M, M, 5*M);
        b = zeros(M,1);

        % i = 1 (uses u_{-1}=u_1)
        A(1,1) = 1 + 7*gamma;
        if M>=2, A(1,2) = -4*gamma; end
        if M>=3, A(1,3) =  gamma;   end
        b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1)); % phi=0

        % interior: i=2..M-2
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

        % i = M-1
        if M>=2
            A(M-1, M-3) = gamma;
            A(M-1, M-2) = -4*gamma;
            A(M-1, M-1) = 1 + 5*gamma;
            A(M-1, M  ) = -2*gamma;
            b(M-1) = dt^2*f0 + 2*u(M, j+1) - u(M, j);
        end

        % i = M: contact vs no-contact
        if ~in_contact
            A(M, M-2) =  2*gamma;
            A(M, M-1) = -4*gamma;
            A(M, M  ) =  1 + 2*gamma;
            b(M)      =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
        else
            % -alpha^2 u_xxx(1) = kappa (u(1)-y_-) + beta*(u(1)^{j+1}-u(1)^j)/dt
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
