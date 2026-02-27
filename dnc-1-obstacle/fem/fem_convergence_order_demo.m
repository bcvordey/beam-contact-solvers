function fem_convergence_order_demo()
% Measure observed spatial order for the FEM tip trace by halving h.
% Uses FEM (C^1 Hermite cubics) + implicit central-difference in time,
% interior viscosity C = epsv*K_b, and a unilateral DNC at the tip.
%
% Rates use: p = log( e(h_{i-1}) / e(h_i) ) / log( h_{i-1} / h_i )

%% Shared parameters (kept consistent with your FDM sample)
L = 1.0;   T = 10.0;   dt = 0.01;
alpha2 = 1;     f0 = -0.2;
y_minus = -0.02; kappa = 1;
beta_dnc = 1;      % pick one beta to test order
epsv = 0.001;      % interior viscosity (set 0 to mimic FDM exactly)

% Mesh sequence (halving each time; last = finest reference)
hs = [0.32, 0.16, 0.08, 0.04, 0.02, 0.01];

% Run solver for each h
Nt = round(T/dt); t = linspace(0,T,Nt+1); %#ok<NASGU>
Utip = cell(numel(hs),1);
for k = 1:numel(hs)
    dx = hs(k);
    [~, u_tip] = fem_tip_trace_single(beta_dnc, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv);
    Utip{k} = u_tip(:);   % column
end

% Reference = finest mesh
Uref = Utip{end};

% Errors vs. reference (∞-norm in time and L2 in time) for all but finest
E_inf = zeros(numel(hs)-1,1);
E_L2  = zeros(numel(hs)-1,1);
for k = 1:numel(hs)-1
    e = Utip{k} - Uref;
    E_inf(k) = max(abs(e));
    E_L2(k)  = sqrt(dt * sum(e.^2));
end

% Observed orders using p = log( e_{i-1}/e_i ) / log( h_{i-1}/h_i )
% Align with rows so that p(1) = NaN (no previous), p(i) uses levels (i-1,i)
p_inf = nan(numel(E_inf),1);
p_L2  = nan(numel(E_L2),1);
for i = 2:numel(E_inf)
    p_inf(i) = pair_rate(E_inf(i-1), E_inf(i), hs(i-1), hs(i));
    p_L2(i)  = pair_rate(E_L2(i-1),  E_L2(i),  hs(i-1), hs(i));
end

% Display (first row has no rate; last error corresponds to hs(end-1))
fprintf('\nConvergence vs finest (reference h=%.5g)\n', hs(end));
fprintf('%10s %12s %12s %10s %10s\n','h','||e||_inf','||e||_L2','p_inf','p_L2');
for k = 1:numel(hs)-1
    pinf_str = val2str(p_inf(k));
    pL2_str  = val2str(p_L2(k));
    fprintf('%10.5g %12.4e %12.4e %10s %10s\n', hs(k), E_inf(k), E_L2(k), pinf_str, pL2_str);
end

% %% Plot error vs h (optional)
% figure; loglog(hs(1:end-1), E_inf, 'o-', 'LineWidth',1.6); hold on; grid on;
% loglog(hs(1:end-1), E_L2,  's-', 'LineWidth',1.6);
% set(gca,'XDir','reverse');  % finer → right
% xlabel('h'); ylabel('error (tip trace)');
% legend('||e||_\infty(t)','||e||_{L^2(t)}','Location','southwest');
% title(sprintf('FEM spatial convergence at tip (\\beta=%.2g, \\epsilon=%.2g)',beta_dnc,epsv));
end

% ================== Helpers below (in the same file) =====================

function r = pair_rate(e_prev, e_curr, h_prev, h_curr)
% p = log(e_prev/e_curr) / log(h_prev/h_curr), robust to zeros/NaNs
    if ~isfinite(e_prev) || ~isfinite(e_curr) || ~isfinite(h_prev) || ~isfinite(h_curr) ...
       || h_prev <= 0 || h_curr <= 0
        r = NaN; return;
    end
    if e_prev==0 && e_curr==0
        r = NaN; return;           % indistinguishable (both zero)
    elseif e_curr==0 && e_prev>0
        r = Inf; return;           % clear super-convergence relative to reference
    elseif e_prev==0 && e_curr>0
        r = -Inf; return;          % pathological (should not happen in refinement)
    end
    r = log(e_prev/e_curr) / log(h_prev/h_curr);
end

function s = val2str(v)
% Format rate for printing; NaN -> '--', Inf -> 'Inf', -Inf -> '-Inf'
    if isnan(v), s = '--';
    elseif isinf(v) && v>0, s = 'Inf';
    elseif isinf(v) && v<0, s = '-Inf';
    else, s = sprintf('%.3f', v);
    end
end

% ======= Your existing solver and assembly helpers (unchanged) ==========
function [t, u_right] = fem_tip_trace_single(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv)
% FEM + implicit central difference + unilateral DNC at the tip (one stop).
% Returns time vector and right-end displacement u(1,t).

% Grids
J  = round(L/dx);   dx = L/J;
Nt = round(T/dt);   dt = T/Nt;
t  = linspace(0,T,Nt+1);

% Assemble FE operators
[M, Kb, F, free, id_tip_w] = assembleHermiteCubic1D(J, dx, f0);
K = alpha2 * Kb;
C = epsv   * Kb;

% Reduced blocks
Mff = M(free,free);
Kff = K(free,free);
Cff = C(free,free);
Ff  = F(free);
e_tip_full = zeros(size(M,1),1); e_tip_full(id_tip_w) = 1;
e_tip = e_tip_full(free);

% Precompute base LHS for implicit central difference
A_base = (1/dt^2)*Mff + Kff + (1/(2*dt))*Cff;

% Initial conditions (u^0, v^0 = 0), consistent a^0 with DNC
ndof = size(M,1);
U0  = zeros(ndof,1);
V0  = zeros(ndof,1);
chi0 = double(U0(id_tip_w) <= y_minus);
bDNC0 = chi0 * ( kappa*(U0(id_tip_w)-y_minus) - beta*V0(id_tip_w) ) * e_tip_full;
A0 = M \ ( F - C*V0 - K*U0 - bDNC0 );

% Start-up u^1
U1 = U0 + dt*V0 + 0.5*dt^2*A0;

% Enforce left BC strongly: w(0)=0, theta(0)=0
U0(1)=0; U0(2)=0;
U1(1)=0; U1(2)=0;

% Store tip
u_right = zeros(1,Nt+1);
u_right(1) = U0(id_tip_w);
u_right(2) = U1(id_tip_w);

% March in time
Unm1 = U0;  % u^{n-1}
Un   = U1;  % u^{n}
for n = 1:Nt-1
    % Active set from current state (matches FDM's use of u^n)
    chi_n = double(Un(id_tip_w) <= y_minus);

    % LHS: add positive (kappa + beta/dt) when in contact (rank-1 at tip)
    A = A_base + chi_n*(kappa + beta/dt)*(e_tip*e_tip.');

    % RHS: mass + viscous CD stencils + tip constant part
    R = Ff ...
        + (2/dt^2)*Mff*Un(free) - (1/dt^2)*Mff*Unm1(free) ...
        + (1/(2*dt))*Cff*Unm1(free);
    R = R + chi_n*( kappa*y_minus + (beta/dt)*Un(id_tip_w) )*e_tip;

    % Solve for u^{n+1} (free), assemble full, apply BC
    Ufp1 = A \ R;
    Up1 = zeros(ndof,1); Up1(free) = Ufp1;
    Up1(1)=0; Up1(2)=0;

    % record and advance
    u_right(n+2) = Up1(id_tip_w);
    Unm1 = Un; Un = Up1;
end
end

function [M, Kb, F, free, id_tip_w] = assembleHermiteCubic1D(J, h, f0)
% Cubic Hermite beam element (C^1), DOFs per node: [w, theta]
% Nodes: 0..J, elements: [x_{e-1}, x_e]
% Essential BCs: w(0)=0, theta(0)=0 handled by DOF elimination at solve time.

nn   = J+1;
ndof = 2*nn;

M  = zeros(ndof); 
Kb = zeros(ndof); 
F  = zeros(ndof,1);

% Local (length h): geometric bending Kb_e and consistent mass M_e
Kb_e = (1/h^3)*[ 12,     6*h,   -12,     6*h;
                 6*h,  4*h^2,   -6*h,  2*h^2;
                -12,    -6*h,    12,    -6*h;
                 6*h,  2*h^2,   -6*h,  4*h^2 ];
M_e  = (h/420)*[ 156,    22*h,    54,   -13*h;
                 22*h,  4*h^2,  13*h,   -3*h^2;
                 54,    13*h,   156,   -22*h;
                -13*h, -3*h^2, -22*h,    4*h^2 ];
% Consistent element load for f = const
F_e  = f0 * [h/2; h^2/12; h/2; -h^2/12];

for e = 1:J
    a = 2*e-1; edofs = [a, a+1, a+2, a+3];
    Kb(edofs,edofs) = Kb(edofs,edofs) + Kb_e;
    M(edofs,edofs)  = M(edofs,edofs)  + M_e;
    F(edofs)        = F(edofs)        + F_e;
end

% DOFs: [w0,th0, w1,th1, ..., wJ,thJ]
fixed = [1,2];              % left end clamped in value and slope
free  = setdiff(1:ndof, fixed);

id_tip_w = 2*nn - 1;        % right-end displacement DOF
end
