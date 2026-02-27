function beam_dnc1obsFEM_compare_beta_visc()
% FEM with C^1 Hermite cubics + implicit central difference in time,
% interior viscosity C = epsv * Kb, and unilateral DNC at the tip (one stop).
% Plots u(1,t) for multiple beta values on one figure, like the FDM sample.

%% Fixed parameters (match the sample)
L = 1.0; T = 10.0; dx = 0.01; dt = 0.01;
alpha2 = 1;                % use the same alpha^2 as your FDM sample
f0 = -0.2;                  % uniform load
y_minus = -0.02;            % obstacle level
kappa = 10;                  % contact stiffness (as in the FDM title)
betas = [1e-3, 1e-2, 1e-1, 1, 5];
epsv  = 0.001;               % interior viscosity (can set to 0 to match FDM)

% Mesh/time grids
J  = round(L/dx);      dx = L/J;
Nt = round(T/dt);      dt = T/Nt;
x  = linspace(0,L,J+1).';
t  = linspace(0,T,Nt+1);

% Build FE operators
[ M, Kb, F, free, id_tip_w ] = assembleHermiteCubic1D(J, dx, f0);
K = alpha2 * Kb;
C = epsv   * Kb;

% Reduced blocks
Mff = M(free,free);
Kff = K(free,free);
Cff = C(free,free);
Ff  = F(free);

% tip selector in reduced space
e_tip_full = zeros(size(M,1),1); e_tip_full(id_tip_w) = 1;
e_tip = e_tip_full(free);

% Precompute time-invariant part of LHS
A_base = (1/dt^2)*Mff + Kff + (1/(2*dt))*Cff;

% Storage for traces
tip_hist = zeros(numel(betas), Nt+1);

figure; hold on; grid on;

for ib = 1:numel(betas)
    beta_dnc = betas(ib);

    % Initial conditions (u^0, v^0, a^0)
    ndof = size(M,1);
    u0  = zeros(ndof,1);
    v0  = zeros(ndof,1);

    % Initial contact and a^0 from semi-discrete ODE
    chi0 = double(u0(id_tip_w) <= y_minus);
    bDNC0 = chi0 * ( kappa*(u0(id_tip_w)-y_minus) - beta_dnc*v0(id_tip_w) ) * e_tip_full; % NOTE sign + (LHS)
    a0 = M \ ( F - C*v0 - K*u0 - bDNC0 );

    % Build u^0, u^1 (central-diff startup)
    U0 = u0;
    U1 = u0 + dt*v0 + 0.5*dt^2*a0;

    % Enforce left BC strongly (w(0)=0, theta(0)=0)
    U0(1)=0; U0(2)=0;
    U1(1)=0; U1(2)=0;

    % Record tip
    tip_hist(ib,1) = U0(id_tip_w);
    tip_hist(ib,2) = U1(id_tip_w);

    % Time stepping: build u^{n+1} from u^n, u^{n-1}
    Unm1 = U0;   % u^{n-1}
    Un   = U1;   % u^{n}

    for n = 1:Nt-1
        % Active-set from current (previous) tip state, like the FDM sample
        chi_n = double(Un(id_tip_w) <= y_minus);

        % LHS (free DOFs): add tip contact stiffness (kappa + beta/dt)*chi_n
        A = A_base + chi_n * (kappa + beta_dnc/dt) * (e_tip*e_tip.');

        % RHS (free DOFs)
        % Central-difference mass & interior viscosity:
        % RHS = F^{n+1} + (2/dt^2)M U^n - (1/dt^2)M U^{n-1} + (1/(2dt)) C U^{n-1}
        R = Ff + (2/dt^2)*Mff*Un(free) - (1/dt^2)*Mff*Unm1(free) + (1/(2*dt))*Cff*Unm1(free);

        % Tip constant RHS from DNC (backward Euler at the boundary)
        % + chi_n * [ kappa*y_- + (beta/dt) * u_tip^n ] * e_tip
        R = R + chi_n * ( kappa*y_minus + (beta_dnc/dt)*Un(id_tip_w) ) * e_tip;

        % (No inhomog BCs at x=0 here; if present, subtract K_fc*known^{n+1}
        %  and (1/(2dt))C_fc*known^{n+1} and add (1/(2dt))C_fc*known^{n-1}.)

        % Solve and assemble full vector
        Ufp1 = A \ R;
        Up1 = zeros(ndof,1);
        Up1(free) = Ufp1;
        Up1(1)=0; Up1(2)=0;     % left BC

        % Write back & record
        tip_hist(ib,n+2) = Up1(id_tip_w);

        % Shift time levels
        Unm1 = Un;
        Un   = Up1;
    end

    % Plot this curve
    plot(t, tip_hist(ib,:), 'LineWidth',1.6, 'DisplayName', sprintf('\\beta = %g', beta_dnc));
end

% Obstacle line
yline(y_minus, '--', 'LineWidth',1.6, 'DisplayName', sprintf('y_- = %.3g', y_minus));

xlim([0 T]); xlabel('$Time$','Interpreter','latex'); ylabel('$u(1,t)$','Interpreter','latex');grid off;
title(sprintf('FEM Solution: Right Tip Displacement vs Time ($\\kappa=%.3g$, $\\varepsilon=%.3g$)', kappa, epsv), ...
      'Interpreter','latex');

legend('Location','southeast');

end

% -------------------------------------------------------------------------
function [M, Kb, F, free, id_tip_w] = assembleHermiteCubic1D(J, h, f0)
% Cubic Hermite beam element (C^1), DOFs per node: [w, theta]
% Nodes: 0..J, elements: [x_{e-1}, x_e]
% Essential BCs: w(0)=0, theta(0)=0 enforced by DOF elimination.

nn   = J+1;
ndof = 2*nn;

M  = zeros(ndof); 
Kb = zeros(ndof); 
F  = zeros(ndof,1);

% Local matrices (length h), geometric bending Kb_e and mass M_e (consistent)
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

% Essential BCs at x=0: remove DOFs 1 (w0) and 2 (theta0)
fixed = [1,2];
free  = setdiff(1:ndof, fixed);

% Tip displacement DOF index (right-end w)
id_tip_w = 2*nn - 1;

end
