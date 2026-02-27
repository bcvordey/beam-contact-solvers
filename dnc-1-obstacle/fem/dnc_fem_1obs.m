% FEM for Euler–Bernoulli beam with DNC tip law (Hermite C1 + Newmark)
% Produces tip-displacement curves for kappa=0.1 and beta in {1e-3,0.1,1,5}
clear; clc;

% ---- Physical/model parameters (match your simulation chapter) ----
L = 1.0;
alpha2 = 1;         % coefficient in alpha^2 * u_xxxx
fconst = -1;          % distributed load (constant)
y_minus = -0.02;        % obstacle level
kappa = 10;            % obstacle stiffness
betas = [1e-3, 1e-2, 1e-1, 1, 5];  % DNC damping values to loop over
phi = @(t) 0.0;         % left-end displacement
% time integration
T = 10.0;
dx = 0.01;  dt = 0.01;  % as requested
J = round(L/dx);        % number of elements
dx = L/J;               % snap to integer mesh
Nt = round(T/dt);       % number of steps
dt = T/Nt;              % snap to grid
tgrid = linspace(0,T,Nt+1);

% Newmark parameters (average acceleration)
gNM = 1/2; bNM = 1/4;

% ---- Mesh and DOFs ----
nn = J+1;         % nodes
ndof = 2*nn;      % [w_i, theta_i] per node
% left BC: w(0)=phi(t), theta(0)=0  -> eliminate DOFs 1 and 2
fixed = [1,2];    % fixed DOFs
free = setdiff(1:ndof, fixed);
% index of tip displacement DOF and slope DOF
id_tip_w = 2*nn-1; 
id_tip_th = 2*nn;

% ---- Assemble M, K, F (consistent) ----
M = zeros(ndof);  K = zeros(ndof);  F = zeros(ndof,1);
h = dx;
Ke = (alpha2/h^3)*[ 12,   6*h, -12,   6*h;
                     6*h, 4*h^2, -6*h, 2*h^2;
                    -12,  -6*h,  12,  -6*h;
                     6*h, 2*h^2, -6*h, 4*h^2];
Me = (h/420)*[ 156,   22*h,  54,  -13*h;
               22*h,  4*h^2, 13*h, -3*h^2;
               54,    13*h, 156,  -22*h;
              -13*h, -3*h^2, -22*h, 4*h^2];
% consistent element load for f = const
Fe = fconst * [h/2; h^2/12; h/2; -h^2/12];

for e = 1:J
    % element dof map [w_i, th_i, w_{i+1}, th_{i+1}]
    a = 2*e-1; edofs = [a, a+1, a+2, a+3];
    K(edofs,edofs) = K(edofs,edofs) + Ke;
    M(edofs,edofs) = M(edofs,edofs) + Me;
    F(edofs)       = F(edofs)       + Fe;
end

% Apply inhomogeneous w(0)=phi(t) via shifting loads at each step
% Build reduced matrices
Kff = K(free,free);
Mff = M(free,free);
Ff  = F(free);

% Unit vector for tip displacement in reduced system
e_tip_full = zeros(ndof,1); e_tip_full(id_tip_w)=1;
e_tip = e_tip_full(free);

% storage for curves
tip_hist = zeros(numel(betas), Nt+1);

% ---- Time loop for each DNC beta ----
for ib = 1:numel(betas)
    beta_dnc = betas(ib);        % DNC damping
    % init states (full then reduced)
    u  = zeros(ndof,1);  v = zeros(ndof,1);  a = zeros(ndof,1);
    % enforce initial left BC
    u(1) = phi(0); u(2) = 0; % slope = 0
    uf = u(free); vf = v(free); af = a(free);

    % record tip
    tip_hist(ib,1) = u(id_tip_w);

    for n = 1:Nt
        tn   = tgrid(n);
        tnp1 = tgrid(n+1);

        % left Dirichlet at n+1 (only displacement DOF)
        w0_np1 = phi(tnp1);  th0_np1 = 0;

        % build effective RHS from inhomogeneous BC shift:
        % split unknown as (free DOFs) + known fixed values
        % here known values only at DOFs 1 (w0) and 2 (th0)
        uc_known = zeros(ndof,1); uc_known(1)=w0_np1; uc_known(2)=th0_np1;
        % predictor on free DOFs
        u_pred_f = uf + dt*vf + dt^2*(0.5 - bNM)*af;
        v_pred_f = vf + (1 - gNM)*dt*af;

        % ---- Contact iteration (fully implicit law) ----
        % start with contact status from previous tip value
        u_tip_guess = u(id_tip_w);
        % smooth indicator (tiny eta to avoid chattering in algebra)
        eta = 1e-8;
        chi = 0.5*(1 - tanh((u_tip_guess - y_minus)/eta)); % ~1 if u <= y_-
        chi = max(0,min(1,chi));

        maxIter = 4;
        for it = 1:maxIter
            % tangent coefficient induced by boundary law at the tip
            c_tan = chi*( kappa + beta_dnc * gNM/(bNM*dt) );

            % Effective matrix and RHS on free DOFs
            K_eff = Kff + (1/(bNM*dt^2))*Mff + c_tan*(e_tip*e_tip.');
            % body load at n+1
            R = Ff + (1/(bNM*dt^2))*Mff*u_pred_f;

            % add effect of known left BCs at n+1:
            % subtract K(:,fixed)*u_known from RHS
            known_np1 = [w0_np1; th0_np1];
            R = R - K(free,fixed)*known_np1;

            % boundary-force constant part
            u_pred_tip = e_tip.'*u_pred_f + 0;           % tip part from free predictors
            v_pred_tip = e_tip.'*v_pred_f + 0;
            g_const = chi*( kappa*y_minus - beta_dnc*v_pred_tip + beta_dnc*(gNM/(bNM*dt))*u_pred_tip );
            R = R + g_const * e_tip;

            % Solve
            u_np1_f = K_eff \ R;

            % compute new contact status
            u_full_np1 = zeros(ndof,1);
            u_full_np1(free) = u_np1_f; u_full_np1(1)=w0_np1; u_full_np1(2)=th0_np1;
            u_tip_new = u_full_np1(id_tip_w);
            chi_new = 0.5*(1 - tanh((u_tip_new - y_minus)/eta));
            chi_new = max(0,min(1,chi_new));

            if abs(chi_new - chi) < 1e-6, break; else, chi = chi_new; end
        end

        % update accelerations/velocities on free DOFs
        a_np1_f = (u_np1_f - u_pred_f)/(bNM*dt^2);
        v_np1_f = v_pred_f + gNM*dt*a_np1_f;

        % write back full vectors
        u = zeros(ndof,1); v = zeros(ndof,1); a = zeros(ndof,1);
        u(free)=u_np1_f; v(free)=v_np1_f; a(free)=a_np1_f;
        u(1)=w0_np1; u(2)=th0_np1;  % fixed DOFs

        % store tip displacement
        tip_hist(ib,n+1) = u(id_tip_w);

        % shift to next step
        uf = u(free);  vf = v(free);  af = a(free);
    end
end

% ---- Plot: tip displacement vs time for the requested betas ----
figure; hold on; box on;
plot(tgrid, tip_hist(1,:), 'r', 'LineWidth',1.5);
plot(tgrid, tip_hist(2,:), 'b', 'LineWidth',1.5);
plot(tgrid, tip_hist(3,:), 'g', 'LineWidth',1.5);
plot(tgrid, tip_hist(4,:), 'y', 'LineWidth',1.5);
yline(y_minus, '--k'); 
xlabel('time'); ylabel('u(L,t)');
title('Tip displacement, \kappa=0.1, varying \beta (FEM+Newmark+DNC)');
legend('\beta=10^{-2}','\beta=10^{-1}','\beta=1','\beta=5','obstacle','Location','best');
