% FEM for Euler–Bernoulli beam with TWO-OBSTACLE DNC tip law (Hermite C1 + Newmark)
% Time-varying body force q(t) = F0*sin(omega*t) to drive bouncing between obstacles.
clear; clc;

% ---- Physical/model parameters ----
L = 1.0;
alpha2 = 1.296;            % coefficient in alpha^2 * u_xxxx
y_minus = -0.2;           % bottom obstacle level
y_plus  =  0.2;           % top obstacle level
kappa   = 10;             % DNC stiffness
betas   = [1e-2, 1, 5];  % sweep of DNC damping values
phi     = @(t) 0.0;        % left-end displacement (clamp at zero)

% ---- Time-varying load (uniform in x) ----
F0    = 5.0;               % amplitude (tune up/down to induce contact)
omega = 10.0;              % angular frequency (rad/s)
qfun  = @(t) F0*sin(omega*t);

% ---- Time integration ----
T  = 20.0;
dx = 0.01;  dt = 0.01;
J  = round(L/dx);          % number of elements
dx = L/J;                  % snap to integer mesh
Nt = round(T/dt);          % number of steps
dt = T/Nt;                 % snap to grid
tgrid = linspace(0,T,Nt+1);

% Newmark parameters (average acceleration)
gNM = 1/2; bNM = 1/4;

% ---- Mesh and DOFs ----
nn   = J+1;                % nodes
ndof = 2*nn;               % [w_i, theta_i] per node
fixed = [1,2];             % clamp: w(0)=phi(t), theta(0)=0
free  = setdiff(1:ndof, fixed);
id_tip_w  = 2*nn-1;        % right-end displacement DOF
id_tip_th = 2*nn;          % (unused, slope DOF)

% ---- Assemble M, K ----
M = zeros(ndof);  K = zeros(ndof);
h = dx;
Ke = (alpha2/h^3)*[ 12,   6*h, -12,   6*h;
                     6*h, 4*h^2, -6*h, 2*h^2;
                    -12,  -6*h,  12,  -6*h;
                     6*h, 2*h^2, -6*h, 4*h^2];
Me = (h/420)*[ 156,   22*h,  54,  -13*h;
               22*h,  4*h^2, 13*h, -3*h^2;
               54,    13*h, 156,  -22*h;
              -13*h, -3*h^2, -22*h,  4*h^2];

for e = 1:J
    a = 2*e-1; edofs = [a, a+1, a+2, a+3];
    K(edofs,edofs) = K(edofs,edofs) + Ke;
    M(edofs,edofs) = M(edofs,edofs) + Me;
end

% ---- Consistent load for unit uniform q=1 (time scaling done in loop) ----
% Element consistent load for q=1:
fe_unit = [h/2; h^2/12; h/2; -h^2/12];
F_unit  = zeros(ndof,1);
for e = 1:J
    a = 2*e-1; edofs = [a, a+1, a+2, a+3];
    F_unit(edofs) = F_unit(edofs) + fe_unit;
end

% Reduced matrices/vectors
Kff = K(free,free);
Mff = M(free,free);
F_unit_f = F_unit(free);

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
        tnp1 = tgrid(n+1);

        % left Dirichlet at n+1 (only displacement DOF)
        w0_np1 = phi(tnp1);  th0_np1 = 0;

        % Newmark predictors on free DOFs
        u_pred_f = uf + dt*vf + dt^2*(0.5 - bNM)*af;
        v_pred_f = vf + (1 - gNM)*dt*af;

        % ---- TWO-OBSTACLE DNC: smooth contact indicators + short iteration ----
        eta = 1e-8;        % regularization for tanh switch
        maxIter = 4;

        % Start indicators from current tip value
        u_tip_guess  = (e_tip.'*uf);
        chi_minus = 0.5*(1 - tanh((u_tip_guess - y_minus)/eta));    % ~1 if u<=y_-
        chi_plus  = 0.5*(1 + tanh((u_tip_guess - y_plus )/eta));    % ~1 if u>=y_+
        chi_minus = max(0,min(1,chi_minus));
        chi_plus  = max(0,min(1,chi_plus));

        % time-varying distributed load at t_{n+1}
        qn1 = qfun(tnp1);
        Ff  = qn1 * F_unit_f;   % scale unit consistent vector

        u_np1_f = u_pred_f;  % placeholder

        for it = 1:maxIter
            s = chi_minus + chi_plus;  % total activation (0 or ~1)

            % Tangent on tip DOF from DNC law (moved to LHS):
            %  + kappa*u + beta*(gNM/(bNM*dt))*u   (active only when in contact)
            c_tan = s * ( kappa + beta_dnc * gNM/(bNM*dt) );

            % Effective matrix and RHS on free DOFs
            K_eff = Kff + (1/(bNM*dt^2))*Mff + c_tan*(e_tip*e_tip.');
            R = Ff + (1/(bNM*dt^2))*Mff*u_pred_f;

            % effect of known left BCs at n+1:
            known_np1 = [w0_np1; th0_np1];
            R = R - K(free,fixed)*known_np1;

            % Boundary-force constant part on RHS:
            % κ(χ_- y_- + χ_+ y_+)  +  s*β( (gNM/(bNM*dt))*u_pred_tip - v_pred_tip )
            u_pred_tip = e_tip.'*u_pred_f;
            v_pred_tip = e_tip.'*v_pred_f;
            g_const = kappa*(chi_minus*y_minus + chi_plus*y_plus) ...
                    + s * beta_dnc * ( (gNM/(bNM*dt))*u_pred_tip - v_pred_tip );
            R = R + g_const * e_tip;

            % Solve for u_{n+1} on free DOFs
            u_np1_f = K_eff \ R;

            % Update indicators with NEW tip and recheck
            u_tip_new = e_tip.'*u_np1_f;
            chi_minus_new = 0.5*(1 - tanh((u_tip_new - y_minus)/eta));
            chi_plus_new  = 0.5*(1 + tanh((u_tip_new - y_plus )/eta));
            chi_minus_new = max(0,min(1,chi_minus_new));
            chi_plus_new  = max(0,min(1,chi_plus_new));

            if abs(chi_minus_new-chi_minus)+abs(chi_plus_new-chi_plus) < 1e-6
                break;
            else
                chi_minus = chi_minus_new; chi_plus = chi_plus_new;
            end
        end

        % Newmark corrector on free DOFs
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
figure; hold on; box on; grid on;
clr = lines(numel(betas));
for ib = 1:numel(betas)
    plot(tgrid, tip_hist(ib,:), 'LineWidth',1.0, 'Color',clr(ib,:), ...
        'DisplayName', sprintf('\\beta = %g', betas(ib)));
end
yline(y_minus, '--k', 'HandleVisibility','off');
yline(y_plus,  '--k', 'HandleVisibility','off');
xlabel('time'); ylabel('u(L,t)');
title(sprintf('Tip displacement, \\kappa=%.3g, two-obstacle DNC, q(t)=%.1fsin(%gt)', ...
      kappa, F0, omega));
legend('Location','best');
