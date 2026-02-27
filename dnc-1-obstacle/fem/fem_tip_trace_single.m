function [t, u_right] = fem_tip_trace_single(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv)
% FEM + implicit central difference in time + unilateral DNC at the tip.
% Returns time vector and right-end displacement.

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

% Precompute base LHS for implicit CD
A_base = (1/dt^2)*Mff + Kff + (1/(2*dt))*Cff;

% Initial conditions
ndof = size(M,1);
U0  = zeros(ndof,1);
V0  = zeros(ndof,1);
chi0 = double(U0(id_tip_w) <= y_minus);
bDNC0 = chi0 * ( kappa*(U0(id_tip_w)-y_minus) - beta*V0(id_tip_w) ) * e_tip_full;
A0 = M \ ( F - C*V0 - K*U0 - bDNC0 );

% Start-up
U1 = U0 + dt*V0 + 0.5*dt^2*A0;
U0(1)=0; U0(2)=0;  U1(1)=0; U1(2)=0;

u_right = zeros(1,Nt+1);
u_right(1) = U0(id_tip_w);
u_right(2) = U1(id_tip_w);

Unm1 = U0;  % u^{n-1}
Un   = U1;  % u^{n}

for n = 1:Nt-1
    chi_n = double(Un(id_tip_w) <= y_minus);

    % LHS
    A = A_base + chi_n*(kappa + beta/dt)*(e_tip*e_tip.');

    % RHS
    R = Ff + (2/dt^2)*Mff*Un(free) - (1/dt^2)*Mff*Unm1(free) + (1/(2*dt))*Cff*Unm1(free);
    R = R + chi_n*( kappa*y_minus + (beta/dt)*Un(id_tip_w) )*e_tip;

    % Solve and assemble
    Ufp1 = A \ R;
    Up1 = zeros(ndof,1);
    Up1(free) = Ufp1;
    Up1(1)=0; Up1(2)=0;

    % record + advance
    u_right(n+2) = Up1(id_tip_w);
    Unm1 = Un; Un = Up1;
end
end
