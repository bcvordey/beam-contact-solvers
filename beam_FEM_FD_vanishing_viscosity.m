function beam_FEM_FD_vanishing_viscosity()
% Compare FEM (with viscosity epsv > 0) against FD (epsv = 0 reference).
% Computes L2 and Linf errors at the tip, and convergence rates.

%% Parameters
epsvals = [0.05,0.025,0.0125,0.00625,0.003125];   % viscosities for FEM
L = 1.0; T = 10.0; dx = 0.01; dt = 0.01;
alpha2 = 1; f0 = -0.2; y_minus = -0.02; kappa = 1; beta = 0.1;

%% Run FD reference (no viscosity)
[ tref, uref ] = simulate_FD_tip(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);

%% Loop over epsv for FEM
numEps = numel(epsvals);
L2_err = zeros(numEps,1);
Linf_err = zeros(numEps,1);

for k = 1:numEps
    epsv = epsvals(k);
    [ tFEM, uFEM ] = simulate_FEM_tip(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv);

    % Align grids (should be identical if dx, dt, T match)
    if ~isequal(tref,tFEM)
        error('Time grids differ between FEM and FD solvers!');
    end

    % Error at tip trace
    e = uFEM - uref;

    % L2 and Linf norms
    L2_err(k)   = sqrt(sum(e.^2)*dt);
    Linf_err(k) = max(abs(e));
end

%% Compute rates
pL2   = nan(numEps,1);
pLinf = nan(numEps,1);
for k = 2:numEps
    pL2(k)   = log(L2_err(k-1)/L2_err(k)) / log(epsvals(k-1)/epsvals(k));
    pLinf(k) = log(Linf_err(k-1)/Linf_err(k)) / log(epsvals(k-1)/epsvals(k));
end

%% Print table
fprintf('Convergence study FEM(eps) -> FD(0)\n');
fprintf(' epsv       L2-error      Linf-error      p(L2)      p(Linf)\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:numEps
    fprintf('%8.4f   %12.4e   %12.4e   %8.3f   %8.3f\n', ...
        epsvals(k), L2_err(k), Linf_err(k), pL2(k), pLinf(k));
end

end

% -------------------------------------------------------------------------
function [t,u_right] = simulate_FD_tip(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Wrapper for FD scheme
[t,u_right] = simulate_right_end(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);
end

% -------------------------------------------------------------------------
function [t,u_tip] = simulate_FEM_tip(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv)
% Wrapper for FEM scheme with viscosity epsv
betas = beta;  % use single beta
% Use your FEM solver but pass epsv in
[t,u_tip] = run_FEM_one_beta(beta,kappa,y_minus,L,T,dx,dt,alpha2,f0,epsv);
end

% -------------------------------------------------------------------------
function [t, u_tip] = run_FEM_one_beta(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0, epsv)
% Adapted from your FEM code, returns time vector and tip trace

J  = round(L/dx); dx = L/J;
Nt = round(T/dt); dt = T/Nt;
t  = linspace(0,T,Nt+1);

% Build FE operators
[M,Kb,F,free,id_tip_w] = assembleHermiteCubic1D(J, dx, f0);
K = alpha2 * Kb;
C = epsv   * Kb;

Mff = M(free,free);
Kff = K(free,free);
Cff = C(free,free);
Ff  = F(free);

e_tip_full = zeros(size(M,1),1); e_tip_full(id_tip_w) = 1;
e_tip = e_tip_full(free);

A_base = (1/dt^2)*Mff + Kff + (1/(2*dt))*Cff;

% Initial conditions
ndof = size(M,1);
u0 = zeros(ndof,1);
v0 = zeros(ndof,1);
chi0 = double(u0(id_tip_w) <= y_minus);
bDNC0 = chi0 * ( kappa*(u0(id_tip_w)-y_minus) - beta*v0(id_tip_w) ) * e_tip_full;
a0 = M \ ( F - C*v0 - K*u0 - bDNC0 );

U0 = u0;
U1 = u0 + dt*v0 + 0.5*dt^2*a0;
U0(1)=0; U0(2)=0; U1(1)=0; U1(2)=0;

u_tip = zeros(1,Nt+1);
u_tip(1) = U0(id_tip_w);
u_tip(2) = U1(id_tip_w);

Unm1 = U0; Un = U1;

for n=1:Nt-1
    chi_n = double(Un(id_tip_w) <= y_minus);
    A = A_base + chi_n*(kappa + beta/dt)*(e_tip*e_tip.');
    R = Ff + (2/dt^2)*Mff*Un(free) - (1/dt^2)*Mff*Unm1(free) + (1/(2*dt))*Cff*Unm1(free);
    R = R + chi_n*(kappa*y_minus + (beta/dt)*Un(id_tip_w))*e_tip;
    Ufp1 = A\R;
    Up1 = zeros(ndof,1);
    Up1(free)=Ufp1;
    Up1(1)=0; Up1(2)=0;
    u_tip(n+2)=Up1(id_tip_w);
    Unm1=Un; Un=Up1;
end
end


% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Runs one FD simulation for a given beta; returns time vector and right-end trace.

M  = round(L/dx);              
N  = round(T/dt);              
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

phi = @(tt) 0.0;               
u0  = @(xx) 0.0;               
v0  = @(xx) 0.0;               

gamma = alpha2*(dt^2)/(dx^4);

u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0;                    

for j = 1:(N-1)
    in_contact = (u(M+1, j+1) <= y_minus);
    A = spalloc(M,M,5*M);
    b = zeros(M,1);

    % Interior equations (same as your FD code) ...
    % [ Paste from your original simulate_right_end definition ]
    % ---------------------------------------------------------

    % i=1 row
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2,j+1) - u(2,j);

    % interior rows
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r,i-2) = A(r,i-2) + gamma;
        end
        A(r,i-1) = A(r,i-1) - 4*gamma;
        A(r,i)   = A(r,i)   + (1 + 6*gamma);
        A(r,i+1) = A(r,i+1) - 4*gamma;
        A(r,i+2) = A(r,i+2) + gamma;
        b(r) = dt^2*f0 + 2*u(i+1,j+1) - u(i+1,j);
    end

    % i=M-1 row
    A(M-1,M-3)=gamma;
    A(M-1,M-2)=-4*gamma;
    A(M-1,M-1)=1+5*gamma;
    A(M-1,M)  =-2*gamma;
    b(M-1)=dt^2*f0 + 2*u(M,j+1)-u(M,j);

    % i=M row: contact or free
    if ~in_contact
        A(M,M-2)=2*gamma; A(M,M-1)=-4*gamma; A(M,M)=1+2*gamma;
        b(M)=dt^2*f0 + 2*u(M+1,j+1)-u(M+1,j);
    else
        c = (dx^3/alpha2)*(kappa+beta/dt);
        A(M,M-2)=2*gamma; A(M,M-1)=-4*gamma; A(M,M)=1+2*gamma*(c+1);
        b(M)=dt^2*f0 + 2*(1+gamma*(beta*dx^3)/(alpha2*dt))*u(M+1,j+1)...
             -u(M+1,j) + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    w=A\b;
    u(1,j+2)=0;
    u(2:end,j+2)=w;
end

u_right = u(end,:);
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