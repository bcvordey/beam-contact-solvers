function [x,t,u,metrics] = beam_contact_implicit_FD(f0,kappa,beta,params)
% Implicit second-order FD for an Euler–Bernoulli beam with unilateral contact
% Unknowns at each step: w = [u_1^{j+1}, ..., u_M^{j+1}]^T
% Geometry: nodes i=0..M (so u has size (M+1) x (N+1)); u(1, :) is i=0 (left end)
%
% Inputs
%   f0, kappa, beta : scalars (constant external forcing, obstacle stiffness, damping)
%   params (optional struct) may contain:
%       L, T, dx, dt, alpha2, y_minus
%       phi_t   : @(t) time-dependent left boundary displacement u(0,t)=phi(t)   [default 0]
%       u0fun   : @(x) initial displacement u(x,0)                               [default 0]
%       v0fun   : @(x) initial velocity    u_t(x,0)                              [default 0]
%
% Outputs
%   x, t, u   : spatial grid, time grid, displacement field
%   metrics   : struct with tip trace and contact metrics

if nargin < 4, params = struct(); end

% ---- Defaults; allow overrides via params ----
L       = getfielddef(params,'L',1.0);
T       = getfielddef(params,'T',5.0);
dx      = getfielddef(params,'dx',0.01);
dt      = getfielddef(params,'dt',0.01);
alpha2  = getfielddef(params,'alpha2',1.296);
y_minus = getfielddef(params,'y_minus',-0.02);

% Left boundary in time: u(0,t) = phi(t)
if isfield(params,'phi_t') && ~isempty(params.phi_t)
    phi = params.phi_t;                   % @(t)
else
    phi = @(tt) 0.0;
end

% Initial conditions in space (optional)
if isfield(params,'u0fun') && ~isempty(params.u0fun)
    u0 = params.u0fun;                    % @(x)
else
    u0 = @(xx) 0.0;
end
if isfield(params,'v0fun') && ~isempty(params.v0fun)
    v0 = params.v0fun;                    % @(x)
else
    v0 = @(xx) 0.0;
end

% ---- Derived sizes/parameters ----
M  = round(L/dx);              % spatial subintervals -> nodes are i=0..M
N  = round(T/dt);              % time steps           -> times are j=0..N
gamma = alpha2*(dt^2)/(dx^4);  % scheme parameter

% ---- Grids & storage (i=0..M, j=0..N -> MATLAB indices +1) ----
x = linspace(0,L,M+1).';       % column vector
t = linspace(0,T,N+1);         % row vector
u = zeros(M+1, N+1);           % u(i+1, j+1) ~ u_i^j

% ---- Initial conditions ----
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);

% Enforce Dirichlet at x=0 for all times (vectorized)
u(1,:) = arrayfun(phi, t);

% ---- Time stepping: build & solve A w = b at each j ----
for j = 1:(N-1)           % have u^0 and u^1; compute up to u^N

    % Contact decision at right end using u_M^j (index M+1, j+1)
    in_contact = (u(M+1, j+1) <= y_minus);

    % Assemble linear system for unknowns [u_1^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % ------- Row for i=1 (ghost u_{-1}^{j+1}=u_1^{j+1} enforces u_x(0,t)=0) -------
    % (1 + 7γ) u_1 - 4γ u_2 + γ u_3 = dt^2 f + 2 u_1^j - u_1^{j-1} + 4γ φ(t_{j+1})
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));

    % ------- Interior: i = 2 .. M-2 (five-point stencil) -------
    % γ u_{i-2} - 4γ u_{i-1} + (1+6γ) u_i - 4γ u_{i+1} + γ u_{i+2}
    %   = dt^2 f + 2u_i^j - u_i^{j-1}
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            % touches u_0^{j+1} = phi(t_{j+1}) (known)
            b(r) = b(r) - gamma*phi(t(j+1));
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % ------- i = M-1 row (uses u_{M+1} eliminated) -------
    % γ u_{M-3} - 4γ u_{M-2} + (1+5γ) u_{M-1} - 2γ u_M
    %   = dt^2 f + 2 u_{M-1}^j - u_{M-1}^{j-1}
    if M >= 2
        A(M-1, M-3) = gamma;
        A(M-1, M-2) = -4*gamma;
        A(M-1, M-1) = 1 + 5*gamma;
        A(M-1, M  ) = -2*gamma;
        b(M-1) = dt^2*f0 + 2*u(M, j+1) - u(M, j);
    end

    % ------- i = M row (right boundary): contact vs no-contact -------
    if ~in_contact
        % No contact:
        % 2γ u_{M-2} - 4γ u_{M-1} + (1+2γ) u_M
        %   = dt^2 f + 2 u_M^j - u_M^{j-1}
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)     =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        % Contact (penalty with linear damping at the tip)
        % effective coefficient c = (dx^3/α^2) * (κ + β/dt)
        c        = (dx^3/alpha2)*(kappa + beta/dt);
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);
        b(M)     =  dt^2*f0 ...
                  + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
                  - u(M+1, j) ...
                  + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    % Solve for w and write back
    w = A \ b;
    u(1,   j+2) = phi(t(j+2));  % enforce left boundary at next time
    u(2:end,j+2) = w;
end

% ---- Metrics for summaries ----
Tol = 1e-12;
tip = u(end,:);
contact = tip <= y_minus + Tol;
penetration = max(0, y_minus - tip);
first_idx = find(contact,1,'first');

metrics.t_first         = iff(isempty(first_idx), NaN, t(first_idx));
metrics.num_impacts     = sum(diff([0 contact 0])==1);
metrics.dmax            = max(penetration);
metrics.contact_frac    = mean(contact);
metrics.tip             = tip;
metrics.contact         = contact;
metrics.y_minus         = y_minus;
metrics.kappa           = kappa;
metrics.beta            = beta;
metrics.f0              = f0;

end

% ---------- small helpers ----------
function v = getfielddef(s,field,default)
if isfield(s,field) && ~isempty(s.(field)), v = s.(field); else, v = default; end
end

function y = iff(cond,a,b)
if cond, y=a; else, y=b; end
end
