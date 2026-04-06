function [x,t,u,metrics] = beam_contact_implicit_FD(f0,kappa,beta,params)
% Implicit second-order FD for Euler–Bernoulli beam with unilateral contact
% Inputs:
%   f0, kappa, beta : scalars (forcing, contact stiffness, contact damping)
%   params (optional struct): overrides for {L,T,dx,dt,alpha2,y_minus}
% Outputs:
%   x,t,u : grids and displacement
%   metrics: struct with tip, contact, etc., useful for summaries

if nargin < 4, params = struct(); end

% ---- Defaults; allow overrides via params ----
L       = getfielddef(params,'L',1.0);
T       = getfielddef(params,'T',5.0);
dx      = getfielddef(params,'dx',0.01);
dt      = getfielddef(params,'dt',0.01);
alpha2  = getfielddef(params,'alpha2',1);
y_minus = getfielddef(params,'y_minus',-0.02);

% Left boundary/ICs
phi = @(tt) 0.0;
u0  = @(xx) 0.0;
v0  = @(xx) 0.0;

M  = round(L/dx);
N  = round(T/dt);
gamma = alpha2*(dt^2)/(dx^4);

% Grids & storage
x = linspace(0,L,M+1).';
t = linspace(0,T,N+1);
u = zeros(M+1,N+1);

% ICs
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = phi(t);  % enforce at x=0

% ---- Time stepping ----
for j = 1:(N-1)
    in_contact = (u(M+1,j+1) <= y_minus);

    A = spalloc(M,M,5*M);
    b = zeros(M,1);

    % i=1 (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M>=2, A(1,2) = -4*gamma; end
    if M>=3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2,j+1) - u(2,j) + 4*gamma*phi(t(j+1));

    % interior i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r,i-2) = A(r,i-2) + gamma;
        else
            b(r) = b(r) - gamma*phi(t(j+1));
        end
        A(r,i-1) = A(r,i-1) - 4*gamma;
        A(r,i  ) = A(r,i  ) + (1 + 6*gamma);
        A(r,i+1) = A(r,i+1) - 4*gamma;
        A(r,i+2) = A(r,i+2) + gamma;
        b(r) = b(r) + dt^2*f0 + 2*u(i+1,j+1) - u(i+1,j);
    end

    % i=M-1
    if M>=2
        A(M-1,M-3) = gamma;
        A(M-1,M-2) = -4*gamma;
        A(M-1,M-1) = 1 + 5*gamma;
        A(M-1,M  ) = -2*gamma;
        b(M-1) = dt^2*f0 + 2*u(M    ,j+1) - u(M    ,j);
    end

    % i=M : no-contact vs contact
    if ~in_contact
        A(M,M-2) =  2*gamma;
        A(M,M-1) = -4*gamma;
        A(M,M  ) =  1 + 2*gamma;
        b(M)     =  dt^2*f0 + 2*u(M+1,j+1) - u(M+1,j);
    else
        c        = (dx^3/alpha2)*(kappa + beta/dt);
        A(M,M-2) =  2*gamma;
        A(M,M-1) = -4*gamma;
        A(M,M  ) =  1 + 2*gamma*(c + 1);
        b(M)     =  dt^2*f0 ...
                  + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1,j+1) ...
                  - u(M+1,j) ...
                  + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    % Solve & write back
    w = A \ b;
    u(1,   j+2) = phi(t(j+2));
    u(2:end,j+2) = w;
end

% ---- Basic figure (optional) ----
% figure; plot(t,u(end,:),'LineWidth',1.2); grid on;
% xlabel('time'); ylabel('u(L,t)');

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

% --- small helpers ---
function v = getfielddef(s,field,default)
if isfield(s,field), v = s.(field); else, v = default; end
end
function y = iff(cond,a,b)
if cond, y=a; else, y=b; end
end
