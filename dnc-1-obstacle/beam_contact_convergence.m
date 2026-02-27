function beam_contact_convergence()
% Convergence check for the right-end trace u(1,t) with the implicit FD+DNC scheme.
% Compares successive meshes with h_{l+1} = h_l/2 and reports L2/Linf-in-time errors.
%
% Problem data (as in the prompt)
L = 1.0;  T = 5.0;  alpha2 = 1.296;  f0 = -1.5;  y_minus = -0.02;
kappa = 1e-2;       % tune if desired
beta  = 1e-2;       % pick a representative beta for the convergence test

% Mesh sequence (h = Δx = Δt), from coarse to fine
hs = [0.04, 0.02, 0.01, 0.005];

% Storage
U = cell(numel(hs),1);
Tgrid = cell(numel(hs),1);

% Run simulations
for k = 1:numel(hs)
    h = hs(k);
    dx = h; dt = h;
    [t, u_right] = run_beam(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0);
    Tgrid{k} = t(:).';         % row
    U{k}     = u_right(:).';   % row
end

% Compute errors between successive levels (fine restricted to coarse times)
E2   = nan(1,numel(hs)-1);     % L2(0,T) in time
Einf = nan(1,numel(hs)-1);     % Linf(0,T) in time
for k = 1:numel(hs)-1
    tC = Tgrid{k};     uC = U{k};
    tF = Tgrid{k+1};   uF = U{k+1};
    % since dt halves exactly, coarse times are every 2nd fine step
    ratio = round((hs(k))/hs(k+1));  % should be 2
    idx = 1:ratio:numel(tF);
    % defensive: align lengths in case of rounding
    m = min(numel(tC), numel(idx));
    diff = uC(1:m) - uF(idx(1:m));
    tuse = tC(1:m);
    % L2 and Linf in time (continuous norms via trapezoid and max)
    E2(k)   = sqrt(trapz(tuse, diff.^2));
    Einf(k) = max(abs(diff));
end

% Observed orders (between successive pairs)
p2   = nan(1,numel(E2)-1);
pinf = nan(1,numel(Einf)-1);
for k = 1:numel(E2)-1
    p2(k)   = log(E2(k)/E2(k+1))   / log(hs(k)/hs(k+1));
    pinf(k) = log(Einf(k)/Einf(k+1)) / log(hs(k)/hs(k+1));
end

% Pretty print table
fprintf('\nConvergence of right-end trace u(1,t)  (beta = %.3g, kappa = %.3g)\n', beta, kappa);
fprintf('----------------------------------------------------------------------------- \n');
fprintf('%10s %12s %14s %14s %10s %10s\n','h=dx=dt','N=T/h','L2 error','Linf error','p(L2)','p(Linf)');
fprintf('----------------------------------------------------------------------------- \n');
for k = 1:numel(hs)
    h = hs(k);
    N = round(T/h);
    if k < numel(hs)
        e2 = E2(k);  einf = Einf(k);
        if k==1
            fprintf('%10.4f %12d %14.4e %14.4e %10s %10s\n', h, N, e2, einf, '-', '-');
        else
            fprintf('%10.4f %12d %14.4e %14.4e %10.3f %10.3f\n', h, N, e2, einf, p2(k-1), pinf(k-1));
        end
    else
        % finest level has no error vs finer mesh
        fprintf('%10.4f %12d %14s %14s %10s %10s\n', h, N, '-', '-', '-', '-');
    end
end
fprintf('----------------------------------------------------------------------------- \n');

% Log–log plot
figure; hold on;
loglog(hs(1:end-1), E2,   'o-','LineWidth',1.6,'DisplayName','L^2_t error');
loglog(hs(1:end-1), Einf, 's-','LineWidth',1.6,'DisplayName','L^\infty_t error');

% slope-2 reference line (anchored to finest E2)
href = hs(1:end-1);
Cref = E2(end) / (href(end)^2 + eps);
loglog(href, Cref*href.^2, ':','LineWidth',1.2,'DisplayName','O(h^2) reference');

grid on; box on;
xlabel('h = \Delta x = \Delta t');
ylabel('error in u(1,t)');
title(sprintf('Convergence of right-end trace (\\beta=%.3g, \\kappa=%.3g)', beta, kappa), 'Interpreter','tex');
legend('Location','south east');

end

% ======================================================================
function [t, u_right] = run_beam(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% One run of the implicit FD scheme with DNC at x=1.
M  = round(L/dx);              % nodes 0..M
N  = round(T/dt);              % times 0..N
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

phi = @(tt) 0.0;               % u(0,t)=0
u0  = @(xx) 0.0;
v0  = @(xx) 0.0;

gamma = alpha2*(dt^2)/(dx^4);

u = zeros(M+1, N+1);           % u(i+1,j+1) ~ u_i^j
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0;                    % Dirichlet at x=0

for j = 1:(N-1)
    % contact decision from previous time level at right end
    in_contact = (u(M+1, j+1) <= y_minus);

    % Unknowns: w = [u_1^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);  % phi=0 eliminates 4γ*phi

    % interior i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i = M-1
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)     = dt^2*f0 + 2*u(M, j+1) - u(M, j);

    % i = M (contact vs no-contact)
    if ~in_contact
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        c = (dx^3/alpha2) * (kappa + beta/dt);
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma*(c + 1);
        b(M)      =  dt^2*f0 ...
                   + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
                   -   u(M+1, j) ...
                   + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    % Solve and write back
    w = A \ b;
    u(1,   j+2) = 0;
    u(2:end,j+2) = w;
end

u_right = u(end,:);
end
