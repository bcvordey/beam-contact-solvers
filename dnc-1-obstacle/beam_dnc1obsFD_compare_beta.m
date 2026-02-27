function beam_dnc1obsFD_compare_beta()
% Compare right-end displacement for multiple beta values on one plot.
% Uses the implicit FD scheme with unilateral contact at x=1.

%% Fixed parameters
L = 1.0; T = 10.0; dx = 0.01; dt = 0.01;
alpha2 = 1;
f0 = -0.2;
y_minus = -0.02;                 % <-- obstacle level requested
kappa = 10;                   % keep your previous choice (shown in title)
betas = [1e-3, 1e-2, 1e-1, 1, 5];     % the requested beta set

figure; hold on;
leg = cell(0,1);

for b = betas
    [t, u_right] = simulate_right_end(b, kappa, y_minus, L, T, dx, dt, alpha2, f0);
    plot(t, u_right, 'LineWidth', 1.6, 'DisplayName', sprintf('\\beta = %g', b));
end

% Obstacle line (broken/dashed) on same graph
hObs = yline(y_minus, '--','LineWidth',1.6, 'DisplayName', sprintf('y_- = %.3g', y_minus));

grid on; xlim([0 T]);
xlabel('$Time$','Interpreter','latex');
ylabel('$u(1,t)$','Interpreter','latex'); grid off;
title(sprintf('Right Tip Displacement vs Time ($\\kappa = %.3g$)', kappa), 'Interpreter','latex');
legend('Location','southeast');
end

% -------------------------------------------------------------------------
function [t, u_right] = simulate_right_end(beta, kappa, y_minus, L, T, dx, dt, alpha2, f0)
% Runs one simulation for a given beta; returns time vector and right-end trace.

M  = round(L/dx);              % nodes: 0..M
N  = round(T/dt);              % times: 0..N
x  = linspace(0,L,M+1).';
t  = linspace(0,T,N+1);

phi = @(tt) 0.0;               % u(0,t) = 0
u0  = @(xx) 0.0;               % initial displacement
v0  = @(xx) 0.0;               % initial velocity

gamma = alpha2*(dt^2)/(dx^4);

% Storage (u(i+1,j+1) ~ u_i^j)
u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = 0;                    % enforce Dirichlet at x=0

for j = 1:(N-1)                % compute u^{j+1} from u^{j}, u^{j-1}
    in_contact = (u(M+1, j+1) <= y_minus);

    % Unknowns are w = [u_1^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % i=1 row (uses u_{-1}=u_1)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);  % phi=0 so the +4γ*phi term vanishes

    % Interior rows: i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            % would add gamma*phi to RHS; phi=0 -> no effect
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i=M-1 row
    A(M-1, M-3) = gamma;
    A(M-1, M-2) = -4*gamma;
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1)     = dt^2*f0 + 2*u(M, j+1) - u(M, j);

    % i=M row: contact vs no-contact
    if ~in_contact
        % No contact:
        A(M, M-2) =  2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) =  1 + 2*gamma;
        b(M)      =  dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        % Contact:
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
    u(1,   j+2) = 0;       % phi = 0
    u(2:end,j+2) = w;
end

u_right = u(end,:);         % trace at x=1
end
