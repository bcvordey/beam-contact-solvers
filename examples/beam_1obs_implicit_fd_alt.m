function [x,t,u] = beam_contact_implicit_FD()
% Implicit second-order FD for the Euler–Bernoulli beam with unilateral contact
% Unknown vector at each step: [u_1^{j+1}, ..., u_M^{j+1}]^T
% Left BC: u(0,t)=phi(t)=0 and u_x(0,t)=0 (ghost: u_{-1}=u_1)
% Right BC: u_xx(1,t)=0 and u_xxx(1,t) given by contact law (or 0 in no-contact)

%% Parameters (from the prompt)
L = 1.0;            % beam length
T = 10.0;            % final time
dx = 0.01;          % space step
dt = 0.01;          % time step
M  = round(L/dx);   % number of spatial subintervals => nodes are 0..M
N  = round(T/dt);   % number of time steps       => times are 0..N

alpha2 = 1;     % alpha^2
f0     = -0.2;      % constant forcing
y_minus = -0.02;    % obstacle
kappa  = 1e-2;      % contact stiffness (adjustable)
beta   = 0.1;      % contact damping   (adjustable)

phi = @(tt) 0.0;    % left displacement boundary
u0  = @(xx) 0.0;    % initial displacement
v0  = @(xx) 0.0;    % initial velocity

gamma = alpha2*(dt^2)/(dx^4);  % scheme parameter

%% Grids & storage (i=0..M, j=0..N -> MATLAB indices +1)
x = linspace(0,L,M+1).';   % column
t = linspace(0,T,N+1);     % row
u = zeros(M+1, N+1);       % u(i+1, j+1) ~ u_i^j

% Initial conditions (u^0, u^1)
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);

% Enforce Dirichlet at x=0 for all times we touch
u(1, :) = phi(t);

%% Time stepping: build and solve A w = b at each j (compute u^{j+1})
for j = 1:(N-1)   % j=1..N-1 (we have u^0 and u^1; compute up to u^N)

    % Decide contact at the RIGHT end using u_M^j (i.e., u(M+1, j+1))
    in_contact = (u(M+1, j+1) <= y_minus);

    % Assemble A (M x M) and b (M x 1)
    % Unknown vector w corresponds to [u_1^{j+1}; u_2^{j+1}; ... ; u_M^{j+1}]
    A = spalloc(M, M, 5*M);  % five-diagonal structure plus a couple special rows
    b = zeros(M,1);

    %------- Row for i=1 (uses u_{-1}^{j+1} = u_1^{j+1} via u_x(0,t)=0) -------
    % (1 + 7γ) u_1 - 4γ u_2 + γ u_3 = dt^2 f + 2 u_1^j - u_1^{j-1} + 4γ φ(t_{j+1})
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(1+1, j+1) - u(1+1, j) + 4*gamma*phi(t(j+1));  % phi=0 here

    %------- Interior: i = 2 .. M-2 using the standard five-point stencil -------
    % γ u_{i-2} - 4γ u_{i-1} + (1+6γ) u_i - 4γ u_{i+1} + γ u_{i+2} = dt^2 f + 2u_i^j - u_i^{j-1}
    for i = 2:(M-2)
        r = i;  % row index in A/b

        % Left two-off-diagonal (i-2). If i=2, this hits u_0^{j+1}=phi(t_{j+1}) (known).
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            % move γ*u_0^{j+1} to RHS
            b(r) = b(r) - gamma*phi(t(j+1));  % zero here
        end

        % The remaining stencil entries are within unknowns
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        % RHS
        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    %------- i = M-1 row (uses u_{M+1} = 2u_M - u_{M-1}) already eliminated) -------
    % γ u_{M-3} - 4γ u_{M-2} + (1+5γ) u_{M-1} - 2γ u_M = dt^2 f + 2 u_{M-1}^j - u_{M-1}^{j-1}
    if M >= 2
        A(M-1, M-3) = gamma;
        A(M-1, M-2) = -4*gamma;
        A(M-1, M-1) = 1 + 5*gamma;
        A(M-1, M  ) = -2*gamma;
        b(M-1) = dt^2*f0 + 2*u(M-1+1, j+1) - u(M-1+1, j);
    end

    %------- i = M row: no-contact vs contact -------
    if ~in_contact
        % No contact:
        % 2γ u_{M-2} - 4γ u_{M-1} + (1+2γ) u_M = dt^2 f + 2 u_M^j - u_M^{j-1}
        A(M, M-2) = 2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) = 1 + 2*gamma;
        b(M) = dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        % Contact:
        % 2γ u_{M-2} - 4γ u_{M-1} + (1 + 2γ*( dx^3/α^2*(κ + β/dt) + 1 )) u_M
        %   = dt^2 f + 2*(1 + γ*(β*dx^3)/(α^2*dt)) u_M^j - u_M^{j-1} + 2γ*(dx^3/α^2)*κ*y_-
        c = (dx^3/alpha2)*(kappa + beta/dt);
        A(M, M-2) = 2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) = 1 + 2*gamma*(c + 1);
        b(M) = dt^2*f0 ...
             + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             - u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    % Solve for w = [u_1^{j+1}; ... ; u_M^{j+1}]
    w = A \ b;

    % Write back: u_0^{j+1} = phi(t_{j+1}), u_i^{j+1} = w(i) for i=1..M
    u(1,   j+2) = phi(t(j+2));     % Dirichlet (zero here)
    u(2:end,j+2) = w;              % unknowns
end

%% Plot: time on x-axis, right end displacement u(1,t) on y-axis (x=1 -> i=M)
figure;
plot(t, u(end,:), 'LineWidth', 1.6); grid on;
hObs = yline(y_minus, '--','LineWidth',1.6, 'DisplayName', sprintf('y_- = %.3g', y_minus));
xlabel('time');
ylabel('u(1,t)  (right end)');
title('Right-end displacement vs time (implicit FD with contact)');
end


