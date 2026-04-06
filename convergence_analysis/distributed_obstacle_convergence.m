function convTable = distributed_obstacle_convergence()

    %% Fixed parameters
    L      = 1.0;
    T      = 2.0;
    alpha2 = 1.0;
    kappa  = 1.0;
    beta   = 1e-3;
    f0     = -0.2;

    % Obstacle profile ψ(x): constant (change if desired)
    psiFun = @(x) -0.02 + 0*x;

    % Boundary/initial data
    phi = @(tt) 0.0;   % u(0,t)
    u0  = @(xx) 0.0;   % initial displacement
    v0  = @(xx) 0.0;   % initial velocity

    %% Reference solution: dx = dt = 2^-12
    k_ref  = 12;
    dx_ref = 2^(-k_ref);
    dt_ref = dx_ref;

    [x_ref, t_ref, u_ref] = solve_distributed_DNC( ...
        L,T,dx_ref,dt_ref, alpha2,kappa,beta,f0, psiFun,phi,u0,v0);

    u_ref_T = u_ref(:, end);  % reference at final time T=2

    %% Test meshes: k = 4,...,10
    kList = 4:10;
    nK    = numel(kList);

    dxList   = zeros(nK,1);
    errLinf  = zeros(nK,1);
    rate     = NaN(nK,1);   % first rate is undefined

    for idx = 1:nK
        k_val = kList(idx);
        dx    = 2^(-k_val);
        dt    = dx;           % dx = dt as required

        dxList(idx) = dx;

        % Solve on coarse mesh
        [x_coarse, t_coarse, u_coarse] = solve_distributed_DNC( ...
            L,T,dx,dt, alpha2,kappa,beta,f0, psiFun,phi,u0,v0);

        % Final time slice at T=2 for coarse mesh
        u_coarse_T = u_coarse(:, end);   % size (M+1,1)

        % Restrict reference solution to coarse grid nodes
        % ratio = dx / dx_ref = 2^(12-k_val) is an integer
        ratio   = round(dx / dx_ref);
        M_coarse = numel(x_coarse) - 1;

        fine_idx = 1:ratio:(M_coarse*ratio + 1);   % matching nodes
        u_ref_on_coarse = u_ref_T(fine_idx);

        % Sanity check (optional)
        % max(abs(x_coarse - x_ref(fine_idx)))

        % Linf error at final time
        errLinf(idx) = max(abs(u_ref_on_coarse - u_coarse_T));
    end

    %% Compute convergence rates
    % rate_j = log( ||~ - u_{dx/2}|| / ||~ - u_dx|| ) / log(0.5)
    % Here dx_j = dx_{j-1}/2, so we use errors at successive mesh sizes.
    for idx = 2:nK
        ratioErr   = errLinf(idx) / errLinf(idx-1);
        rate(idx)  = log(ratioErr) / log(0.5);
    end

    %% Build and display table
    convTable = table(dxList, errLinf, rate, ...
        'VariableNames', {'dx', 'LinfError', 'Rate'});

    fprintf('\nConvergence table for distributed DNC scheme (T = 2):\n');
    fprintf('kappa = %.3g, beta = %.3g, f = %.3g, alpha^2 = %.3g\n', ...
            kappa, beta, f0, alpha2);
    fprintf('%12s  %15s  %10s\n', 'dx', 'LinfError', 'Rate');
    fprintf('-----------------------------------------------\n');
    for idx = 1:nK
        if isnan(rate(idx))
            fprintf('%12.4e  %15.8e      %6s\n', ...
                dxList(idx), errLinf(idx), '---');
        else
            fprintf('%12.4e  %15.8e  %10.4f\n', ...
                dxList(idx), errLinf(idx), rate(idx));
        end
    end
    fprintf('\n');

end

% ======================================================================
%   Implicit FD solver for the distributed DNC obstacle problem
% ======================================================================
function [x, t, u] = solve_distributed_DNC( ...
    L,T,dx,dt, alpha2,kappa,beta,f0, psiFun,phi,u0f,v0f)

    M = round(L/dx);
    N = round(T/dt);

    x = linspace(0,L,M+1).';
    t = linspace(0,T,N+1);

    psi = psiFun(x);

    u = zeros(M+1, N+1);
    gamma = alpha2 * (dt^2) / (dx^4);

    % Initial data (second-order start)
    u(:,1) = u0f(x);
    u(:,2) = u0f(x) + dt*v0f(x);
    u(1,:) = phi(t);  % Dirichlet at x=0

    for j = 1:(N-1)
        % Contact indicator at time level j+1
        I = (u(:, j+1) <= psi);
        I(1) = false;  % left boundary not in contact

        Ivec = double(I(2:end));  % length M

        % DNC contributions on interior nodes i=1..M
        diag_add = dt^2 * (kappa .* Ivec + (beta/dt) .* Ivec);
        rhs_add  = dt^2 * (kappa .* psi(2:end) .* Ivec ...
                          + (beta/dt) .* Ivec .* u(2:end, j+1));

        A = spalloc(M, M, 5*M);
        b = zeros(M,1);

        % i = 1 row (left closure with ghost)
        A(1,1) = 1 + 7*gamma;
        if M >= 2, A(1,2) = -4*gamma; end
        if M >= 3, A(1,3) =  gamma;   end

        b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));
        A(1,1) = A(1,1) + diag_add(1);
        b(1)   = b(1)   + rhs_add(1);

        % Interior rows i = 2..M-2
        for i = 2:(M-2)
            r = i;
            if i-2 >= 1
                A(r, i-2) = A(r, i-2) + gamma;
            else
                b(r) = b(r) - gamma*phi(t(j+1));
            end
            A(r, i-1) = A(r, i-1) - 4*gamma;
            A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
            A(r, i+1) = A(r, i+1) - 4*gamma;
            A(r, i+2) = A(r, i+2) + gamma;

            b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);

            A(r, i) = A(r, i) + diag_add(i);
            b(r)    = b(r)    + rhs_add(i);
        end

        % i = M-1 row
        if M >= 2
            A(M-1, M-3) = A(M-1, M-3) + gamma;
            A(M-1, M-2) = A(M-1, M-2) - 4*gamma;
            A(M-1, M-1) = A(M-1, M-1) + (1 + 5*gamma);
            A(M-1, M  ) = A(M-1, M  ) - 2*gamma;

            b(M-1) = b(M-1) + dt^2*f0 + 2*u(M, j+1) - u(M, j);

            A(M-1, M-1) = A(M-1, M-1) + diag_add(M-1);
            b(M-1)      = b(M-1)      + rhs_add(M-1);
        end

        % i = M row (right end)
        A(M, M-2) = A(M, M-2) + 2*gamma;
        A(M, M-1) = A(M, M-1) - 4*gamma;
        A(M, M  ) = A(M, M  ) + (1 + 2*gamma);

        b(M) = b(M) + dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);

        A(M, M) = A(M, M) + diag_add(M);
        b(M)    = b(M)    + rhs_add(M);

        % Solve and update
        w = A \ b;
        u(1,   j+2)   = phi(t(j+2));
        u(2:end,j+2) = w;
    end
end
