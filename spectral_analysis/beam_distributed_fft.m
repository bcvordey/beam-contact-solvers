function beam_contact_distributed_DNC_fft()


%% ---- Parameters (edit as needed) ----
L  = 1.0;   T  = 20.0;
dx = 0.01;  dt = 0.01;

alpha2  = 1;          % α^2
f0      = -0.2;       % constant forcing f(x,t)=f0


kappaVec = 10; %[0.001, 1, 10, 50, 100, 200, 500, 1000];    

betaVec = [0.001, 0.01, 0.1, 0.5, 1, 5]; 

x_probe = 1;  

% Obstacle profile ψ(x)
psiFun  = @(x) -0.02 + 0*x;   % constant obstacle: ψ(x) = -0.02

phi  = @(tt) 0.0;             % boundary excitation at x=0
u0f  = @(xx) 0.0;             % initial displacement
v0f  = @(xx) 0.0;             % initial velocity

% Force both to row vectors
kappaVec = kappaVec(:).';
betaVec  = betaVec(:).';

nk = numel(kappaVec);
nb = numel(betaVec);

%% ---- Solve & Plot depending on (nk, nb) configuration ----

% First, solve once to get grids and probe index
Sref = solve_single_kappa_psi_minimal(kappaVec(1), betaVec(1), ...
    L,T,dx,dt, alpha2,psiFun,f0, phi,u0f,v0f);
x  = Sref.x; 
t  = Sref.t;
dt = t(2) - t(1);

[~, idx_probe] = min(abs(x - x_probe));
x_probe_eff = x(idx_probe);

% =========================
% CASE A: vary kappa, fixed beta
% =========================
if nk > 1 && nb == 1
    beta = betaVec(1);

    figure('Color','w','Name','FFT vs kappa at selected x');
    hold on; grid on;

    for ik = 1:nk
        kappa = kappaVec(ik);

        S = solve_single_kappa_psi_minimal(kappa, beta, ...
            L,T,dx,dt, alpha2,psiFun,f0, phi,u0f,v0f);

        u_probe = S.u(idx_probe,:);
        [f, P1] = single_sided_fft_plain(u_probe, dt);

        plot(f, P1, 'LineWidth', 1, ...
            'DisplayName', sprintf('\\kappa = %.3g', kappa));
    end

    xlabel('$Frequency \ (Hz)$', 'Interpreter','latex');
    ylabel('$Amplitude \ |U(f)|$', 'Interpreter','latex');
    title(sprintf(['Single-Sided Amplitude Spectrum at $x = %.3f$, ' ...
                   'Varying $\\kappa$ ($\\beta = %.3g$)'], ...
                   x_probe_eff, beta), ...
          'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex');
    grid off; box on;
    xlim([0, 5]);
    ylim ([0, 0.014]);
    return;
end

% =========================
% CASE B: vary beta, fixed kappa
% =========================
if nk == 1 && nb > 1
    kappa = kappaVec(1);

    figure('Color','w','Name','FFT vs beta at selected x');
    hold on; grid on;

    for ib = 1:nb
        beta = betaVec(ib);

        S = solve_single_kappa_psi_minimal(kappa, beta, ...
            L,T,dx,dt, alpha2,psiFun,f0, phi,u0f,v0f);

        u_probe = S.u(idx_probe,:);
        [f, P1] = single_sided_fft_plain(u_probe, dt);

        plot(f, P1, 'LineWidth', 1, ...
            'DisplayName', sprintf('\\beta = %.3g', beta));
    end

    xlabel('$Frequency \ (Hz)$', 'Interpreter','latex');
    ylabel('$Amplitude \ |U(f)|$', 'Interpreter','latex');
    title(sprintf(['Single-Sided Amplitude Spectrum at $x = %.3f$, ' ...
                   'Varying $\\beta$ ($\\kappa = %.3g$)'], ...
                   x_probe_eff, kappa), ...
          'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex');
    grid off; box on;
    xlim([0, 5]);
    ylim([0, 0.01]);
    return;
end

% =========================
% CASE C: both vectors (nk>1, nb>1)
% For each kappa(i), plot all betas on its own figure
% =========================
for ik = 1:nk
    kappa = kappaVec(ik);

    figure('Color','w','Name', ...
        sprintf('FFT at selected x for kappa=%.3g, varying beta', kappa));
    hold on; grid on;

    for ib = 1:nb
        beta = betaVec(ib);

        S = solve_single_kappa_psi_minimal(kappa, beta, ...
            L,T,dx,dt, alpha2,psiFun,f0, phi,u0f,v0f);

        u_probe = S.u(idx_probe,:);
        [f, P1] = single_sided_fft_plain(u_probe, dt);

        plot(f, P1, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('\\beta = %.3g', beta));
    end

    xlabel('Frequency (Hz)', 'Interpreter','latex');
    ylabel('$|U(f)|$ at selected $x$', 'Interpreter','latex');
    title(sprintf(['Distributed DNC beam: FFT at $x \\approx %.3f$ ' ...
                   'for varying $\\beta$ (with $\\kappa = %.3g$)'], ...
                   x_probe_eff, kappa), ...
          'Interpreter','latex');
    legend('Location','northeast', 'Interpreter','latex');
    xlim([0, 8]);
    ylim([0, 0.01]);
end

end

function S = solve_single_kappa_psi_minimal(kappa, beta, ...
    L,T,dx,dt, alpha2,psiFun,f0, phi,u0f,v0f)
% Solve the distributed DNC beam with obstacle ψ(x) for a single (kappa, beta).
% Returns x, t, u, psi.

M = round(L/dx);  
N = round(T/dt);

x = linspace(0,L,M+1).';  
t = linspace(0,T,N+1);

psi = psiFun(x);                        % obstacle profile at nodes

u = zeros(M+1, N+1);
gamma = alpha2*(dt^2)/(dx^4);

% Initial data (2nd-order start)
u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);  % Dirichlet at x=0

for j = 1:(N-1)
    % Contact test against ψ(x) at time level j+1
    I = (u(:, j+1) <= psi);  
    I(1) = false;   % left boundary not in contact (Dirichlet)

    Ivec = double(I(2:end));  % column vector of 0/1

    % DNC diagonal contributions (elementwise)
    diag_add = dt^2 * (kappa .* Ivec + (beta/dt) .* Ivec);
    rhs_add  = dt^2 * (kappa .* psi(2:end) .* Ivec ...
                      + (beta/dt) .* Ivec .* u(2:end, j+1));

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % Row i=1 (left boundary closure)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));

    A(1,1) = A(1,1) + diag_add(1);
    b(1)   = b(1)   + rhs_add(1);

    % Interior i=2..M-2
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

    % Right end i=M-1
    if M >= 2
        A(M-1, M-3) = A(M-1, M-3) + gamma;
        A(M-1, M-2) = A(M-1, M-2) - 4*gamma;
        A(M-1, M-1) = A(M-1, M-1) + (1 + 5*gamma);
        A(M-1, M  ) = A(M-1, M  ) - 2*gamma;
        b(M-1) = b(M-1) + dt^2*f0 + 2*u(M, j+1) - u(M, j);

        A(M-1, M-1) = A(M-1, M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % Right end i=M
    A(M, M-2) = A(M, M-2) + 2*gamma;
    A(M, M-1) = A(M, M-1) - 4*gamma;
    A(M, M  ) = A(M, M  ) + (1 + 2*gamma);
    b(M) = b(M) + dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);

    A(M, M) = A(M, M) + diag_add(M);
    b(M)    = b(M)    + rhs_add(M);

    % Solve and write back
    w = A \ b;
    u(1,   j+2)   = phi(t(j+2));
    u(2:end,j+2) = w;
end

S.x    = x;
S.t    = t;
S.u    = u;
S.psi  = psi;
S.kappa = kappa;
S.beta  = beta;

end


function [f, P1] = single_sided_fft_plain(y, dt)
% Plain FFT (no detrend, no window), single-sided amplitude spectrum.

    % Ensure column vector
    if size(y,1) == 1
        y = y.';
    end

    N  = numel(y);
    Fs = 1/dt;

    Y  = fft(y);
    P2 = abs(Y/N);                % two-sided amplitude
    K  = floor(N/2) + 1;          % number of positive frequencies
    P1 = P2(1:K);                 % single-sided

    if K > 2
        P1(2:end-1) = 2*P1(2:end-1);  % double interior amplitudes
    end

    f = Fs*(0:K-1)/N;             % frequency axis in Hz
end
