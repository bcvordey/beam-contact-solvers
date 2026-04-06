function beam_dnc_distributed_fft_2D_suite(kappa, beta)
% 2D FFT / spectral analysis for the distributed obstacle DNC beam
% with time-dependent sinusoidal forcing f(t).
%
% Usage:
%   beam_dnc_distributed_fft_2D_suite();                 % default
%   beam_dnc_distributed_fft_2D_suite(1e3, 0.5);         % custom kappa, beta

%% ---- Parameters ----
L  = 1.0;   T  = 10.0;
dx = 0.01;  dt = 0.01;

alpha2 = 1.0;

if nargin < 1 || isempty(kappa), kappa = 1; end
if nargin < 2 || isempty(beta),  beta  = 0.01;  end

% ---- Time-dependent forcing f(t) = sum_k F_k sin(omega_k t) ----
Fvec    = [0.4, 0.6];        % amplitudes
omegavec = [11.0, 12.0];      % angular frequencies (rad/s)
f_fun   = @(tt) Fvec(1)*sin(omegavec(1)*tt) + ...
                Fvec(2)*sin(omegavec(2)*tt);
% You can edit Fvec and omegavec to pick any combination of frequencies.

% Obstacle profile psi(x)
psiFun = @(x) -0.02 + 0*x;   % constant floor; modify if needed

% Boundary/initial conditions
phi  = @(tt) 0.0;
u0f  = @(xx) 0.0;
v0f  = @(xx) 0.0;

%% ---- Solve PDE for given (kappa, beta) ----
S = solve_distributed_dnc(kappa, beta, ...
    L,T,dx,dt, alpha2,psiFun,f_fun, phi,u0f,v0f);

x   = S.x;
t   = S.t;
u   = S.u;
psi = S.psi;
dt  = t(2) - t(1);
Fs  = 1/dt;

fprintf('Solved distributed DNC beam (sin forcing): kappa=%.3g, beta=%.3g, M=%d, N=%d\n', ...
        kappa, beta, numel(x)-1, numel(t)-1);

%% 1) 2D FFT of u(x,t): space-time spectrum
plot_2d_fft_ut(x, t, u, kappa, beta);

%% 2) Contact vs free region FFT (temporal)
plot_fft_contact_vs_free(t, u, psi, dt, kappa, beta);

%% 3) Spectrogram at a chosen spatial point
x_probe = 1.0;  % tip (change to 0.5, etc., if desired)
plot_spectrogram_point(x, t, u, x_probe, Fs, kappa, beta);

%% 4) Mode-energy vs temporal frequency
plot_mode_energy_vs_freq(u, x, dt, kappa, beta);

%% 5) Mode energy vs spatial mode index (modal suppression)
plot_mode_energy_vs_mode(u, x, t, kappa, beta);

end


% =======================================================================
%                      SOLVER: distributed DNC beam
% =======================================================================
function S = solve_distributed_dnc(kappa, beta, ...
    L,T,dx,dt, alpha2,psiFun,f_fun, phi,u0f,v0f)
% Distributed DNC solver with time-dependent forcing f(t).

M = round(L/dx);
N = round(T/dt);

x = linspace(0,L,M+1).';
t = linspace(0,T,N+1);

psi = psiFun(x);
u   = zeros(M+1, N+1);

gamma = alpha2*(dt^2)/(dx^4);

% Initial data
u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);  % Dirichlet at x=0

for j = 1:(N-1)
    % Time-dependent load at time level t_{j+1}
    fjp1 = f_fun(t(j+1));

    % Contact indicator at time level j+1
    I = (u(:, j+1) <= psi);
    I(1) = false;  % left boundary excluded

    % DNC diagonal contributions
    diag_add = dt^2*(kappa*double(I(2:end)) + (beta/dt)*double(I(2:end)));
    rhs_add  = dt^2*(kappa*psi(2:end).*double(I(2:end)) + ...
                     (beta/dt)*double(I(2:end)).*u(2:end,j+1));

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % Row i=1 (left closure)
    A(1,1) = 1 + 7*gamma;
    if M>=2, A(1,2) = -4*gamma; end
    if M>=3, A(1,3) =  gamma;   end

    b(1) = dt^2*fjp1 + 2*u(2,j+1) - u(2,j) + 4*gamma*phi(t(j+1));
    A(1,1) = A(1,1) + diag_add(1);
    b(1)   = b(1)   + rhs_add(1);

    % Interior rows i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r,i-2) = A(r,i-2) + gamma;
        end
        A(r,i-1) = A(r,i-1) - 4*gamma;
        A(r,i)   = A(r,i)   + (1 + 6*gamma);
        A(r,i+1) = A(r,i+1) - 4*gamma;
        A(r,i+2) = A(r,i+2) + gamma;

        b(r) = b(r) + dt^2*fjp1 + 2*u(i+1,j+1) - u(i+1,j);

        A(r,i) = A(r,i) + diag_add(i);
        b(r)   = b(r)   + rhs_add(i);
    end

    % Right end i=M-1
    if M>=2
        A(M-1,M-3) = A(M-1,M-3) + gamma;
        A(M-1,M-2) = A(M-1,M-2) - 4*gamma;
        A(M-1,M-1) = A(M-1,M-1) + (1 + 5*gamma);
        A(M-1,M)   = A(M-1,M)   - 2*gamma;

        b(M-1) = b(M-1) + dt^2*fjp1 + 2*u(M,j+1) - u(M,j);

        A(M-1,M-1) = A(M-1,M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % Right end i=M
    A(M,M-2) = A(M,M-2) + 2*gamma;
    A(M,M-1) = A(M,M-1) - 4*gamma;
    A(M,M)   = A(M,M)   + (1 + 2*gamma);

    b(M) = b(M) + dt^2*fjp1 + 2*u(M+1,j+1) - u(M+1,j);

    A(M,M) = A(M,M) + diag_add(M);
    b(M)   = b(M)   + rhs_add(M);

    % Solve
    w = A \ b;

    u(1,j+2)     = phi(t(j+2));
    u(2:end,j+2) = w;
end

S.x     = x;
S.t     = t;
S.u     = u;
S.psi   = psi;
S.kappa = kappa;
S.beta  = beta;
end


% =======================================================================
%                          1) 2D FFT of u(x,t)
% =======================================================================
function plot_2d_fft_ut(x, t, u, kappa, beta)
% 2D FFT in space (mode index) and time (frequency).

[Nx, Nt] = size(u);
dx = x(2)-x(1);
dt = t(2)-t(1);
Fs = 1/dt;

% Remove mean in time to focus on oscillatory content (optional but helpful)
u0 = u - mean(u,2);

% 2D FFT
Uhat = fft2(u0);

% Frequencies (time) – keep positive half
Kf  = floor(Nt/2)+1;
f_t = (0:Kf-1)*Fs/Nt;

% Spatial "modes" (index) – we use 0..Nx-1 here
Km   = floor((Nx)/2)+1;
midx = 0:(Km-1);   % mode index (0 = spatial average)

% Take magnitude
P2 = abs(Uhat/Nt);      % not too worried about exact scaling here
P  = P2(1:Km,1:Kf);     % positive time freq, first half spatial modes

figure('Color','w','Name','2D FFT u(x,t)');
imagesc(f_t, midx, log10(P + 1e-12));
set(gca,'YDir','normal');
xlabel('Temporal frequency f (Hz)');
ylabel('Spatial mode index');
title(sprintf('2D FFT of u(x,t), \\kappa=%.3g, \\beta=%.3g', kappa, beta), ...
      'Interpreter','tex');
colorbar; ylabel(colorbar,'log_{10}|U(k,\omega)|');
end

% =======================================================================
%                 2) Contact vs Free: temporal FFT comparison
% =======================================================================
function plot_fft_contact_vs_free(t, u, psi, dt, kappa, beta)
% Average temporal FFT over contact vs free nodes.

[Nx, Nt] = size(u);
Fs = 1/dt;

contact = (u <= psi);
ever_contact = any(contact,2);
free_nodes   = ~ever_contact;

Uc = u(ever_contact,:);
Uf = u(free_nodes,:);

[f_c, Pc_mean] = avg_fft_over_rows(Uc, dt);
[f_f, Pf_mean] = avg_fft_over_rows(Uf, dt);

figure('Color','w','Name','FFT Contact vs Free');
hold on; grid on;

if ~isempty(Pc_mean)
    plot(f_c, Pc_mean, 'LineWidth',1.8, 'DisplayName','Contact region');
end
if ~isempty(Pf_mean)
    plot(f_f, Pf_mean, '--', 'LineWidth',1.8, 'DisplayName','Free region');
end

xlabel('Frequency (Hz)');
ylabel('Average |U(f)|');
title(sprintf('Contact vs Free FFT (\\kappa=%.3g, \\beta=%.3g)', kappa, beta), ...
      'Interpreter','tex');
legend('Location','northeast');
xlim([0, Fs/2]);
end

function [f, meanP1] = avg_fft_over_rows(U, dt)
if isempty(U)
    f = [];
    meanP1 = [];
    return;
end

[M, N] = size(U);
Fs = 1/dt;
K  = floor(N/2)+1;

meanP1 = zeros(1,K);

for i = 1:M
    y  = U(i,:);
    Y  = fft(y);
    P2 = abs(Y/N);
    P1 = P2(1:K);
    if K>2
        P1(2:end-1) = 2*P1(2:end-1);
    end
    meanP1 = meanP1 + P1;
end

meanP1 = meanP1 / M;
f = (0:K-1)*Fs/N;
end

% =======================================================================
%                    3) Spectrogram at a chosen point
% =======================================================================
function plot_spectrogram_point(x, t, u, x_probe, Fs, kappa, beta)
% STFT / spectrogram at a given spatial location x_probe.

[~, idx] = min(abs(x - x_probe));
y = u(idx,:);

figure('Color','w','Name','Spectrogram at x probe');

% Use Signal Processing Toolbox spectrogram if available
if exist('spectrogram','file') == 2
    win      = 256;
    noverlap = round(0.6*win);
    nfft     = 512;
    [S,F,Tstft] = spectrogram(y - mean(y), win, noverlap, nfft, Fs, 'yaxis');
    imagesc(Tstft, F, 20*log10(abs(S)+1e-12));
    axis xy;
    xlabel('Time'); ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram at x=%.2f (\\kappa=%.3g, \\beta=%.3g)', ...
          x(idx), kappa, beta), 'Interpreter','tex');
    cb = colorbar; ylabel(cb,'Magnitude (dB)');
else
    plot(t, y, 'LineWidth',1.5);
    xlabel('Time'); ylabel('u(x,t)');
    title(sprintf('No spectrogram(): time trace at x=%.2f', x(idx)));
end

end

% =======================================================================
%           4) Mode-energy heatmap: E(omega) vs frequency
% =======================================================================
function plot_mode_energy_vs_freq(u, x, dt, kappa, beta)
% Integrate |Û(k,omega)|^2 over space → energy vs frequency.

[Nx, Nt] = size(u);
dx = x(2)-x(1);
Fs = 1/dt;

% 2D FFT
Uhat = fft2(u);
Kf   = floor(Nt/2)+1;
Upos = Uhat(:,1:Kf);

% Energy vs frequency: sum over space modes
E = sum(abs(Upos).^2, 1) * dx;

f = (0:Kf-1)*Fs/Nt;

figure('Color','w','Name','Energy vs frequency');
plot(f, E, 'LineWidth',1.8);
xlabel('Frequency (Hz)');
ylabel('E(\omega) (arbitrary units)');
title(sprintf('Energy vs frequency (\\kappa=%.3g, \\beta=%.3g)', ...
      kappa, beta), 'Interpreter','tex');
grid on;
xlim([0, Fs/2]);
end

% =======================================================================
%   5) Mode energy vs spatial mode index (modal suppression picture)
% =======================================================================
function plot_mode_energy_vs_mode(u, x, t, kappa, beta)
% Spatial FFT at each time, then integrate |Û(k,t)|^2 over time.

[Nx, Nt] = size(u);
dt = t(2)-t(1);

Ux = fft(u, [], 1);   % FFT in x only, each column is time

% Energy per spatial mode: integrate over time
E_mode = sum(abs(Ux).^2, 2) * dt;   % size Nx x 1
Km = floor(Nx/2)+1;
midx = 0:(Km-1);
E_plot = E_mode(1:Km);

figure('Color','w','Name','Mode energy vs mode index');
stem(midx, E_plot, 'filled');
xlabel('Spatial mode index');
ylabel('Mode energy (arbitrary units)');
title(sprintf('Mode energy vs mode index (\\kappa=%.3g, \\beta=%.3g)', ...
      kappa, beta), 'Interpreter','tex');
grid on;
end
