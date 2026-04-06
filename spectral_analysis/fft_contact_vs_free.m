function beam_dnc_distributed_fft_contact_vs_free()
% Compare FFT of the distributed DNC beam between CONTACT and FREE regions.
%
% Usage:
%   beam_dnc_distributed_fft_contact_vs_free(1);
%   beam_dnc_distributed_fft_contact_vs_free([0.1 0.5 1 5]);

%% ---- Parameters ----
L  = 1.0;   T  = 10.0;
dx = 0.01;  dt = 0.01;

alpha2  = 1;          
kappa   = 1;       
f0      = -0.2;       


betaVec = [0.001,0.01,0.1,1,5];

betaVec = betaVec(:).';

psiFun  = @(x) -0.02 + 0*x;  

phi  = @(tt) 0.0;
u0f  = @(xx) 0.0;
v0f  = @(xx) 0.0;

%% ---- Loop over beta ----
for k = 1:numel(betaVec)
    beta = betaVec(k);

    % Solve PDE
    S = beam_dnc_distributed_solver(beta, ...
        L,T,dx,dt, alpha2,kappa,psiFun,f0, phi,u0f,v0f);

    x   = S.x;
    t   = S.t;
    u   = S.u;
    psi = S.psi;
    dt  = t(2) - t(1);

    % Contact classification
    contact      = (u <= psi);
    ever_contact = any(contact, 2);
    free_nodes   = ~ever_contact;

    Uc = u(ever_contact,:);
    Uf = u(free_nodes,:);

    % FFT for contact vs free
    [f_c, Pc_mean] = fft_average_over_nodes(Uc, dt);
    [f_f, Pf_mean] = fft_average_over_nodes(Uf, dt);

    % ---- Plot comparison ----
    figure('Color','w','Name',sprintf('FFT Contact vs Free (beta=%.3g)', beta));
    hold on; grid on;

    if ~isempty(Pc_mean)
        plot(f_c, Pc_mean, 'LineWidth', 1.6, 'DisplayName','Contact region');
    end
    if ~isempty(Pf_mean)
        plot(f_f, Pf_mean, '--', 'LineWidth', 1.6, 'DisplayName','Free region');
    end

    xlabel('Frequency (Hz)');
    ylabel('|U(f)| (average)');
    title(sprintf('Distributed DNC Beam: FFT Contact vs Free (\\beta=%.3g, \\kappa=%.3g)', ...
                  beta, kappa), 'Interpreter','tex');
    legend('Location','northeast');
    xlim([0 10]); ylim([0 inf]);
end

end

function S = beam_dnc_distributed_solver(beta, L,T,dx,dt, alpha2,kappa,psiFun,f0, phi,u0f,v0f)
% Minimal distributed DNC solver used for FFT analysis.

M = round(L/dx);
N = round(T/dt);

x = linspace(0,L,M+1).';
t = linspace(0,T,N+1);
psi = psiFun(x);

u = zeros(M+1, N+1);
gamma = alpha2*(dt^2)/(dx^4);

u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);

for j = 1:(N-1)
    I = (u(:, j+1) <= psi);
    I(1) = false;

    diag_add = dt^2*(kappa*double(I(2:end)) + (beta/dt)*double(I(2:end)));
    rhs_add  = dt^2*(kappa*psi(2:end).*double(I(2:end)) + ...
                     (beta/dt)*double(I(2:end)).*u(2:end,j+1));

    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % Row i=1
    A(1,1) = 1 + 7*gamma;
    if M>=2, A(1,2) = -4*gamma; end
    if M>=3, A(1,3) =  gamma; end

    b(1) = dt^2*f0 + 2*u(2,j+1) - u(2,j) + 4*gamma*phi(t(j+1));
    A(1,1) = A(1,1) + diag_add(1);
    b(1)   = b(1)   + rhs_add(1);

    % Interior rows
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r,i-2) = A(r,i-2) + gamma;
        end
        A(r,i-1) = A(r,i-1) - 4*gamma;
        A(r,i)   = A(r,i)   + (1 + 6*gamma);
        A(r,i+1) = A(r,i+1) - 4*gamma;
        A(r,i+2) = A(r,i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1,j+1) - u(i+1,j);

        A(r,i) = A(r,i) + diag_add(i);
        b(r)   = b(r)   + rhs_add(i);
    end

    % Right end i=M-1
    if M>=2
        A(M-1,M-3) = A(M-1,M-3) + gamma;
        A(M-1,M-2) = A(M-1,M-2) - 4*gamma;
        A(M-1,M-1) = A(M-1,M-1) + (1 + 5*gamma);
        A(M-1,M)   = A(M-1,M)   - 2*gamma;

        b(M-1) = dt^2*f0 + 2*u(M,j+1) - u(M,j);

        A(M-1,M-1) = A(M-1,M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % Right end i=M
    A(M,M-2) = A(M,M-2) + 2*gamma;
    A(M,M-1) = A(M,M-1) - 4*gamma;
    A(M,M)   = A(M,M)   + (1 + 2*gamma);

    b(M) = dt^2*f0 + 2*u(M+1,j+1) - u(M+1,j);

    A(M,M) = A(M,M) + diag_add(M);
    b(M)   = b(M)   + rhs_add(M);

    % Solve
    w = A \ b;

    u(1,j+2) = phi(t(j+2));
    u(2:end,j+2) = w;
end

% Pack
S.x = x;
S.t = t;
S.u = u;
S.psi = psi;
S.beta = beta;

end

function [f, meanP1] = fft_average_over_nodes(U, dt)
% Average single-sided FFT amplitude over rows of U.

if isempty(U)
    f = [];
    meanP1 = [];
    return;
end

[M, N] = size(U);
Fs = 1/dt;
K  = floor(N/2) + 1;

meanP1 = zeros(1, K);

for i = 1:M
    y = U(i,:);
    Y = fft(y);
    P2 = abs(Y/N);
    P1 = P2(1:K);
    if K>2
        P1(2:end-1) = 2*P1(2:end-1);
    end
    meanP1 = meanP1 + P1;
end

meanP1 = meanP1 / M;
f = Fs * (0:K-1) / N;

end
