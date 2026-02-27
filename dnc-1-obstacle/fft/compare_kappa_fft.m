function compare_kappa_fft()
% Plain FFT of u at x=1 for kappa = [1, 10, 100, 1000]
% No detrend, no window. Just fft + single-sided amplitude spectrum.

kappas = [0.01, 1, 10, 100, 1000];

figure('Name','FFT at x=1 for various'); hold on; grid off;

for k = 1:numel(kappas)
    kappa = kappas(k);
    [beta, u, dt] = simulate_with_kappa(kappa);   % run your solver with this kappa

    y = u(end,:).';               % u at x=1 (last node), column vector
    N = numel(y);                 % number of time samples
    Fs = 1/dt;                    % sampling frequency (Hz)
    
    % --- plain FFT, single-sided amplitude spectrum
    Y  = fft(y);
    P2 = abs(Y/N);                % two-sided amplitude
    P1 = P2(1:floor(N/2)+1);      % single-sided
    if numel(P1) > 2
        P1(2:end-1) = 2*P1(2:end-1);
    end
    f = Fs*(0:floor(N/2))/N;      % frequency axis (Hz)

    plot(f, P1, 'DisplayName', sprintf('\\kappa = %g', kappa),'LineWidth', 1.0);
end

xlabel('$Frequency \ (Hz)$','Interpreter','latex');
ylabel('$Amplitude \ |U(f)|$','Interpreter','latex');
xlim([0, 5]);  % adjust if you want to zoom in
ylim([0, 0.014]);
% title('Single-Sided Amplitude Spectrum at x=1');
title(sprintf('Single-Sided Amplitude Spectrum at $x = 1$, $\\beta = %.3g$', beta), ...
      'Interpreter','latex');

legend show;
 
end

% ==========================================================
% Minimal wrapper of your solver with kappa as an argument
% (based directly on your posted function, but parameterized)
% ==========================================================
function [beta,u,dt] = simulate_with_kappa(kappa)
    % --- your original parameters (kept the same) ---
    L = 1.0;             % beam length
    T = 20.0;            % final time
    dx = 0.01;           % space step
    dt = 0.01;           % time step
    M  = round(L/dx);    % nodes 0..M
    N  = round(T/dt);    % times 0..N

    alpha2 = 1;     
    f0     = -0.2;       
    y_minus = -0.02;     
    beta   = 0.001;        

    phi = @(tt) 0.0;     
    u0  = @(xx) 0.0;     
    v0  = @(xx) 0.0;     

    gamma = alpha2*(dt^2)/(dx^4);

    % grids and storage
    x = linspace(0,L,M+1).';   
    t = linspace(0,T,N+1);     
    u = zeros(M+1, N+1);       

    % initial conditions
    u(:,1) = u0(x);
    u(:,2) = u0(x) + dt*v0(x);
    u(1, :) = phi(t);

    % time stepping
    for j = 1:(N-1)
        in_contact = (u(M+1, j+1) <= y_minus);

        A = spalloc(M, M, 5*M);
        b = zeros(M,1);

        % i=1
        A(1,1) = 1 + 7*gamma;
        if M >= 2, A(1,2) = -4*gamma; end
        if M >= 3, A(1,3) =  gamma;   end
        b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));

        % interior i=2..M-2
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
        end

        % i = M-1
        if M >= 2
            A(M-1, M-3) = gamma;
            A(M-1, M-2) = -4*gamma;
            A(M-1, M-1) = 1 + 5*gamma;
            A(M-1, M  ) = -2*gamma;
            b(M-1) = dt^2*f0 + 2*u(M, j+1) - u(M, j);
        end

        % i = M (right boundary)
        if ~in_contact
            A(M, M-2) = 2*gamma;
            A(M, M-1) = -4*gamma;
            A(M, M  ) = 1 + 2*gamma;
            b(M) = dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
        else
            c = (dx^3/alpha2)*(kappa + beta/dt);
            A(M, M-2) = 2*gamma;
            A(M, M-1) = -4*gamma;
            A(M, M  ) = 1 + 2*gamma*(c + 1);
            b(M) = dt^2*f0 ...
                 + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
                 - u(M+1, j) ...
                 + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
        end

        w = A \ b;
        u(1,   j+2) = phi(t(j+2));
        u(2:end,j+2) = w;
    end
end
