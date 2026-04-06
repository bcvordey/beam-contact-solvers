% fft_vs_alpha2.m — Dominant frequency vs beam stiffness α^2 (clean FFT)
clc; close all;

% ---- Sets ----
alpha2Set = [4 6 8 10];     % beam stiffness values
kappa     = 1;              % fixed obstacle stiffness
beta      = 0.01;           % fixed damping
Fconst    = -0.2;           % forcing
yminus    = -0.02;          % obstacle level

% Solver/grid defaults
params = struct('T',5,'dt',0.01,'dx',0.001,'y_minus',yminus);

% FFT options
useFreeSegment = true;      % use pre-contact segment (cleaner, fewer impact harmonics)
padFactor      = 4;         % zero-padding factor for finer frequency grid

nA = numel(alpha2Set);
cols = lines(nA);
f_peak = nan(nA,1);

figure('Color','w','Name','Tip spectrum vs alpha2'); hold on; grid on;

for ia = 1:nA
    A2 = alpha2Set(ia);
    p  = params; p.alpha2 = A2;

    % Solve and get the right-tip signal
    [~, t, ~, met] = beam_contact_implicit_FD(Fconst, kappa, beta, p);
    tip = met.tip;                 % u(L,t)

    % Sampling
    dt = t(2) - t(1);
    Fs = 1/dt;

    % Choose segment: free vibration before first contact (if desired)
    if useFreeSegment && any(met.contact)
        j1  = find(met.contact, 1, 'first');
        idx = 1:max(j1-1, 64);     % ensure a minimal length
    else
        idx = 1:numel(t);
    end
    x  = detrend(tip(idx),'linear');           % remove DC and linear trend

    % Window + zero pad
    N   = numel(x);
    if exist('hann','file'), w = hann(N,'periodic'); else, n=(0:N-1)'; w=0.5-0.5*cos(2*pi*n/(N-1)); end
    xw  = x(:).*w;
    Nfft= 2^nextpow2(N*padFactor);
    X   = fft(xw, Nfft);

    % One-sided amplitude spectrum with proper scaling
    N2  = floor(Nfft/2)+1;
    f   = (0:N2-1)*(Fs/Nfft);                 % Hz
    A   = abs(X(1:N2)) / (sum(w)/2);          % window + one-sided amplitude scaling
    A(1)= 0;                                  % ignore DC for normalization

    % Normalize each curve to its own peak for readability
    Aplot = A / max(A + eps);
    plot(f, Aplot, 'LineWidth',1.6, 'Color', cols(ia,:));

    % Dominant frequency (exclude DC)
    [~, kmax] = max(A(2:end));  kmax = kmax+1;
    f_peak(ia) = f(kmax);

    fprintf('alpha^2 = %g: dominant freq ≈ %.3f Hz (segment %d samples)\n', A2, f_peak(ia), N);
end

xlabel('Frequency (Hz)');
ylabel('Normalized |U(f)|');
ylim([0,2])
ylim([1.1,1.8])
title(sprintf('Tip spectrum (clean segment) — \\kappa=%g, \\beta=%g, f=%g', kappa, beta, Fconst), 'Interpreter','tex');
legend(arrayfun(@(a)sprintf('\\alpha^2=%g',a), alpha2Set,'UniformOutput',false), 'Location','northeast');
hold off;
