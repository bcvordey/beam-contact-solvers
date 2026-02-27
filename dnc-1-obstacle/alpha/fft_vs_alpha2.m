% fft_vs_alpha2.m — Show that dominant frequency increases with beam stiffness α^2
clc; close all;

% ---- Choose sets ----
alpha2Set = [4 6 8 10];  % beam stiffness proxy
kappa     = 1;               % fix obstacle stiffness
beta      = 0.01;               % fix damping
Fconst    = -0.2;             % fix force (sign as in your model)
yminus    = -0.02;             % obstacle level (just for reference)

% ---- Solver grid defaults (override as needed) ----
params = struct('T',5,'dt',0.01,'dx',0.001,'y_minus',yminus);

% ---- Storage ----
nA = numel(alpha2Set);
specF   = cell(nA,1);   % frequency vectors
specA   = cell(nA,1);   % amplitude spectra
f_peak  = nan(nA,1);    % dominant frequencies
labelA2 = strings(nA,1);
cols    = lines(nA);

for ia = 1:nA
    A2 = alpha2Set(ia);
    p = params; p.alpha2 = A2;

    % Run your solver (tip is u(L,t) = u(end,:))
    [~, t, u, met] = beam_contact_implicit_FD(Fconst, kappa, beta, p);
    tip = met.tip;    % same as u(end,:)

    % Sampling info
    dt = t(2) - t(1);
    Fs = 1/dt;

    % ---- Pick FREE segment (before first contact) to avoid impact harmonics ----
    if any(met.contact)
        j1 = find(met.contact, 1, 'first');
        idx = 1:max(j1-1, 16);  % ensure at least a few samples
    else
        idx = 1:numel(t);       % no contact: use entire signal
    end
    sig = tip(idx);

    % Detrend + Hann window to reduce leakage
    sig = detrend(sig, 'linear');
    N   = numel(sig);
    w   = hann(N);
    xw  = sig(:) .* w;

    % FFT and one-sided amplitude spectrum
    X   = fft(xw);
    N2  = floor(N/2) + 1;
    f   = (0:N2-1)*(Fs/N);
    A   = abs(X(1:N2));
    % Amplitude normalization (window + one-sided)
    A   = A / (sum(w)/2);

    % Store
    specF{ia} = f;
    specA{ia} = A;
    labelA2(ia) = sprintf('\\alpha^2=%g', A2);

    % Dominant frequency (ignore DC bin)
    if N2 >= 2
        [~, kmax] = max(A(2:end));   % exclude f=0
        f_peak(ia) = f(kmax+1);
    else
        f_peak(ia) = NaN;
    end
end

% % ---------------- Figure 1: Spectra overlay (colors by α^2) ----------------
% f1 = figure('Color','w','Name','Tip_spectra_vs_alpha2');
% hold on; grid on;
% for ia = 1:nA
%     plot(specF{ia}, specA{ia}, 'LineWidth', 1.3, 'Color', cols(ia,:));
% end
% xlabel('Frequency (Hz)');
% ylabel('|U(f)| (a.u.)');
% title(sprintf('Tip Spectrum (free segment) — \\kappa=%g, \\beta=%g, F=%g', kappa, beta, Fconst), ...
%       'Interpreter','tex');
% ylim([0, 0.0033]);
% xlim([0, 5]);
% legend(cellstr(labelA2), 'Location','northeast');
% % saveas(f1, 'Tip_spectrum_vs_alpha2.jpg');
% 
% % ---------------- Figure 2: Peak frequency vs α^2 ----------------
% f2 = figure('Color','w','Name','Peak_freq_vs_alpha2');
% plot(alpha2Set, f_peak, '-o', 'LineWidth', 1.4);
% grid on;
% xlabel('\alpha^2');
% ylabel('Dominant frequency (Hz)');
% title(sprintf('Dominant Frequency vs. Beam Stiffness \\alpha^2 (\\kappa=%g, \\beta=%g, F=%g)', ...
%       kappa, beta, Fconst),'Interpreter','tex');

%---------- FFT ----------%
uFFT = fft(tip);

fig_fft = figure;
plot(abs(uFFT),'LineWidth',1.6,'MarkerSize', 1.5);
xlabel('Frequency (Hz)');
xlim([-10, 505]);
ylim([0, 1.3]);
ylabel('|FFT[u(M,t)]|');
% title({'FFT of Contact End of the Beam', sprintf('\\alpha=%.2f, \\beta=%.2f, \\kappa=%.3f, obstacle=%.2f, force=%.2f', alpha, beta, kappa, obstacle, force)});

%---------- SHIFT FFT ----------%
uFFTshift = fftshift(uFFT);

fig_fftshift = figure;
plot(abs(uFFTshift),'LineWidth', 1,'MarkerSize', 1.5);
% title({'FFT SHIFT of Contact End of the Beam', sprintf('\\alpha=%.2f, \\beta=%.2f, \\kappa=%.3f, obstacle=%.2f, force=%.2f', alpha, beta, kappa, obstacle, force)});
xlabel('Frequency (Hz)');
ylabel('|FFTshift[u(M,t)]|');