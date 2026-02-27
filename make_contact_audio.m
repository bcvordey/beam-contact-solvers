% make_contact_audio.m — listen to the beam contacts for soft/stiff/hard obstacles
clc; close all;

% ---------- Simulation settings ----------
alpha2   = 6;                 % beam stiffness
kappaSet = [0.01, 100, 1000]; % soft / stiff / hard obstacle
beta     = 0.1;               % damping
Fconst   = -0.25;             % external force (or set 0 to hear pure impacts)
yminus   = -0.02;             % obstacle level
params   = struct('T',5,'dt',0.002,'dx',0.01,'y_minus',yminus);  % Fs_sim = 1/dt = 500 Hz

% Optional boundary excitation example:
% Aexc = 0.001; fHz = 2;
% params.phi_t = @(tt) Aexc*sin(2*pi*fHz*tt);

% ---------- Audio mapping ----------
speedUp   = 60;                % play 'speedUp' times faster -> all freqs × speedUp
hpCutHz   = 40;                % remove DC/very low rumble before playback
doWrite   = true;              % save WAVs
outFs     = [];                % [] = just use Fs_play; or set 44100 to hard-set

for K = kappaSet
    % Run the simulation
    [~, t, ~, met] = beam_contact_implicit_FD(Fconst, K, beta, params);
    tip = met.tip;                     % u(L,t)
    dt  = t(2)-t(1);
    Fs_sim = 1/dt;

    % ---- Build a contact-force proxy (spring + dashpot only when in contact) ----
    d      = max(0, yminus - tip);     % penetration
    ddot   = [diff(d)/dt, 0];          % time-derivative (simple finite diff)
    Fcont  = K.*d + beta.*ddot;        % contact force proxy
    Fcont(~(d>0)) = 0;                 % only when in contact

    % Choose which to listen to: tip or contact force
    xsig = Fcont;                      % try tip instead: xsig = tip;

    % ---- Audio pre-processing ----
    % Time-compress by playing at a higher sample rate:
    Fs_play = Fs_sim * speedUp;        % pitch & speed × speedUp (no resampling needed)

    % High-pass to remove DC/low drift (keep transients/clicks)
    xhp = highpass(xsig, hpCutHz/speedUp, Fs_sim);  % compensate cutoff by speedUp

    % Normalize safely for playback
    xhp = xhp - mean(xhp);
    xhp = xhp / max(1e-12, max(abs(xhp)));

    % Optional: resample to a standard audio rate (e.g., 44100 Hz)
    if ~isempty(outFs)
        [p,q] = rat(outFs / Fs_play, 1e-9);
        y = resample(xhp, p, q);  Fs_out = outFs;
    else
        y = xhp;                  Fs_out = Fs_play;
    end

    % ---- Play and/or save ----
    fprintf('Playing K=%.3g  at Fs_play = %.1f Hz (speedUp=%gx)\n', K, Fs_out, speedUp);
    soundsc(y, Fs_out);   % scales to [-1,1] and plays loudly but unclipped

    if doWrite
        fname = sprintf('contact_audio_kappa_%g_speedup_%gx.wav', K, speedUp);
        audiowrite(fname, y, Fs_out);
        fprintf('Wrote %s\n', fname);
    end
end
