% dis-obs.m
% Solve a distributed-obstacle beam problem and save a GIF animation.
% Running this file creates one figure (the animation) and writes one GIF.

%% ---- Parameters (edit as needed) ----
L  = 1.0;
T  = 10.0;
dx = 0.01;
dt = 0.01;

alpha2  = 1;      % alpha^2
kappa   = 100;    % obstacle stiffness
f0      = -0.2;   % constant forcing f(x,t)=f0
betaVec = 0.001;  % one or more beta values

% Obstacle profile psi(x)
psiFun  = @(x) -0.02 + 0*x;

phi  = @(tt) 0.0;  % boundary excitation at x=0
u0f  = @(xx) 0.0;  % initial displacement
v0f  = @(xx) 0.0;  % initial velocity

%% ---- Solve for each beta ----
nb  = numel(betaVec);
sol = cell(nb,1);
for k = 1:nb
    sol{k} = solve_single_beta_psi(betaVec(k), L, T, dx, dt, ...
        alpha2, kappa, psiFun, f0, phi, u0f, v0f);
end

%% ---- GIF animation (single figure only) ----
kAnim = 1;  % which beta solution to animate (1..nb)
numtag = @(v) strrep(strrep(sprintf('%.3g', v), '.', 'p'), '-', 'm');
gif_filename = sprintf('beam_dist_anim_k%s_b%s.gif', ...
    numtag(kappa), numtag(betaVec(kAnim)));

animate_beam_distributed_gif(sol{kAnim}, kappa, gif_filename);
fprintf('GIF saved to: %s\n', fullfile(pwd, gif_filename));

% ------------ Helper: solve for a single beta with psi(x) ------------
function S = solve_single_beta_psi(beta, L, T, dx, dt, alpha2, kappa, psiFun, f0, phi, u0f, v0f)
M = round(L/dx);
N = round(T/dt);
x = linspace(0, L, M+1).';
t = linspace(0, T, N+1);
psi = psiFun(x);
u = zeros(M+1, N+1);
gamma = alpha2*(dt^2)/(dx^4);

% Initial data (2nd-order start)
u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);  % Dirichlet

for j = 1:(N-1)
    I = (u(:, j+1) <= psi);
    I(1) = false;
    diag_add = dt^2 * (kappa * double(I(2:end)) + (beta/dt) * double(I(2:end)));
    rhs_add  = dt^2 * (kappa * psi(2:end) .* double(I(2:end)) ...
                     + (beta/dt) * double(I(2:end)) .* u(2:end, j+1));

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

    % Right end i=M-1 (free end)
    if M >= 2
        A(M-1, M-3) = A(M-1, M-3) + gamma;
        A(M-1, M-2) = A(M-1, M-2) - 4*gamma;
        A(M-1, M-1) = A(M-1, M-1) + (1 + 5*gamma);
        A(M-1, M  ) = A(M-1, M  ) - 2*gamma;
        b(M-1) = b(M-1) + dt^2*f0 + 2*u(M, j+1) - u(M, j);
        A(M-1, M-1) = A(M-1, M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % Right end i=M (free end)
    A(M, M-2) = A(M, M-2) + 2*gamma;
    A(M, M-1) = A(M, M-1) - 4*gamma;
    A(M, M  ) = A(M, M  ) + (1 + 2*gamma);
    b(M) = b(M) + dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    A(M, M) = A(M, M) + diag_add(M);
    b(M)    = b(M)    + rhs_add(M);

    % Solve and write back
    w = A \ b;
    u(1, j+2) = phi(t(j+2));
    u(2:end, j+2) = w;
end

% Derived fields
contact = (u <= psi);
penetration = max(psi - u, 0);
frac_contact = mean(contact,1);
max_pen = max(penetration, [], 1);

% Pack
S.x = x;
S.t = t;
S.u = u;
S.psi = psi;
S.contact = contact;
S.penetration = penetration;
S.frac_contact = frac_contact;
S.max_pen = max_pen;
S.beta = beta;
end

function animate_beam_distributed_gif(S, kappa, gif_filename)
% Single-figure animation with automatic GIF export.
% Tip marker is green when free and red when penetrating/contacting obstacle.

x   = S.x;
t   = S.t;
u   = S.u;
psi = S.psi;
L   = x(end);
dt  = t(2) - t(1);

fig = figure('Name', sprintf('Beam with DNC: animation (kappa = %g, beta = %g)', kappa, S.beta), ...
             'Color', 'w');
ax = axes('Parent', fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

plot(ax, x, psi, 'k--', 'LineWidth', 1.4, 'DisplayName', 'Obstacle \psi(x)');
hBeam = plot(ax, x, u(:,1), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Beam');
hPen = plot(ax, nan, nan, 'r.', 'MarkerSize', 10, 'DisplayName', 'Penetration');
hTip = plot(ax, x(end), u(end,1), 'go', 'MarkerFaceColor', 'g', ...
            'MarkerSize', 7, 'DisplayName', 'Tip');
plot(ax, 0, 0, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'DisplayName', 'Clamp');

xlabel(ax, 'x');
ylabel(ax, 'u(x,t)');
ymin = min([u(:); psi(:)]);
ymax = max([u(:); psi(:)]);
ypad = 0.05 * max(1e-12, ymax - ymin);
axis(ax, [0, L, ymin - ypad, ymax + ypad]);

hTitle = title(ax, sprintf('Distributed obstacle at t = %.2f  (kappa = %.3g, beta = %.3g)', ...
    t(1), kappa, S.beta));
legend(ax, 'Location', 'best');

% Animation/GIF controls
skip = 3;               % show every skip-th frame
playback_speed = 1000;  % sim-seconds per real-second
frame_delay = (skip * dt) / max(playback_speed, eps);
firstFrame = true;

for j = 1:skip:numel(t)
    uj = u(:,j);
    inContact = (uj <= psi);

    set(hBeam, 'YData', uj);
    set(hPen, 'XData', x(inContact), 'YData', uj(inContact));
    set(hTip, 'YData', uj(end));

    if inContact(end)
        set(hTip, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    else
        set(hTip, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
    end

    set(hTitle, 'String', sprintf('Distributed obstacle at t = %.2f  (kappa = %.3g, beta = %.3g)', ...
        t(j), kappa, S.beta));
    drawnow;

    frame = getframe(fig);
    [imind, cmap] = rgb2ind(frame2im(frame), 256);
    if firstFrame
        imwrite(imind, cmap, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frame_delay);
        firstFrame = false;
    else
        imwrite(imind, cmap, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end
end
