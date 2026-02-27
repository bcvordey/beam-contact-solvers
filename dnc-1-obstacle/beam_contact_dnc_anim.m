function [x,t,u,pen,reaction,contact] = beam_contact_dnc_anim()
% BEAM_CONTACT_DNC_ANIM
% Implicit second-order FD for Euler–Bernoulli beam with unilateral (DNC) contact
% at the right end; left end clamped with u(0,t)=0 and u_x(0,t)=0.
% Adds: animation of u(x,t), tip displacement vs time, penetration depth, and reaction.
%
% Outputs:
%   x, t      : grids
%   u         : displacement field (size: (M+1) x (N+1))
%   pen       : penetration depth p(t) = max(0, y_- - u(L,t))
%   reaction  : normal reaction R(t) = kappa*p + beta*dp/dt (zero off-contact)
%   contact   : logical vector marking contact states

%% ---------------- Parameters ----------------
L = 1.0;            % beam length
T = 5.0;            % final time
dx = 0.01;          % space step
dt = 0.01;          % time step
M  = round(L/dx);   % spatial subintervals => nodes i = 0..M (MATLAB rows 1..M+1)
N  = round(T/dt);   % time steps            => j = 0..N  (MATLAB cols 1..N+1)

alpha2 = 1;     % alpha^2
f0     = -0.2;      % constant forcing
y_minus = -0.02;    % obstacle level
kappa  = 10;      % contact stiffness
beta   = 0.1;      % contact damping

phi = @(tt) 0.0;    % left displacement boundary u(0,t)=0
u0  = @(xx) 0.0;    % initial displacement
v0  = @(xx) 0.0;    % initial velocity

gamma = alpha2*(dt^2)/(dx^4);  % scheme parameter

%% ---------------- Grids & storage ----------------
x = linspace(0,L,M+1).';   % column (i = 0..M)
t = linspace(0,T,N+1);     % row    (j = 0..N)
u = zeros(M+1, N+1);       % u(i+1, j+1) ~ u_i^j

% Initial conditions: u^0, u^1
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);

% Enforce Dirichlet at x=0 for all times touched
u(1, :) = phi(t);

%% ---------------- Time stepping ----------------
for j = 1:(N-1)   % j=1..N-1 (we have u^0 and u^1; compute up to u^N)

    % Contact decision at right end uses tip from time level j (u_M^j)
    in_contact = (u(M+1, j+1) <= y_minus);

    % Unknown w := [u_1^{j+1}; ... ; u_M^{j+1}]  (size M)
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % ---- Row for i=1 (uses u_{-1}^{j+1} = u_1^{j+1} via u_x(0,t)=0) ----
    % (1 + 7γ) u_1 - 4γ u_2 + γ u_3 = dt^2 f + 2 u_1^j - u_1^{j-1} + 4γ φ(t_{j+1})
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(1+1, j+1) - u(1+1, j) + 4*gamma*phi(t(j+1));  % phi=0 here

    % ---- Interior: i = 2 .. M-2 (standard 5-pt stencil) ----
    % γ u_{i-2} - 4γ u_{i-1} + (1+6γ) u_i - 4γ u_{i+1} + γ u_{i+2} = dt^2 f + 2u_i^j - u_i^{j-1}
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            b(r) = b(r) - gamma*phi(t(j+1));  % zero here
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % ---- i = M-1 row ----
    % γ u_{M-3} - 4γ u_{M-2} + (1+5γ) u_{M-1} - 2γ u_M = dt^2 f + 2 u_{M-1}^j - u_{M-1}^{j-1}
    if M >= 3, A(M-1, M-3) = gamma;     end
    if M >= 2, A(M-1, M-2) = -4*gamma;  end
    A(M-1, M-1) = 1 + 5*gamma;
    A(M-1, M  ) = -2*gamma;
    b(M-1) = dt^2*f0 + 2*u(M-1+1, j+1) - u(M-1+1, j);

    % ---- i = M row: free vs contact ----
    if ~in_contact
        % No contact: u_xx(L,t)=0, u_xxx(L,t)=0
        % 2γ u_{M-2} - 4γ u_{M-1} + (1+2γ) u_M = dt^2 f + 2 u_M^j - u_M^{j-1}
        A(M, M-2) = 2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) = 1 + 2*gamma;
        b(M) = dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        % Contact (DNC): R = kappa*p + beta*dp/dt, p = max(0, y_- - u(L,t))
        c = (dx^3/alpha2)*(kappa + beta/dt);
        A(M, M-2) = 2*gamma;
        A(M, M-1) = -4*gamma;
        A(M, M  ) = 1 + 2*gamma*(c + 1);
        b(M) = dt^2*f0 ...
             + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             - u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    % ---- Solve for w and write back ----
    w = A \ b;
    u(1,   j+2) = phi(t(j+2));  % left Dirichlet
    u(2:end,j+2) = w;           % unknowns
end

%% ---------------- Post-processing: contact, penetration, reaction ----------------
tip = u(end,:);                 % u(L,t)
gap = tip - y_minus;            % >0 separation, <0 raw interpenetration
pen = max(0, -gap);             % p(t) = max(0, y_- - u(L,t))
pen_rate = [0, diff(pen)/dt];   % backward diff, same length
contact = pen > 0;              % logical
reaction = kappa*pen + beta*pen_rate;
reaction(~contact) = 0;         % zero off-contact for clarity

%% ---------------- Animation ----------------
%% ---------------- Animation (with pacing controls) ----------------
fig1 = figure('Name','Beam with DNC: animation','Color','w');
ymax = 1.05*max([1e-12, max(abs(u),[],'all'), abs(y_minus)]);
plt = plot(x, u(:,1), 'b-', 'LineWidth', 3); hold on; grid on;
line([1, 1], [-1.1*ymax, y_minus], 'Color', 'k', 'LineWidth', 5); % Vertical bar
% yline(y_minus, 'k--', 'Obstacle (y_-)', 'LabelHorizontalAlignment','left');
tipMarker = plot(x(end), u(end,1), 'go', 'MarkerFaceColor','g', 'DisplayName','Tip');
xlabel('x'); ylabel('u(x,t)');
title('Vibration with Damped Normal Compliance at x=L');
axis([0, x(end), -ymax, ymax]);

% --- Animation controls ---
skip = 1;                 % 1 = show every step; increase to skip steps
playback_speed = 10;     % sim-seconds per real-second:
                          % 1.0 = real-time, 0.5 = 2x slower (longer), 0.25 = 4x slower
% Alternatively, target a fixed FPS instead of a sim/real mapping:
% target_fps = 25; frame_real_dt = 1/target_fps;

for j = 1:skip:numel(t)
    set(plt, 'YData', u(:,j));
    set(tipMarker, 'YData', u(end,j));
    if u(end,j) < y_minus   % in contact => red tip
        set(tipMarker, 'MarkerEdgeColor','r','MarkerFaceColor','r');
    else
        set(tipMarker, 'MarkerEdgeColor','g','MarkerFaceColor','g');
    end
    drawnow;
    % Pace the playback:
    pause( (skip*dt)/max(playback_speed, eps) );   % comment this line if you prefer max FPS
    % or: pause(frame_real_dt);                    % use fixed FPS pacing instead
end


%% ---------------- Plots: tip trace, penetration, reaction ----------------
% Helper: make a tag for filenames using kappa,beta (replace '.' with 'p')
fmtVal = @(v) regexprep(regexprep(sprintf('%.6f', v),'0+$',''),'\.$','');
tag = sprintf('kappa-%s_beta-%s', strrep(fmtVal(kappa),'.','p'), strrep(fmtVal(beta),'.','p'));

% Ensure output folder
outDir = fullfile(pwd, 'plots');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% % Tip vs time with obstacle and contact markers
fig2 = figure('Name','Tip vs time & obstacle','Color','w'); hold on; grid on;
plot(t, tip, 'b-', 'LineWidth', 1.6, 'DisplayName','u(L,t)');
yline(y_minus, 'k--', 'DisplayName','Obstacle y_-');
plot(t(contact), tip(contact), 'ro', 'MarkerSize', 3, 'DisplayName','Contact');
xlabel('time'); ylabel('displacement at tip'); title('Right-end displacement vs time (DNC)');
legend('Location','best');
% save_pdf(gcf, fullfile(outDir, sprintf('tip_vs_time__%s.pdf', tag)));

% Penetration depth
% fig3 = figure('Name','Penetration depth','Color','w'); grid on;
% plot(t, pen, 'm-', 'LineWidth', 1.6);
% xlabel('time'); ylabel('penetration p(t) = max(0, y_- - u(L,t))');
% title('Penetration into obstacle (DNC)');
% % save_pdf(gcf, fullfile(outDir, sprintf('penetration__%s.pdf', tag)));

% Contact reaction
% fig4 = figure('Name','Contact reaction','Color','w'); grid on;
% plot(t, reaction, 'c-', 'LineWidth', 1.6);
% xlabel('time'); ylabel('$R(t) = \kappa\,p(t) + \beta\,\dot p(t)$','Interpreter','latex');
% title('Normal reaction at the tip (DNC)');
% % save_pdf(gcf, fullfile(outDir, sprintf('reaction__%s.pdf', tag)));

end

% --------- Helper: save a figure as vector PDF (robust across MATLAB versions) ---------
function save_pdf(figHandle, filename)
    try
        % R2020a+: exportgraphics keeps vector content when ContentType='vector'
        exportgraphics(figHandle, filename, 'ContentType','vector');
    catch
        % Fallback for older MATLAB: painters ensures vector for line/patch
        set(figHandle, 'Renderer', 'painters');
        print(figHandle, filename, '-dpdf');
    end
end
