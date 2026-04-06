function [x,t,u,pen,reaction,contact] = beam_contact_dnc_gif()
% BEAM_CONTACT_DNC_ANIM
% Implicit FD Euler–Bernoulli beam with DNC at x=L, animation + exports:
% - Animated GIF (and optional MP4) of the beam motion
% - PDFs of tip trace, penetration, and reaction, tagged with kappa/beta

%% ---------------- Parameters ----------------
L = 1.0;   T = 10.0;   dx = 0.01;   dt = 0.01;
M = round(L/dx);   N = round(T/dt);

alpha2 = 1;
f0 = -5;
y_minus = -0.02;
kappa  = 1000;
beta   = 0.1;

phi = @(tt) 0.0;   u0 = @(xx) 0.0;   v0 = @(xx) 0.0;
gamma = alpha2*(dt^2)/(dx^4);

%% ---------------- Grids & storage ----------------
x = linspace(0,L,M+1).';
t = linspace(0,T,N+1);
u = zeros(M+1, N+1);
u(:,1) = u0(x);
u(:,2) = u0(x) + dt*v0(x);
u(1,:) = phi(t);

%% ---------------- Time stepping ----------------
for j = 1:(N-1)
    in_contact = (u(M+1, j+1) <= y_minus);
    A = spalloc(M, M, 5*M); b = zeros(M,1);

    % i = 1 row (left end with u_x=0 via ghost)
    A(1,1) = 1 + 7*gamma;  if M>=2, A(1,2) = -4*gamma; end
    if M>=3, A(1,3) = gamma; end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j);

    % interior rows i = 2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1, A(r,i-2) = A(r,i-2) + gamma; end
        A(r,i-1) = A(r,i-1) - 4*gamma;
        A(r,i  ) = A(r,i  ) + (1 + 6*gamma);
        A(r,i+1) = A(r,i+1) - 4*gamma;
        A(r,i+2) = A(r,i+2) + gamma;
        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);
    end

    % i = M-1
    if M>=3, A(M-1,M-3) = gamma; end
    if M>=2, A(M-1,M-2) = -4*gamma; end
    A(M-1,M-1) = 1 + 5*gamma;  A(M-1,M) = -2*gamma;
    b(M-1) = dt^2*f0 + 2*u(M, j+1) - u(M, j);

    % i = M (right end: free vs contact)
    if ~in_contact
        A(M,M-2) = 2*gamma;  A(M,M-1) = -4*gamma;  A(M,M) = 1 + 2*gamma;
        b(M) = dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    else
        c = (dx^3/alpha2)*(kappa + beta/dt);
        A(M,M-2) = 2*gamma;  A(M,M-1) = -4*gamma;  A(M,M) = 1 + 2*gamma*(c+1);
        b(M) = dt^2*f0 ...
             + 2*(1 + gamma*(beta*dx^3)/(alpha2*dt))*u(M+1, j+1) ...
             - u(M+1, j) ...
             + 2*gamma*(dx^3/alpha2)*kappa*y_minus;
    end

    w = A\b;
    u(1, j+2)   = phi(t(j+2));
    u(2:end,j+2) = w;
end

%% ---------------- Post-processing ----------------
tip      = u(end,:);
gap      = tip - y_minus;
pen      = max(0, -gap);
pen_rate = [0, diff(pen)/dt];
contact  = pen > 0;
reaction = kappa*pen + beta*pen_rate;  reaction(~contact) = 0;

% Tag + output dir for files
fmtVal = @(v) regexprep(regexprep(sprintf('%.6f', v),'0+$',''),'\.$','');
tag = sprintf('kappa-%s_beta-%s', strrep(fmtVal(kappa),'.','p'), strrep(fmtVal(beta),'.','p'));
outDir = fullfile(pwd, 'plots');  if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---------------- Animation + EXPORT (GIF/MP4) ----------------
fig1 = figure('Name','Beam with DNC: animation','Color','w');
ymax = 1.05*max([1e-12, max(abs(u),[],'all'), abs(y_minus)]);
plt = plot(x, u(:,1), 'b-', 'LineWidth', 3); hold on; grid on;
line([1, 1], [-1.1*ymax, y_minus], 'Color', 'k', 'LineWidth', 5); % Vertical bar
tipMarker = plot(x(end), u(end,1), 'go', 'MarkerFaceColor','g', 'DisplayName','Tip');
xlabel('x'); ylabel('u(x,t)'); title(sprintf('One Boundary Obstacle at x=1 ($\\kappa = %.3g, \\beta = %.3g$)',kappa,beta), 'Interpreter','latex');
axis([0, x(end), -ymax, ymax]);

% --- playback controls ---
skip = 5;                 % show every step (increase to skip frames)
playback_speed = 1000;     % sim-sec per real-sec (1.0=real time, 0.5 = 2x slower)
delayTime = (skip*dt)/max(playback_speed, eps);    % seconds per frame

% --- export controls ---
writeGif = true;   gifFile = fullfile(outDir, sprintf('beam_dnc_anim__%s.gif', tag));
writeMp4 = false;  mp4File = fullfile(outDir, sprintf('beam_dnc_anim__%s.mp4', tag));

if writeMp4
    v = VideoWriter(mp4File, 'MPEG-4');
    v.FrameRate = max(1, min(60, round(1/delayTime)));
    open(v);
end
if writeGif
    firstGIF = true;
end

for j = 1:skip:numel(t)
    set(plt, 'YData', u(:,j));
    set(tipMarker, 'YData', u(end,j));
    if contact(j)
        set(tipMarker,'MarkerEdgeColor','r','MarkerFaceColor','r');
    else
        set(tipMarker,'MarkerEdgeColor','g','MarkerFaceColor','g');
    end
    drawnow;

    % Capture frame once, feed to MP4 and GIF
    F = getframe(fig1);
    if writeMp4, writeVideo(v, F); end
    if writeGif
        [im, cm] = rgb2ind(frame2im(F), 256);
        if firstGIF
            imwrite(im, cm, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
            firstGIF = false;
        else
            imwrite(im, cm, gifFile, 'gif', 'WriteMode','append', 'DelayTime', delayTime);
        end
    end

    % Optional real-time pacing in the live preview:
    pause(delayTime);
end
if writeMp4, close(v); end

%% ---------------- Time-series plots (saved as PDFs, tagged) ----------------
% % Tip vs time
% fig2 = figure('Name','Tip vs time & obstacle','Color','w'); hold on; grid on;
% plot(t, tip, 'b-', 'LineWidth', 1.5, 'DisplayName','u(L,t)');
% yline(y_minus, 'k--', 'DisplayName','Obstacle y_-');
% plot(t(contact), tip(contact), 'ro', 'MarkerSize', 3, 'DisplayName','Contact');
% xlabel('time'); ylabel('displacement at tip'); title('Right-end displacement vs time (DNC)');
% legend('Location','best');
% save_pdf(gcf, fullfile(outDir, sprintf('tip_vs_time__%s.pdf', tag)));
% 
% % Penetration
% fig3 = figure('Name','Penetration depth','Color','w'); grid on;
% plot(t, pen, 'm-', 'LineWidth', 1.6);
% xlabel('time'); ylabel('penetration p(t) = max(0, y_- - u(L,t))');
% title('Penetration into obstacle (DNC)');
% save_pdf(gcf, fullfile(outDir, sprintf('penetration__%s.pdf', tag)));
% 
% % Reaction
% fig4 = figure('Name','Contact reaction','Color','w'); grid on;
% plot(t, reaction, 'c-', 'LineWidth', 1.5);
% xlabel('time'); ylabel('R(t) = \kappa\,p(t) + \beta\,\dot p(t)');
% title('Normal reaction at the tip (DNC)');
% save_pdf(gcf, fullfile(outDir, sprintf('reaction__%s.pdf', tag)));
end

% ---------- helper: save a figure as vector PDF ----------
function save_pdf(figHandle, filename)
    try
        exportgraphics(figHandle, filename, 'ContentType','vector');
    catch
        set(figHandle, 'Renderer', 'painters');
        print(figHandle, filename, '-dpdf');
    end
end
