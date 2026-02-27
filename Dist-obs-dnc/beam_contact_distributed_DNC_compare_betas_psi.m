function beam_contact_distributed_DNC_compare_betas_psi()


%% ---- Parameters (edit as needed) ----
L  = 1.0;   T  = 10.0;
dx = 0.01;  dt = 0.01;

alpha2  = 1;            % α^2
kappa   = 100;           % κ
f0      = -0.2;           % constant forcing f(x,t)=f0

betaVec = 0.001;   % Define the range of beta values for comparison

% psiFun = @(x) -0.02 + 0.01*sin(pi*x).^2;   % range: [-0.02, -0.01]
% psiFun  = @(x) -0.02 - 0.01*sin(pi*x).^2;                 % sinusoidal trench
  % psiFun  = @(x) -(0.02 + 0.01*cos(50*pi*x));                % wavy floor
psiFun  = @(x) -0.02 + 0*x;                               % back to constant

phi  = @(tt) 0.0;         % boundary excitation at x=0
u0f  = @(xx) 0.0;         % initial displacement
v0f  = @(xx) 0.0;         % initial velocity

% Solve for each beta and store results
nb   = numel(betaVec);
sol  = cell(nb,1);
for k = 1:nb
    sol{k} = solve_single_beta_psi(betaVec(k), L,T,dx,dt, alpha2,kappa,psiFun,f0, phi,u0f,v0f);
end

% Shared grids (identical across beta runs)
x = sol{1}.x;   t = sol{1}.t;   psi_vec = sol{1}.psi;  psi_tip = psi_vec(end);
% ===== Minimal 2D animation (choose which beta to animate) =====
     kAnim = 1;           % 1..nb
     make_video = true;  % set true to save MP4
     animate_beam_distributed_2d(sol{kAnim}, kappa, make_video, ...
     sprintf('beam_dist_anim_k%.3g_b%.3g.mp4', kappa, betaVec(kAnim)));


cols = lines(nb);

%% ====== Overlaid line plots ======
numtag = @(x) strrep(strrep(sprintf('%.3g', x), '.', 'p'), '-', 'm');
kappa_tag = numtag(kappa);
beta_tags = arrayfun(numtag, betaVec, 'UniformOutput', false);
beta_tag  = strjoin(beta_tags, '_');





% 1) Tip displacement vs time
% fig1 = figure('Color','w','Name','Tip displacement vs time'); hold on;
% for k = 1:nb
%     plot(t, sol{k}.u(end,:), 'LineWidth', 1.8, 'Color', cols(k,:));
% end
% yline(psi_tip,'--','LineWidth',1.4,'DisplayName','\psi(1)');
% xlabel('Time','Interpreter','latex'); ylabel('$u(1,t)$','Interpreter','latex');
% title('Right-Tip displacement','Interpreter','latex');
% lg = [arrayfun(@(b)sprintf('$\\beta=%.3g$',b), betaVec, 'uni',0), {'$\psi(1)$'}];
% legend(lg, 'Interpreter','latex', 'Location','best');grid off;
% prettify(gca);

%% ====== Tiled image/surface plots (one panel per beta) ======
% 2) Contact maps (white=no, black=yes)
% fig2 = figure('Color','w','Name','Contact map');
[nrow,ncol] = tilesize(nb);
tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
% for k = 1:nb
%     nexttile;
%     imagesc(t, x, sol{k}.contact); axis xy tight;
%     colormap(gca, flipud(gray)); caxis([0 1]);
%     title(sprintf('$Contact: \\beta = %.3g$', betaVec(k)), 'Interpreter','latex');grid off;
%     xlabel('Time','Interpreter','latex'); ylabel('x','Interpreter','latex'); prettify(gca);
% end
% cb = colorbar; cb.Layout.Tile = 'east'; cb.Ticks=[0 1]; cb.TickLabels={'No','Contact'};

% 3) Penetration heatmaps
% fig3 = figure('Color','w','Name','Penetration heatmap');
% tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
% for k = 1:nb
%     nexttile;
%     imagesc(t, x, sol{k}.penetration); axis xy tight;
%     colormap(gca, parula);
%     title(sprintf('$Penetration: (\\beta = %.3g, \\kappa = %.3g)$', betaVec(k),kappa), 'Interpreter','latex');
%     xlabel('$Time$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); grid off; prettify(gca);
% end
% cb = colorbar; cb.Layout.Tile = 'east';
% ylabel(cb,'penetration $(u - \psi(x))_-$','Interpreter','latex');
 
% fname3 = sprintf('penetration_heatmaps_k%s_b%s.pdf', kappa_tag, beta_tag);
% exportgraphics(fig3, fname3, 'Resolution', 600, 'BackgroundColor','white');

% 4) 3D surfaces with obstacle ψ(x)
% fig4 = figure('Color','w','Name','Displacement surface');
% tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
% for k = 1:nb
%     nexttile;
%     surf(t, x, sol{k}.u, 'EdgeColor','none'); hold on;
%     surf(t, x, repmat(psi_vec,1,numel(t)), ...
%      'FaceColor','k', 'EdgeColor','r', 'FaceAlpha',0.5);
%     shading interp; view(45,30); axis tight vis3d; grid off;
%     title(sprintf('Displacement over $\\psi(x) = -0.02$ \\, ($\\beta=%.3g,\\kappa=%.3g$)', betaVec(k),kappa), 'Interpreter','latex');
%     xlabel('$Time$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); zlabel('$u(x,t$)','Interpreter','latex'); prettify(gca);
% end

% fname4 = sprintf('displacement_surfaces_k%s_b%s.pdf', kappa_tag, beta_tag);
% exportgraphics(fig4, fname4, 'Resolution', 600, 'BackgroundColor','white');

%% ====== Overlaid contact metrics ======
% fig5 = figure('Color','w','Name','Contact metrics');
% tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
% 
% nexttile; hold on;
% for k = 1:nb, plot(t, sol{k}.frac_contact, 'LineWidth',1.8, 'Color', cols(k,:)); end
% ylabel('Fraction in Contact','Interpreter','latex');
% title('Contact Fraction vs Time','Interpreter','latex');grid off;
% legend(arrayfun(@(b)sprintf('$\\beta=%.3g$',b), betaVec, 'uni',0), ...
%        'Interpreter','latex', 'Location','best');
% prettify(gca);
% 
% nexttile; hold on;
% for k = 1:nb, plot(t, sol{k}.max_pen, 'LineWidth',1.8, 'Color', cols(k,:)); end
% ylabel('Max penetration','Interpreter','latex'); xlabel('Time','Interpreter','latex');
% title('Maximum Penetration vs Time','Interpreter','latex');
% prettify(gca);

end  % main

% ------------ Helper: solve for a single beta with ψ(x) ------------
function S = solve_single_beta_psi(beta, L,T,dx,dt, alpha2,kappa,psiFun,f0, phi,u0f,v0f)
M = round(L/dx);  N = round(T/dt);
x = linspace(0,L,M+1).';  t = linspace(0,T,N+1);
psi = psiFun(x);                        % obstacle profile at nodes
u = zeros(M+1, N+1);
gamma = alpha2*(dt^2)/(dx^4);

% Initial data (2nd-order start)
u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);  % Dirichlet

for j = 1:(N-1)
    I = (u(:, j+1) <= psi);  I(1) = false;   % contact test against ψ(x)
    diag_add = dt^2 * (kappa * double(I(2:end)) + (beta/dt) * double(I(2:end)));
    rhs_add  = dt^2 * (kappa * psi(2:end) .* double(I(2:end)) ...
                      + (beta/dt) * double(I(2:end)) .* u(2:end, j+1));

    A = spalloc(M, M, 5*M);  b = zeros(M,1);

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
    u(1,   j+2) = phi(t(j+2));
    u(2:end,j+2) = w;
end

% Derived fields (vs ψ)
contact     = (u <= psi);
penetration = max(psi - u, 0);
frac_contact = mean(contact,1);
max_pen     = max(penetration,[],1);

% Pack
S.x=x; S.t=t; S.u=u; S.psi=psi;
S.contact=contact; S.penetration=penetration;
S.frac_contact=frac_contact; S.max_pen=max_pen;
S.beta = beta;
end

% ------------ Utilities ------------
function prettify(ax)
set(ax, 'FontSize', 12, 'LineWidth', 1.0);
box(ax, 'on'); grid(ax, 'on');
end

function [nr,nc] = tilesize(n)
nr = ceil(sqrt(n)); nc = ceil(n/nr);
end


function animate_beam_distributed_2d(S, kappa, make_video, video_name)
% Minimal 2D animation of u(x,t) over a distributed obstacle ψ(x)
% Colors: blue (free), red (in contact). Matches beam_contact_dnc_anim style.

x   = S.x;          % column (M+1)x1
t   = S.t;          % row 1x(N+1)
u   = S.u;          % (M+1)x(N+1)
psi = S.psi;        % (M+1)x1
L   = x(end);
dt  = t(2)-t(1);

% Figure
fig = figure('Color','w','Name','Beam–Obstacle Animation (Distributed)');
hold on; grid on; box on;

% Obstacle curve ψ(x)
plot(x, psi, 'k--', 'LineWidth', 1.5, 'DisplayName','\psi(x)');

% Beam lines: free (blue) and contact (red) parts
hFree = plot(x, nan(size(x)), 'b-', 'LineWidth', 2.5, 'DisplayName','beam (free)');
hCont = plot(x, nan(size(x)), 'r-', 'LineWidth', 2.5, 'DisplayName','beam (contact)');

% Clamp indicator at x=0 (u(0,t)=0 here)
plot(0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName','clamp');

% Tip marker
hTip = plot(x(end), u(end,1), 'o', 'MarkerSize', 6, ...
            'MarkerFaceColor','b', 'MarkerEdgeColor','b', 'DisplayName','tip');

xlabel('x'); ylabel('u(x,t)');
ymin = min([u(:); psi(:)]) - 0.05;
ymax = max([u(:); psi(:)]) + 0.05;
axis([0 L ymin ymax]);

ttl = title(sprintf('Distributed DNC  (\\beta=%.3g, \\kappa=%.3g, t=%.2f)', S.beta, kappa, t(1)));
set(ttl, 'Interpreter', 'tex');
legend('Location','best');

% Optional MP4
if make_video
    vw = VideoWriter(video_name, 'MPEG-4'); vw.FrameRate = 25; open(vw);
end

% Animate (blue where u>ψ, red where u≤ψ)
skip = 1;              % show every frame
playback = 1.0;        % sim-sec per real-sec (1.0 ≈ real-time)

for j = 1:skip:numel(t)
    uj   = u(:,j);
    mask = (uj <= psi);           % contact at nodes

    yCont = uj;  yCont(~mask) = NaN;
    yFree = uj;  yFree(mask)  = NaN;

    set(hCont, 'YData', yCont);
    set(hFree, 'YData', yFree);

    set(hTip, 'XData', x(end), 'YData', uj(end));
    if mask(end)
        set(hTip, 'MarkerFaceColor','r','MarkerEdgeColor','r');
    else
        set(hTip, 'MarkerFaceColor','b','MarkerEdgeColor','b');
    end

    set(ttl, 'String', sprintf('Distributed DNC  (\\beta=%.3g, \\kappa=%.3g, t=%.2f)', S.beta, kappa, t(j)));
    drawnow;

    if make_video, writeVideo(vw, getframe(fig)); end
    pause((skip*dt)/max(playback, eps));
end

if make_video, close(vw); end
end
