% montage_uLt.m — 2x2 time-trace montage for (alpha^2, kappa) at fixed beta & force
clc; close all;

% Fixed parameters
beta   = 0.1;         % damping
Fconst = -1.2;        % force
yminus = -0.02;       % obstacle level

% Soft/stiff beam and soft/hard obstacle
alpha2_soft = 2;   alpha2_stiff = 10;
kappa_soft  = 0.01; kappa_hard  = 1000;

% Solver defaults
base = struct('T',5,'dt',0.01,'dx',0.01,'y_minus',yminus);

% Define the four cases (row = beam stiffness, col = obstacle stiffness)
cases = { ...
  struct('a2',alpha2_soft, 'kap',kappa_soft, 'name','Soft beam v Soft obstacle'); ...
  struct('a2',alpha2_soft, 'kap',kappa_hard, 'name','Soft beam v Hard obstacle'); ...
  struct('a2',alpha2_stiff,'kap',kappa_soft, 'name','Stiff beam v Soft obstacle'); ...
  struct('a2',alpha2_stiff,'kap',kappa_hard, 'name','Stiff beam v Hard obstacle') ...
};

tips  = cell(4,1); metas = cell(4,1); tgrid = []; 
for k = 1:4
    p = base; p.alpha2 = cases{k}.a2;
    [~,t,u,m] = beam_contact_implicit_FD(Fconst, cases{k}.kap, beta, p);
    tips{k}   = m.tip; metas{k} = m; tgrid = t;
end

% Shared axes limits for fair comparison
allVals = cell2mat(tips');
ymin = min(allVals); ymax = max(allVals);
pad  = 0.05*(ymax - ymin + eps);
yl   = [ymin - pad, ymax + pad];

% Choose a color set for the four plots
colors = lines(4);   % or use hsv(4), parula(4), etc.

% Plot 2x2 montage
f = figure('Color','w','Name','uLt_montage_alpha2_kappa');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for k = 1:4
    ax = nexttile; hold(ax,'on'); grid(ax,'off');
    plot(ax, tgrid, tips{k}, 'LineWidth', 1.6, 'Color', colors(k,:));
    yline(ax, yminus, 'k--', 'LineWidth', 1);
    ylim(ax, yl); xlabel(ax,'$Time$','Interpreter','latex'); ylabel(ax,'$u(1,t)$','Interpreter','latex');
    title(ax, sprintf('\\ %s ($\\alpha^2=%g$, $\\kappa=%g$)', cases{k}.name, cases{k}.a2, cases{k}.kap), ...
         'Interpreter','latex');
end

sgtitle(sprintf('Beam Tip Displacement: Soft vs Stiff Beam and Obstacle: $\\beta=%g$,  $f$=%g,  $y_-=%g$', beta, Fconst, yminus), ...
        'Interpreter','latex');

exportgraphics(f,'uLt_montage_alpha2_kappa_beta0p1_F-0p25.pdf','ContentType','vector');
% saveas(f,'uLt_montage_alpha2_kappa_beta0p1_F-0p25.jpg');
