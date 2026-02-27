% run_force_influence.m
clc; close all;

kappaSet = [0.01,100,1000];
betaSet  = [0.001, 1];
forceVec = [-0.05, -0.10, -0.15, -0.20, -0.25, -0.30];

params = struct('T',10);     % optional overrides
tag = @(x) strrep(sprintf('%.3g',x),'.','p');  % filename-safe

% --- fixed color map: one color per force, reused across figures ---
nF   = numel(forceVec);
cols = lines(nF);  % or parula(nF), turbo(nF), etc.

for K = kappaSet
  for B = betaSet
    tips   = cell(nF,1);
    metas  = cell(nF,1);
    tgrid  = [];
    yobs   = NaN;

    % --- run all forces for this (K,B) ---
    for m = 1:nF
        f0 = forceVec(m);
        [x,t,u,met] = beam_contact_implicit_FD(f0,K,B,params);
        tips{m}  = met.tip;
        metas{m} = met;
        tgrid    = t;
        yobs     = met.y_minus;
    end

    % -------- Figure 1: tip overlays (match colors by force) --------
    f1 = figure('Color','w','Name',sprintf('Tip_vs_Time_k%s_b%s',tag(K),tag(B)));
    hold on;
    for m = 1:nF
        plot(tgrid, tips{m}, 'LineWidth', 1.2, ...
            'Color', cols(m,:), ...
            'DisplayName', sprintf('f=%.3g', forceVec(m)));
    end
    hobs = yline(yobs,'k--', 'LineWidth',1.2, 'LabelHorizontalAlignment','left', 'Interpreter','latex');
    hobs.DisplayName = sprintf('y_- = %.3g', yobs);
    xlabel('$Time$','Interpreter','latex'); ylabel('$u(1,t)$','Interpreter','latex');grid off;
    title(sprintf('Tip Displacement $u(1,t)$ vs. Time - Force Sweep ($\\kappa=%.6g$, $\\beta=%.6g$)', K, B), ...
      'Interpreter','latex');
    legend('Location','southeast'); hold off;
    % saveas(f1, sprintf('Tip_vs_Time_k%s_b%s.pdf', tag(K), tag(B)));
    exportgraphics(f1, sprintf('Tip_vs_Time_k%s_b%s.pdf', tag(K), tag(B)), 'ContentType','vector'); 

    % -------- collect summaries --------
    t_first_vec      = cellfun(@(s)s.t_first, metas);
    num_impacts_vec  = cellfun(@(s)s.num_impacts, metas);
    dmax_vec         = cellfun(@(s)s.dmax, metas);
    contact_frac_vec = cellfun(@(s)s.contact_frac, metas);

    % -------- Figure 2: summaries vs force (colored markers by force) --------
    f2 = figure('Color','w','Name',sprintf('Force_Summary_k%s_b%s',tag(K),tag(B)));
    tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    % helper for per-force colored scatter
    mk = 60;  % marker size
    lab = arrayfun(@(v) sprintf('f=%.3g',v), forceVec, 'UniformOutput', false);

    % Tile 1: t_c
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'off');
    plot(ax1, forceVec, t_first_vec, '-*','LineWidth',1.2,'Color','m');
    xlabel('$f$','Interpreter','latex'); ylabel('$t_c$','Interpreter','latex'); title('Time to first contact','Interpreter','latex');

    % Tile 2: contact events
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'off');
    plot(ax2, forceVec, num_impacts_vec, '-*','LineWidth',1.2,'Color','g');
    xlabel('$f$','Interpreter','latex'); ylabel('count'); title('Impact count','Interpreter','latex');

    % Tile 3: d_max
    ax3 = nexttile; hold(ax3,'on'); grid(ax3,'off');
    plot(ax3, forceVec, dmax_vec, '-*','LineWidth',1.2,'Color','r');
    xlabel('$f$','Interpreter','latex'); ylabel('$d_{max}$','Interpreter','latex'); title('Max penetration','Interpreter','latex');

    % Tile 4: contact fraction
    ax4 = nexttile; hold(ax4,'on'); grid(ax4,'off');
    plot(ax4, forceVec, contact_frac_vec, '-*','LineWidth',1.2,'Color','b');
    xlabel('$f$','Interpreter','latex'); ylabel('T_c/T'); title('Contact fraction','Interpreter','latex');

    % One legend, placed for the whole layout:
    % lgd = legend(ax1, ph, lab, 'Location','northwest');
    % lgd.Title.String = 'Force color map';
    sgtitle(sprintf('Contact Metrics vs. Force: $\\kappa=%.6g$, $\\beta=%.6g$', K, B), 'Interpreter','latex');

    % saveas(f2, sprintf('Force_Summary_k%s_b%s.pdf', tag(K), tag(B)));
    exportgraphics(f2, sprintf('Force_Summary_k%s_b%s.pdf', tag(K), tag(B)), 'ContentType','vector');

  end
end
