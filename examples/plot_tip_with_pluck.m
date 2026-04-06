% run_alpha2_kappa.m  — sweep alpha^2 for multiple kappa, fixed beta & force, with time excitation phi(t)
clc; close all;

% --- sets to study ---
alpha2Set = [2,4,6,8,10];   % beam stiffness proxy
kappaSet  = 0.1;           % obstacle stiffness (can be scalar or vector)
betaSet   = 5;            % damping (can be scalar or vector)
Fconst    = -0.2;           % constant external forcing (sign/magnitude as needed)

% --- boundary time excitation: u(0,t) = phi(t) = A*sin(2*pi*f*t) ---
Aexc = -0.005;               % displacement amplitude at x=0
fHz  = 1.5;                 % excitation frequency (Hz)
phi_t = @(tt) Aexc * sin(2*pi*fHz*tt);

% optional overrides into your solver (grid, horizon, obstacle, boundary)
params = struct('T', 10, 'dt', 0.01, 'dx', 0.01, 'y_minus', -0.02, 'phi_t', phi_t);

% helpers
tag = @(x) strrep(sprintf('%.3g',x),'.','p');   % filename-safe numeric
fmt = @(x) sprintf('%.6g',x);                   % clean numeric in titles/labels

for K = kappaSet
  for B = betaSet

    % --- storage over alpha2 ---
    nA   = numel(alpha2Set);
    tips = cell(nA,1);
    metas= cell(nA,1);
    tgrid= [];
    yobs = NaN;

    % consistent colors per alpha2 across both figures
    colsA = lines(nA);

    % --- run all alpha2 for this (K,B) ---
    for ia = 1:nA
        A2 = alpha2Set(ia);
        p  = params; 
        p.alpha2 = A2;                 % set beam stiffness
        [x,t,u,met] = beam_contact_implicit_FD(Fconst, K, B, p);
        tips{ia}  = met.tip;
        metas{ia} = met;
        tgrid     = t;                 % same grid each run
        yobs      = met.y_minus;
    end

    %====================== Figure 1: tip overlays (colors by alpha^2) ======================
    f1Name = sprintf('Tip_vs_Time_k%s_b%s_f%s_A%s_fHz%s', ...
                     tag(K), tag(B), tag(Fconst), tag(Aexc), tag(fHz));
    f1 = figure('Color','w','Name', f1Name); hold on; grid off;

    lgdLbl = cell(nA,1);
    for ia = 1:nA
        plot(tgrid, tips{ia}, 'LineWidth', 1.4, 'Color', colsA(ia,:));
        lgdLbl{ia} = ['\alpha^2 = ', fmt(alpha2Set(ia))];
    end
    hobs = yline(yobs,'k--','LineWidth',1.2,'LabelHorizontalAlignment','left','Interpreter','latex');
    hobs.DisplayName = sprintf('y_- = %.3g', yobs);

    xlabel('$Time$', 'Interpreter','latex'); ylabel('$u(1,t)$', 'Interpreter','latex');
    line1 = sprintf('Right Tip Displacement $\\alpha^2$ sweep ($\\kappa=%s$, $\\beta=%s$, $f=%.3g$) with', ...
                fmt(K), fmt(B), Fconst);
    line2 = sprintf('boundary excitation: $\\varphi(t)=A\\,\\sin(2\\pi f_q t)$ with $A=%.3g$, $f_q=%.3g\\,\\mathrm{Hz}$', ...
                Aexc, fHz);

    title({line1; line2}, 'Interpreter','latex');



    legend([lgdLbl; {sprintf('y_- = %.3g', yobs)}], 'Location','northeast');
    hold off;
    % saveas(f1, [f1Name '.jpg']);
    % exportgraphics(f1, [f1Name '.pdf'], 'ContentType','vector');

    %====================== Figure 2: summaries vs alpha^2 ======================
    % collect vectors
    t_first_vec      = cellfun(@(m)m.t_first,      metas);
    num_impacts_vec  = cellfun(@(m)m.num_impacts,  metas);
    dmax_vec         = cellfun(@(m)m.dmax,         metas);
    contact_frac_vec = cellfun(@(m)m.contact_frac, metas);

    f2Name = sprintf('Metrics_vs_alpha2_k%s_b%s_f%s_A%s_fHz%s', ...
                     tag(K), tag(B), tag(Fconst), tag(Aexc), tag(fHz));
    f2 = figure('Color','w','Name', f2Name);
    tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    % tile 1: t_c
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'off');
    plot(ax1, alpha2Set, t_first_vec, '-o','LineWidth',1.2,'Color','m');
    xlabel('\alpha^2'); ylabel('t_c'); title('First-contact time, t_c');

    % tile 2: impact count
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'off');
    plot(ax2, alpha2Set, num_impacts_vec, '-o','LineWidth',1.2,'Color','g');
    xlabel('\alpha^2'); ylabel('count'); title('Impact count');

    % tile 3: d_{max}
    ax3 = nexttile; hold(ax3,'on'); grid(ax3,'off');
    plot(ax3, alpha2Set, dmax_vec, '-o','LineWidth',1.2,'Color','r');
    xlabel('\alpha^2'); ylabel('d_{max}'); title('Max penetration, d_{max}');

    % tile 4: contact fraction
    ax4 = nexttile; hold(ax4,'on'); grid(ax4,'off');
    plot(ax4, alpha2Set, contact_frac_vec, '-o','LineWidth',1.2,'Color','b');
    xlabel('\alpha^2'); ylabel('T_c/T'); title('Contact fraction, T_c/T');

    line1 = sprintf('Contact Metrics vs. Beam Stiffness $\\alpha^2$ ($\\kappa=%s$, $\\beta=%s$, $f=%.3g$) with', ...
                fmt(K), fmt(B), Fconst);
    line2 = sprintf('boundary excitation: $\\varphi(t)=A\\,\\sin(2\\pi f_q t)$ with $A=%.3g$, $f_q=%.3g\\,\\mathrm{Hz}$', ...
                Aexc, fHz);

    sgtitle({line1; line2}, 'Interpreter','latex');

    % exportgraphics(f2, [f2Name '.pdf'], 'ContentType','vector');

  end
end
