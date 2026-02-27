% stiffness_maps.m
clc; close all;

alpha2Set = logspace(1,10);     % beam stiffness proxy
kappaSet  = logspace(1,1000);     % obstacle stiffness
beta      = 0.01;                    % pick light (e.g., 0.001) or strong (e.g., 1)
force    = -1.5;                    % one representative force (or iterate a few)

% matrices: rows=alpha2, cols=kappa
tc   = NaN(numel(alpha2Set), numel(kappaSet));
dmax = tc; cnt = tc; frac = tc;

base = struct(); base.T=5; base.dt=0.01; base.dx=0.01; base.y_minus=-0.02;

for ia = 1:numel(alpha2Set)
  for ik = 1:numel(kappaSet)
    params = base; params.alpha2 = alpha2Set(ia);
    [~,t,~,m] = beam_contact_implicit_FD(force, kappaSet(ik), beta, params);
    tc(ia,ik)   = m.t_first;
    dmax(ia,ik) = m.dmax;
    cnt(ia,ik)  = m.num_impacts;
    frac(ia,ik) = m.contact_frac;
  end
end

% --- helper to draw a log-log heatmap-like image ---
f = figure('Color','w','Name','Stiffness maps');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
A = log10(alpha2Set); K = log10(kappaSet);

nexttile; imagesc(K, A, tc);   axis xy; colorbar; title('t_c');
set(gca,'XTick',K, 'XTickLabel',compose('10^{%g}',K), ...
        'YTick',A, 'YTickLabel',compose('10^{%g}',A));
xlabel('\kappa'); ylabel('\alpha^2');

nexttile; imagesc(K, A, dmax); axis xy; colorbar; title('d_{max}');
set(gca,'XTick',K,'XTickLabel',compose('10^{%g}',K),'YTick',A,'YTickLabel',compose('10^{%g}',A));
xlabel('\kappa'); ylabel('\alpha^2');

nexttile; imagesc(K, A, cnt);  axis xy; colorbar; title('Impact count');
set(gca,'XTick',K,'XTickLabel',compose('10^{%g}',K),'YTick',A,'YTickLabel',compose('10^{%g}',A));
xlabel('\kappa'); ylabel('\alpha^2');

nexttile; imagesc(K, A, frac); axis xy; colorbar; title('Contact fraction T_c/T');
set(gca,'XTick',K,'XTickLabel',compose('10^{%g}',K),'YTick',A,'YTickLabel',compose('10^{%g}',A));
xlabel('\kappa'); ylabel('\alpha^2');

sgtitle(sprintf('Maps at \\beta=%.3g, F=%.3g', beta, force), 'Interpreter','tex');
