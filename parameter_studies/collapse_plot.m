% collapse_plot.m
clc; close all;
beta = 0.1; force = -1.5;
alpha2Set = [0.1, 1, 10];     % four beams
kappaSet  = logspace(1, 10, 100);
dx=0.01; base = struct('T',5,'dt',0.01,'dx',dx,'y_minus',-0.02);

figure('Color','w'); hold on; grid on;
cols = lines(numel(alpha2Set));
for ia = 1:numel(alpha2Set)
  params = base; params.alpha2 = alpha2Set(ia);
  chi = nan(size(kappaSet)); Dm = chi;
  for ik = 1:numel(kappaSet)
    [~,~,~,m] = beam_contact_implicit_FD(force, kappaSet(ik), beta, params);
    chi(ik) = (dx^3/params.alpha2)*kappaSet(ik);
    Dm(ik)  = m.dmax;
  end
  plot(chi, Dm, '-o', 'Color', cols(ia,:), 'DisplayName', sprintf('\\alpha^2=%.3g', alpha2Set(ia)));
end
set(gca,'XScale','log'); xlabel('\chi = \kappa dx^3/\alpha^2'); ylabel('d_{max}');
title('Collapse of penetration vs. non-dimensional stiffness \chi'); legend('Location','best');
