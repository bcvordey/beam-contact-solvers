% threshold_surface.m
clc; close all;

% --- grids to scan ---
alpha2Set = logspace(-1, 2, 7);   % beam stiffness proxy
kappaSet  = logspace(-2, 6, 8);   % obstacle stiffness
beta      = 0.1;                  % choose a damping level to study
Fgrid     = -logspace(-3, 1, 20); % try forces from small to large magnitude (negative sign if needed)

% --- solver defaults (override here if you want) ---
base = struct('T',5,'dt',0.01,'dx',0.01,'y_minus',-0.02);

% Fstar will hold the minimal |F| that yields contact; NaN if none within scan
Fstar = NaN(numel(alpha2Set), numel(kappaSet));

for ia = 1:numel(alpha2Set)
  for ik = 1:numel(kappaSet)
    params = base; params.alpha2 = alpha2Set(ia);

    firstHit = NaN;  % default: no contact within scanned forces
    for m = 1:numel(Fgrid)
      [~,~,~,met] = beam_contact_implicit_FD(Fgrid(m), kappaSet(ik), beta, params);
      if any(met.contact)                % contact occurred at this force
        firstHit = abs(Fgrid(m));        % minimal |F| found (since Fgrid is ordered small->large)
        break;
      end
    end

    Fstar(ia,ik) = firstHit;
  end
end

% --- plot: contour over log10 axes (rows=alpha2, cols=kappa) ---
[KK, AA] = meshgrid(log10(kappaSet), log10(alpha2Set));  % sizes match Fstar (alpha x kappa)

figure('Color','w','Name','Threshold force for first contact');
contourf(KK, AA, Fstar, 12, 'LineColor','none'); 
cb = colorbar; ylabel(cb,'|F| threshold');
xlabel('log_{10} \kappa'); ylabel('log_{10} \alpha^2');
title(sprintf('Threshold |F| for first contact (\\beta=%.3g)', beta), 'Interpreter','tex');
grid on;
