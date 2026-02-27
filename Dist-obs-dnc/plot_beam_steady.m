function plot_beam_steady()
% plot_beam_steady
% Steady solution of:  alpha^2 * u''''(x) - kappa * u(x) = f - kappa*y_-
% with BCs: u(0)=0, u'(0)=0, u''(1)=0, u'''(1)=0.
%
% This version accepts vectors for kappa and f, makes one figure per kappa,
% and plots a family of curves for the different f values.

    % ---------------- User parameters ----------------
    alpha2  = 1;                                 % α^2 (> 0)
    kappa   = [0.1,1];                     % κ values (vector, > 0)
    y_minus = -0.2;                              % obstacle level (constant)
    f       = [1.0, 0.0, -1.0, -1.6, -2.0, ...   % force values (vector)
               -2.5, -3.0, -3.5, -4.0, -5.0];
    % -------------------------------------------------

    % Discretization in x
    nx = 1200;
    x  = linspace(0,1,nx);

    % Loop over kappa
    for ik = 1:numel(kappa)
        kap = kappa(ik);

        figure('Color','w');
        hold on; grid on;

        % Loop over f values for this kappa
        for jf = 1:numel(f)
            ff = f(jf);

            threshold = ff/(8*alpha2);

            if (y_minus <= threshold)
                u = (ff/(24*alpha2)) .* ( (x.^2 - 4*x + 6) .* x.^2 );
                % tag = 'no-contact';
            else
                u = steady_solution_profile(alpha2, kap, ff, y_minus, x);
                % tag = 'contact';
            end

            % Compute coefficients and solution
            

            plot(x, u, 'LineWidth', 1.6, 'DisplayName', sprintf('f = %.3g', ff));

            % Quick BC residual printout
            % fprintf("kappa=%g, f=%g | BC residuals: u(0)=%.2e, u''(0)=%.2e, u''''(1)=%.2e, u'''''(1)=%.2e,  Delta=%.3e\n", kap, ff, meta.u0, meta.up0, meta.uxx1, meta.uxxx1, meta.Delta);
        end

        xlabel('x','Interpreter','tex');
        ylabel('$\bar{u}(x)$','Interpreter','latex');
        title(sprintf('Steady solution family for \\kappa=%.3g (\\alpha^2=%.3g, y_-=%.3g)', ...
              kap, alpha2, y_minus), 'Interpreter','tex');
        legend('Location','bestoutside'); box on; hold off;
    end
end

% ======================================================================
% Helper: compute steady solution for given (alpha2, kappa, f, y_minus)
% ======================================================================
function [u, meta] = steady_solution_profile(alpha2, kappa, f, y_minus, x)
    % Parameters
    lam  = (kappa/alpha2)^(1/4);     % lambda = (kappa/alpha^2)^(1/4)
    up   = y_minus - f/kappa;        % constant particular solution

    % Shorthands at x=1 (for coefficients)
    ch = cosh(lam);  sh = sinh(lam);
    co = cos(lam);   si = sin(lam);

    a = ch + co;
    b = sh + si;
    c = sh - si;
    Delta = 2*(1 + ch*co);           % system determinant

    if abs(Delta) < 1e-10
        warning('Near-singular system for lambda=%.6g: Delta ~ 0. Results may be ill-conditioned.', lam);
    end

    % Solve for A, B
    A = up * ( -a*co - b*si ) / Delta;
    B = up * (  a*si + c*co ) / Delta;

    % Recover C, D from left-end BCs
    C = -(up + A);
    D = -B;

    % Evaluate u(x)
    lamx = lam * x;
    u = up ...
      + A*cosh(lamx) ...
      + B*sinh(lamx) ...
      + C*cos(lamx)  ...
      + D*sin(lamx);

    % (Optional) BC residuals for quick verification
    u0    = up + A + C;
    up0   = lam*(B + D);
    uxx1  = lam^2*( A*ch + B*sh - C*co - D*si );
    uxxx1 = lam^3*( A*sh + B*ch + C*si - D*co );

    meta = struct('lambda',lam,'u_p',up,'Delta',Delta, ...
                  'u0',u0,'up0',up0,'uxx1',uxx1,'uxxx1',uxxx1, ...
                  'A',A,'B',B,'C',C,'D',D);
end
