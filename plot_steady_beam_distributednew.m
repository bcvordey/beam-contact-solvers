function plot_steady_beam_distributednew()
% Steady-state beam over a distributed obstacle: piecewise (no-contact/contact)
% Uses regime = check_contact_consistency_alpha2(alpha2,f,y_minus)
% and piecewise_u(...) to pick the correct closed form.

% ----- Parameters -----
alpha2  = 1;          
kappa   = 80;          
y_minus = -0.2;       
fvals   = [1.0, 0.0, -1.0, -1.6, -2.0, -2.5, -3.0, -3.5];

% ----- Grid -----
x = linspace(0,1,1000);

% ----- Figure -----
figure('Color','w'); hold on; grid on;

for k = 1:numel(fvals)
    f = fvals(k);

    % piecewise steady solution based on regime check
    u = piecewise_u(x, alpha2, kappa, f, y_minus);

    % regime label for legend
    regime = check_contact_consistency_alpha2(alpha2, f, y_minus);
    plot(x, u, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('f = %.3g (%s)', f, strrep(regime,'_','-')));
end

% Obstacle line
yline(y_minus, '--k', 'LineWidth', 1.6, 'DisplayName','obstacle');

xlabel('$x$', 'Interpreter','latex'); 
ylabel('$\bar{u}(x)$', 'Interpreter','latex');
title(sprintf('Steady-State Beam over Distributed Obstacle ($\\kappa = %.3g$)', kappa), ...
      'Interpreter','latex');
legend('Location','best','Interpreter','latex'); 
box on; grid on;

end

% ==================== Helpers ====================

function regime = check_contact_consistency_alpha2(alpha2, f, y_minus)
% Decide regime using the α²-based threshold: threshold = f/(8 α²)
% If y_- is above the threshold (>=), the beam contacts; otherwise no-contact.
    threshold = f / (8*alpha2);
    if y_minus < threshold
        regime = 'no_contact';
    else
        regime = 'contact';
    end
end

function val = piecewise_u(x, alpha2, kappa, f, y_minus)
% Wrapper that switches between the two closed forms
    regime = check_contact_consistency_alpha2(alpha2, f, y_minus);
    switch regime
        case 'no_contact'
            val = u_no_contact(x, alpha2, f);
        otherwise % 'contact'
            val = u_contact(x, alpha2, kappa, f, y_minus);
    end
end

function u = u_no_contact(x, alpha2, f)
% No-contact steady solution: α² u'''' = f with clamp-free BC
% \bar{u}(x) = (f/(24 α²)) (x^2 - 4x + 6) x^2
    u = (f/(24*alpha2)) .* ((x.^2 - 4*x + 6) .* x.^2);
end

function u = u_contact(x, alpha2, kappa, f, y_minus)
% Contact steady solution via 4x4 system, solved with backslash
% α² u'''' - κ u = f - κ y_-, clamp-free BCs
    if alpha2 <= 0 || kappa <= 0
        error('alpha2 and kappa must be positive.');
    end

    lambda = (kappa/alpha2)^(1/4);
    d = y_minus - f/kappa;

    % shorthands at x=1
    ch = cosh(lambda); sh = sinh(lambda);
    c  = cos(lambda);  s  = sin(lambda);

    % System for [A;B;C;D]:
    % 1) A + C = -d
    % 2) B + D = 0
    % 3) ch*A + sh*B - c*C - s*D = 0
    % 4) sh*A + ch*B + s*C - c*D = 0
    M = [ 1   0    1    0;
          0   1    0    1;
          ch  sh  -c   -s;
          sh  ch   s   -c ];
    b = [-d; 0; 0; 0];

    X = M \ b;  A = X(1); B = X(2); C = X(3); D = X(4);

    % Evaluate u(x)
    u = d ...
      + A*cosh(lambda.*x) ...
      + B*sinh(lambda.*x) ...
      + C*cos (lambda.*x) ...
      + D*sin (lambda.*x);
end
