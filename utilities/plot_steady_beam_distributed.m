function plot_steady_beam_distributed()
% Parameters
alpha2  = 1;          
kappa   = 1;          
y_minus = -0.2;       
fvals   = [1.0, 0.0, -1.0, -1.6, -2.0, -2.5, -3.0, -3.5];

% Grid
x = linspace(0,1,1000);

% Figure
figure(); hold on; grid on;

for k = 1:numel(fvals)
    f = fvals(k);
    threshold = f/(8*alpha2);

    if (y_minus <= threshold)
        u = (f/(24*alpha2)) .* ( (x.^2 - 4*x + 6) .* x.^2 );
    else
        u = contact_solution_closed_form(x, alpha2, kappa, f, y_minus);
        
    end
    
    plot(x, u, 'LineWidth', 1.6, 'DisplayName', sprintf('f = %.3g', f));
end

% Obstacle line
yline(y_minus, '--k', 'LineWidth', 1.6,'DisplayName','obstacle');

xlabel('$x$', 'Interpreter','latex'); 
ylabel('$\bar{u}(x)$', 'Interpreter','latex');
title(sprintf('Steady-State Beam over Distributed Obstacle $(\\kappa = %.3g)$', kappa), 'Interpreter','latex');
legend('Location', 'best', 'Interpreter', 'latex'); 
box on; grid off;

end

%===================== Helpers =====================%
% function u = contact_solution_closed_form(x, alpha2, kappa, f, y_minus)
% lambda = (kappa/alpha2)^(1/4);
% d = y_minus - f/kappa;
% 
% ch = cosh(lambda);
% sh = sinh(lambda);
% c  = cos(lambda);
% s  = sin(lambda);
% 
% den = -2 * (1 + ch*c);
% A   = d * (1 + ch*c + sh*s) / den;
% D   = d * (sh*c + ch*s) / den;
% B   = -D;
% C   = -d - A;
% 
% u = d + A*cosh(lambda*x) + B*sinh(lambda*x) + C*cos(lambda*x) + D*sin(lambda*x);
% 
function u = contact_solution_closed_form(x, alpha2, kappa, f, y_minus)

    if alpha2 <= 0 || kappa <= 0
        error('alpha2 and kappa must be positive.');
    end

    % parameters
    lambda = (kappa/alpha2)^(1/4);
    d = y_minus - f/kappa;

    % shorthands at x=1 for the boundary system
    ch = cosh(lambda);
    sh = sinh(lambda);
    c  = cos(lambda);
    s  = sin(lambda);

    M = [ 1   0    1    0;
          0   1    0    1;
          ch  sh  -c   -s;
          sh  ch   s   -c ];

    b = [-d; 0; 0; 0];

    X = M \ b;
    A = X(1); B = X(2); C = X(3); D = X(4);
    u = d ...
        + A*cosh(lambda.*x) ...
        + B*sinh(lambda.*x) ...
        + C*cos (lambda.*x) ...
        + D*sin (lambda.*x);
end


