function steadystate_1obs()
    % Parameters
    alpha   = 1;
    kappa   = 1000;
    y_minus = -0.2;   % obstacle 0.2 below x-axis
    
    % Try various f
    % (Though the derivation assumes f=const, we'll just see how the code behaves.)
    f_values = [1,0,-1, -1.6,-2, -2.5, -3.0, -3.5, -25];
    
    xgrid = linspace(0, 1, 200);
    
    figure; hold on; grid off;
    xlabel('$x$',Interpreter='latex'); ylabel('$u(x)$',Interpreter='latex');
    title(sprintf('$\\kappa = %.0f$', kappa), ...
      'Interpreter','latex');

    
   
    
    for f = f_values
        % Evaluate piecewise solution
        uvals = arrayfun(@(x) piecewise_u(x, alpha, kappa, f, y_minus), xgrid);
        plot(xgrid, uvals, 'DisplayName', sprintf('f = %.1f', f),'LineWidth', 1.5);
    end

    yl = ylim;                            % current limits
    
    % ensure there is *some* room below y_minus and above y_plus
    tol = max(1e-10, 0.01*diff(yl));       % small, automatic epsilon
    ymin = min(yl(1), y_minus - tol);
    % ymax = max(yl(2), y_plus  + tol);
    % ylim([ymin, ymax]);  

     % Draw y_minus line
    % ymax = 1.05*max([1e-12, max(abs(uvals),[],'all'), abs(y_minus)]);
    % yline(y_minus, 'k--','Label',sprintf('y_{-} = %.2f', y_minus),'DisplayName', 'obstacle','LineWidth', 1.3);
    % line([1 1], [-1.1*ymax, y_minus],'Color','k', 'LineWidth',5, 'LineStyle','--','DisplayName', 'obstacle');   % dashed % Vertical bar

    line([1 1], [y_minus, ymin], 'Color','k', 'LineWidth',5,'LineStyle','-', 'DisplayName', 'obstacle');

    legend('Location','southwest');
    hold off;box on;
end

function regime = check_contact_consistency(alpha, ~, f, y_minus)
    % The critical "end displacement" is f/(8*alpha^2).
    % If y_minus is still above that, the beam does NOT contact.
    threshold = f/(8*alpha^2);
    if y_minus < threshold
        regime = 'no_contact';
    else
        regime = 'contact';
    end
end

function val = piecewise_u(x, alpha, kappa, f, y_minus)
    regime = check_contact_consistency(alpha, kappa, f, y_minus);
    if strcmp(regime, 'no_contact')
        val = u_no_contact(x, alpha, f);
    else
        val = u_contact(x, alpha, kappa, f, y_minus);
    end
end

%--- No-contact formula ----------------------
function val = u_no_contact(x, alpha, f)
   
    val = (f./(24*alpha.^2)) .* ((x.^2 - 4*x + 6) .* x.^2);
end





%--- Contact formula -------------------------
function val = u_contact(x, alpha, kappa, f, y_minus)

    A = a_contact(alpha, kappa, f, y_minus);
    val = (f/(24*alpha^2))*x.^4 + (A/6)*x.^3 - ((f/(4*alpha^2))+A/2)*x.^2;
end

function A = a_contact(alpha, kappa, f, y_minus)
   
   A = (-3*kappa / (kappa + 3*alpha^2)) * ...
    ( y_minus + f * (5/(24*alpha^2) + 1/kappa) );
end

