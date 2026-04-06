function steadystatetwoobstacles()
    alpha   = 1;
    kappa   = 1000;
    y_minus = -0.1;
    y_plus  =  0.1;
    assert(y_minus < y_plus, 'Require y_- < y_+');

    % f_values = [3, 2, 1.6, 1, 0, -1, -1.6, -2, -3];
    f_values = [25, 20, 10, 0.8, 0, -0.8, -10, -20, -25];
    xgrid = linspace(0, 1, 200);

    % --- Build color map keyed by |f| ---
    mags_all   = abs(f_values);
    mags_uniq  = unique(mags_all, 'stable');   % preserve first-seen order
    nmag       = numel(mags_uniq);
    cmap       = lines(max(nmag, 3));          % get at least a few distinct colors

    % helper to map |f| -> color index (with tolerance for floats)
    tol = 1e-12;
    mag_to_idx = @(m) find(abs(mags_uniq - m) <= tol, 1, 'first');

    % --- Plot ---
    figure; hold on; grid off; box on
    xlabel('$x$','Interpreter','latex'); ylabel('$\bar{u}(x)$','Interpreter','latex');
    title(sprintf('Exact Steady-State Solutions with Two Obstacles $(\\kappa = %.3g)$', kappa), ...
          'Interpreter','latex');

    % compute and plot curves, color by |f|
    Uall = nan(numel(f_values), numel(xgrid));
    plt_handles = gobjects(numel(f_values),1);
    for i = 1:numel(f_values)
        f = f_values(i);
        color_idx = mag_to_idx(abs(f));
        this_color = cmap(color_idx, :);

        [uvals, ~, ~] = eval_u_two_obstacles(xgrid, alpha, kappa, f, y_minus, y_plus);
        Uall(i,:) = uvals;

        plt_handles(i) = plot(xgrid, uvals, 'LineWidth', 1.6, ...
                              'Color', this_color, ...
                              'DisplayName', sprintf('f = %.1f', f));
    end

    % y-limits with margin around obstacles
    yl = ylim;
    tol_ylim = max(1e-6, 0.01*diff(yl));
    ymin = min(yl(1), y_minus - tol_ylim);
    ymax = max(yl(2), y_plus  + tol_ylim);
    ylim([ymin, ymax]);

    % draw the two stops at x=1
    line([1 1], [y_minus, ymin], 'Color','k', 'LineWidth',5, ...
         'LineStyle','-', 'HandleVisibility','off');
    line([1 1], [y_plus,  ymax], 'Color','k', 'LineWidth',5, ...
         'LineStyle','-', 'DisplayName','obstacle');

    legend('Location','northwest');
    hold off;
end

% ================= Helpers (unchanged logic) =================
function [uvals, regime, A, u1] = eval_u_two_obstacles(xvec, alpha, kappa, f, y_minus, y_plus)
    regime = determine_regime(alpha, f, y_minus, y_plus);
    switch regime
        case 'no_contact'
            A = -f/(alpha^2);
        case 'contact_bottom'
            A = a_contact(alpha, kappa, f, y_minus);
        case 'contact_top'
            A = a_contact(alpha, kappa, f, y_plus);
        otherwise
            error('Unknown regime');
    end
    uvals = poly_u(xvec, alpha, A, f);
    u1 = poly_u(1, alpha, A, f);
end

function regime = determine_regime(alpha, f, y_minus, y_plus)
    u_free_tip = f/(8*alpha^2);
    if u_free_tip <= y_minus
        regime = 'contact_bottom';
    elseif u_free_tip >= y_plus
        regime = 'contact_top';
    else
        regime = 'no_contact';
    end
end

function u = poly_u(x, alpha, A, f)
    u = (f/(24*alpha^2))*x.^4 + (A/6)*x.^3 - ((f/(4*alpha^2)) + A/2).*x.^2;
end

function A = a_contact(alpha, kappa, f, y_wall)
    A = - ( kappa*y_wall + f*(1 + (5*kappa)/(24*alpha^2)) ) ...
        / ( alpha^2 + kappa/3 );
end
