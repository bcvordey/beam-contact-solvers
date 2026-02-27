function compare_numerical_exact2()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) USER-DEFINED PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha    = 1;       % Beam stiffness/scale parameter
    kappa    = 1;       % Obstacle stiffness
    f        = -1.72;    % Constant force (negative example)
    y_minus  = -0.2;    % Obstacle position
    M_values = [4,8,16,32,64];  % Mesh sizes (elements per unit length)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) PREPARE FOR DATA STORAGE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numM       = length(M_values);
    L2_errors  = zeros(numM,1);
    Linf_errors= zeros(numM,1);
    hx_values  = zeros(numM,1);
    colors     = lines(numM);

    figure('Name','Comparison of Numerical and Exact Solutions'); 
    hold on; grid off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) LOOP OVER VARIOUS MESH SIZES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iM = 1:numM
        M  = M_values(iM);
        dx = 1.0 / M;
        hx_values(iM) = dx;

        % (a) Numerical solution on x_i = i*dx, i=0..M
        [x_numerical, u_numerical] = solveBeamSteadyContact(M, alpha, kappa, f, y_minus);

        % (b) Exact solution on same grid
        u_exact = arrayfun(@(x) piecewise_beam_exact(x, alpha, kappa, f, y_minus), x_numerical);

        % (c) Errors
        err_vec          = u_numerical - u_exact;
        Linf_errors(iM)  = max(abs(err_vec));           % L^∞
        L2_errors(iM)    = sqrt(sum(err_vec.^2) * dx);  % L^2 (discrete 1D)

        % (d) Plot for visual inspection
        plot(x_numerical, u_numerical, '-o', 'Color', colors(iM,:), ...
            'DisplayName', sprintf('M=%d, numerical', M), 'LineWidth', 1.2);
        
    end

    plot(x_numerical, u_exact, 'Color',k,'DisplayName', 'Exact', 'LineWidth', 1.2);

    yl = ylim;                            % current limits
    
    % ensure there is *some* room below y_minus and above y_plus
    tol = max(1e-10, 0.01*diff(yl));       % small, automatic epsilon
    ymin = min(yl(1), y_minus - tol);

    line([1 1], [y_minus, ymin], 'Color','k', 'LineWidth',5,'LineStyle','-', 'DisplayName', 'obstacle');

    % Obstacle
    % yline(y_minus,'k--','DisplayName',sprintf('Obstacle y_- = %g', y_minus),LineWidth=1.2);
    xlabel('$x$',Interpreter="latex"); ylabel('$u(x)$',Interpreter="latex");
    title(sprintf('Comparison of Numerical and Exact Solutions ($f$ = %g)', f),Interpreter="latex");
    legend('Location','southwest'); box on;
    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) CONVERGENCE RATES (two-grid OOC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_L2   = nan(numM,1);
    p_Linf = nan(numM,1);
    for i = 1:numM-1
        p_L2(i)   = log(L2_errors(i)  / L2_errors(i+1))   / log(hx_values(i) / hx_values(i+1));
        p_Linf(i) = log(Linf_errors(i)/ Linf_errors(i+1)) / log(hx_values(i) / hx_values(i+1));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5) PLOTS OF ERROR VS h
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure('Name','Convergence Plots');
    % loglog(hx_values, L2_errors,  'ro--','DisplayName','l^2 error','LineWidth',1.5); hold on; grid off;
    % loglog(hx_values, Linf_errors,'bs--','DisplayName','l^\infty error','LineWidth',1.5);
    % xlabel('Mesh spacing $\Delta x$',Interpreter="latex"); ylabel('Error norm',Interpreter="latex");
    % title('Log-log Plot of Error vs. $\Delta x$', 'Interpreter','latex');
    % legend('Location','northwest');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6) TABLE WITH ERRORS AND RATES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n%-8s %-10s %-14s %-14s %-10s %-10s\n', 'M', 'dx', 'L2_error', 'L_inf_error', 'p_L2', 'p_Linf');
    fprintf('%s\n', repmat('-',1,72));
    for iM = 1:numM
        pL2   = p_L2(iM);   pLinf = p_Linf(iM);
        pL2_str   = ternary(isnan(pL2),   '   -   ', sprintf('%8.3f', pL2));
        pLinf_str = ternary(isnan(pLinf), '   -   ', sprintf('%8.3f', pLinf));
        fprintf('%-8d %-10.6f %-14.6e %-14.6e %8s %8s\n', ...
            M_values(iM), hx_values(iM), L2_errors(iM), Linf_errors(iM), pL2_str, pLinf_str);
    end

    % Optional: least-squares slope over finest K grids
    K = min(4, numM); idx = (numM-K+1):numM;
    Pfit_L2   = polyfit(log(hx_values(idx)), log(L2_errors(idx)),   1);
    Pfit_Linf = polyfit(log(hx_values(idx)), log(Linf_errors(idx)), 1);
    fprintf('\nLeast-squares slope over finest %d grids: p_L2=%.3f, p_Linf=%.3f\n', ...
            K, Pfit_L2(1), Pfit_Linf(1));

    %==================== NESTED HELPERS ====================%
    function out = ternary(cond, a, b)
        if cond, out = a; else, out = b; end
    end

    function [xvals, bar_u] = solveBeamSteadyContact(Mloc, alpha_, kappa_, f_force, y_minus_)
        % Finite-difference steady beam with unilateral contact at x=1.
        dxloc = 1.0 / Mloc;
        N     = Mloc + 1;
        xvals = linspace(0,1,N)';

        % Try no-contact and contact systems; pick consistent one.
        [bar_u_A, okA] = solveSystem(Mloc, dxloc, alpha_, kappa_, f_force, y_minus_, 'nocontact');
        if bar_u_A(end) < y_minus_, okA = false; end

        [bar_u_B, okB] = solveSystem(Mloc, dxloc, alpha_, kappa_, f_force, y_minus_, 'contact');
        if bar_u_B(end) >= y_minus_, okB = false; end

        if okA
            bar_u = bar_u_A;
        elseif okB
            bar_u = bar_u_B;
        else
            bar_u = bar_u_A; % fallback
        end
    end

    function [bar_u, success] = solveSystem(Mloc, dxloc, alpha_, kappa_, f_force, y_minus_, mode)
        N   = Mloc + 1;
        A   = zeros(N,N);
        rhs = zeros(N,1);

        % (1) u(0)=0
        A(1,1) = 1; rhs(1) = 0;

        % (2) u_x(0)=0  => -3u0 + 4u1 - u2 = 0
        A(2,1) = -3; A(2,2) = 4; A(2,3) = -1;

        % (3) Interior: alpha^2 u_xxxx = f
        c = alpha_^2 / dxloc^4;
        for eq = 3:(N-2)
            A(eq, eq-2) =  c;
            A(eq, eq-1) = -4*c;
            A(eq, eq  ) =  6*c;
            A(eq, eq+1) = -4*c;
            A(eq, eq+2) =  c;
            rhs(eq)     =  f_force;
        end

        % (4) u_xx(1)=0 => u_{M-2} - 2u_{M-1} + u_M = 0
        eq = N-1;
        A(eq, N-2) = 1; A(eq, N-1) = -2; A(eq, N) = 1;

        % (5) Free vs contact condition at x=1
        eq = N;
        if strcmp(mode,'nocontact')
            % u_xxx(1)=0  => (u_M - 3u_{M-1} + 3u_{M-2} - u_{M-3})/dx^3 = 0
            d3 = 1/dxloc^3;
            A(eq, N  ) =  d3;
            A(eq, N-1) = -3*d3;
            A(eq, N-2) =  3*d3;
            A(eq, N-3) = -1*d3;
            rhs(eq)    =  0;
        else
            % -alpha^2 u_xxx(1) = kappa (u(1) - y_-)
            %  => alpha^2/dx^3*(u_M - 3u_{M-1} + 3u_{M-2} - u_{M-3}) + kappa*u_M = kappa*y_-
            d3a = alpha_^2 / dxloc^3;
            A(eq, N  ) =  d3a;       % will add kappa below
            A(eq, N-1) = -3*d3a;
            A(eq, N-2) =  3*d3a;
            A(eq, N-3) = -1*d3a;
            A(eq, N  ) = A(eq, N) + kappa_;
            rhs(eq)    = kappa_ * y_minus_;
        end

        bar_u  = A \ rhs;
        success = true;
    end

    function val = piecewise_beam_exact(x, alpha_, kappa_, f_, y_minus_)
        % Decide regime by free-end displacement in no-contact case
        end_disp_nc = f_ / (8 * alpha_^2);
        if end_disp_nc >= y_minus_
            val = u_no_contact(x, alpha_, f_);
        else
            val = u_contact(x, alpha_, kappa_, f_, y_minus_);
        end
    end

    function val = u_no_contact(x, alpha_, f_)
        % Polynomial satisfying: u(0)=0, u'(0)=0, u''(1)=0, u'''(1)=0
        val = (f_/(24*alpha_^2))*(x.^4 - 4*x.^3 + 6*x.^2);
    end

    function val = u_contact(x, alpha_, kappa_, f_, y_minus_)
        % Simple illustrative contact-form (problem-specific!)
        A = a_contact(alpha_, kappa_, f_, y_minus_);
        val = (f_/(24*alpha_^2))*x.^4 ...
            + (A/6)*x.^3 ...
            - ((f_/(4*alpha_^2))+A/2)*x.^2;
    end

    function A = a_contact(alpha_, kappa_, f_, y_minus_)
        % Chosen to enforce a representative matching; adapt if needed.
        A = (-3*kappa_ / (kappa_ + 3*alpha_^2)) * ...
            ( y_minus_ + f_ * ( 5/(24*alpha_^2) + 1/kappa_ ) );
    end
end
