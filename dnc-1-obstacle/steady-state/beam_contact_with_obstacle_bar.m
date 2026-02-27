function beam_contact_with_obstacle_bar()
    %==============================================================
    %  Beam bending with unilateral contact at x=1
    %  Fourth-order ODE: alpha^2 u''''(x) = f,   0 < x < 1
    %
    %  Boundary conditions:
    %     u(0)   = 0
    %     u'(0)  = 0
    %     u''(1) = 0
    %     u'''(1)= -kappa * (u(1) - y_minus)_+
    %
    %  This code also draws a thick black vertical bar at x=1
    %  from y=0 to y=y_minus, representing the obstacle.
    %==============================================================

    % ------------------
    % Parameters
    % ------------------
    alpha   = 1.0;    % beam stiffness scale
    kappa   = 2;    % contact stiffness
    y_minus = -0.2;   % obstacle location below the x-axis
    % You can set y_minus > 0 if obstacle is above the x-axis.

    % Forcing values (some negative, some positive)
    f_values = [-3.0, -2.0, -1.6, -1.0, -0.5, 0];

    % Discretize x in [0, 1] for plotting
    xgrid = linspace(0, 1, 200);

    % ------------------
    % Make figure
    % ------------------
    figure('Name','Beam Contact with Obstacle Bar','NumberTitle','off');
    hold on;  grid on;
    xlabel('x');
    ylabel('u(x)');
    title('Beam displacement with obstacle at x=1');

    % ------------------------------------------------
    % 1) Plot the vertical obstacle bar at x=1
    %    from (1,0) down to (1,y_minus).
    % ------------------------------------------------
    % If y_minus < 0, this line goes below the x-axis.
    % If y_minus > 0, it extends above the x-axis.
    line([1, 1], [y_minus,-0.4], 'Color', 'k', 'LineWidth', 2);

    % ------------------------------------------------
    % 2) Plot the piecewise beam solution for each f
    % ------------------------------------------------
    for f = f_values
        % Evaluate the piecewise solution (no-contact/contact)
        uvals = arrayfun(@(x) piecewise_u(x, alpha, kappa, f, y_minus), xgrid);
        plot(xgrid, uvals, 'DisplayName', sprintf('f = %.1f', f));
    end

    legend('Location','best');
    hold off;
end


%=================================================
%    Piecewise solution: no-contact or contact
%=================================================
function val = piecewise_u(x, alpha, kappa, f, y_minus)
    regime = check_contact_consistency(alpha, kappa, f, y_minus);
    if strcmp(regime, 'no_contact')
        val = u_no_contact(x, alpha, f);
    else
        val = u_contact(x, alpha, kappa, f, y_minus);
    end
end

%=================================================
%   Decide if the solution is no-contact or contact
%=================================================
function regime = check_contact_consistency(alpha, kappa, f, y_minus)
    % Threshold for switching from no-contact to contact
    threshold = 8 * alpha^2 * y_minus;
    
    if f <= threshold
        % Expect no contact
        regime = 'no_contact';
    else
        % Expect contact; check consistency:
        uc = u_contact(1.0, alpha, kappa, f, y_minus);
        if uc > y_minus
            regime = 'contact';
        else
            regime = 'no_contact';
        end
    end
end

%=================================
%   No-contact (free-end) solution
%=================================
function val = u_no_contact(x, alpha, f)
    % From the derivation: 
    %    u(x) = (f/(24 alpha^2)) [6 x^2 - 4 x^3 + x^4]
    val = (f/(24*alpha^2)) * (6*x.^2 - 4*x.^3 + x.^4);
end

%=================================
%   Contact solution
%=================================
function val = u_contact(x, alpha, kappa, f, y_minus)
    A3 = a3_contact(alpha, kappa, f, y_minus);
    A2 = a2_contact(alpha, f, A3);
    val = A2*x.^2 + A3*x.^3 + (f/(24*alpha^2))*x.^4;
end

%==============================
%   Polynomial coefficients
%==============================
function A3 = a3_contact(alpha, kappa, f, y_minus)
    % a3 = [kappa*y_minus + f(5kappa -24)/(24 alpha^2)] / (6 - 2kappa)
    numerator   = kappa*y_minus + (f*(5*kappa - 24)) / (24*alpha^2);
    denominator = 6 - 2*kappa;
    A3 = numerator / denominator;
end

function A2 = a2_contact(alpha, f, A3)
    % a2 = -3 a3 - f/(4 alpha^2)
    A2 = -3*A3 - f/(4*alpha^2);
end
