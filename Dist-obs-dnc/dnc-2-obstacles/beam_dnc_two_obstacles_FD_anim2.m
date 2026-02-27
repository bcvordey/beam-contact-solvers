function beam_dnc_two_obstacles_FD_anim2()
% Implicit FD for Euler–Bernoulli beam with two-sided DNC at the tip.
% u_tt + alpha^2 u_xxxx = F0*sin(omega*t),  x in (0,L), t in (0,T)
% Left clamp: u(0,t)=0, u_x(0,t)=0
% Right: u_xx(1,t)=0 and -alpha^2 u_{xxx}(1,t)
%        = kappa[(u(1,t)-y_-)_- - (u(1,t)-y_+)_+] - beta*u_t(1,t)
% (delta = 1). Plots u(1,t) for multiple beta values and animates u(x,t).

    %% Parameters
    L     = 1.0;
    T     = 10.0;
    M     = 20;                  % spatial subintervals (nodes = M+1)
    N     = 1000;                % time steps
    alpha = 1.296;
    kappa = 0.1;                  % contact stiffness
    y_minus = -0.2;              % bottom obstacle
    y_plus  =  0.2;              % top obstacle
    assert(y_minus < y_plus, 'Need y_minus < y_plus');

    % Time-varying body force: F0*sin(omega*t)
    F0    = 10;
    omega = 12.0;                % rad/s

    % Beta sweep  (legend shows only beta values)
    betas = 5;

    % Which beta index to ANIMATE (1..numel(betas))
    anim_beta_idx = 1;

    %% Discretization
    dx = L / M;
    dt = T / N;
    x  = linspace(0, L, M+1);
    t  = linspace(0, T, N+1);
    gamma = alpha^2 * (dt^2) / (dx^4);
    s = (dx^3) / (alpha^2);      % factor in tip equations

    %% Data (phi == 0 at clamp)
    phi   = @(tt) 0.0;
    ffun  = @(xx,tt) F0 * sin(omega * tt);
    u0fun = @(xx) 0*xx;
    v0fun = @(xx) 0*xx;

    %% Run for each beta and collect tip traces
    Utip = zeros(N+1, numel(betas));
    for k = 1:numel(betas)
        beta = betas(k);
        Utip(:,k) = simulate_one_beta(beta, false);  % no full field recording
    end

    %% Plot tip traces and obstacles
    figure; hold on; grid on;
    for k = 1:numel(betas)
        plot(t, Utip(:,k), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('\\beta = %g', betas(k)));
    end
    yline(y_minus, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
    yline(y_plus,  'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
    xlabel('time'); ylabel('u(1,t)');
    title(sprintf('Right-tip displacement for multiple \\beta (\\kappa = %.3g)', kappa));
    legend('Location','best');
    hold off;

    %% Run again for the chosen beta and CAPTURE the full field for animation
    beta_anim = betas(anim_beta_idx);
    [~, U_anim] = simulate_one_beta(beta_anim, true);  % record U(x,t)

    %% ---------------- Animation ----------------
    fig1 = figure('Name', sprintf('Beam with DNC: animation (alpha = %g)', alpha), ...
                  'Color','w');
    % axis limits based on obstacles and solution
    umax = max(abs(U_anim), [], 'all');
    ymax = 1.05 * max([1e-12, umax, abs(y_minus), abs(y_plus)]);
    plt = plot(x, U_anim(:,1), 'b-', 'LineWidth', 2.5); hold on; grid on;

    % Tip marker that switches color on contact
    tipMarker = plot(x(end), U_anim(end,1), 'go', 'MarkerFaceColor','g', ...
                     'DisplayName','Tip');

    % Obstacle lines (dashed) and vertical bars at the tip
    line([x(end), x(end)], [-1.1*ymax, y_minus], 'Color','k', 'LineWidth',5);  % bottom
    line([x(end), x(end)], [y_plus,  1.1*ymax], 'Color','k', 'LineWidth',5);   % top
    yline(y_minus, 'k--', 'LineWidth', 1.2, 'HandleVisibility','off');
    yline(y_plus,  'k--', 'LineWidth', 1.2, 'HandleVisibility','off');

    xlabel('x'); ylabel('u(x,t)');
    % title('Two Boundary Obstacles at x=1');
    title(sprintf('Two Boundary Obstacles at x=1 ($\\kappa = %.3g, \\beta = %.3g$)',kappa,beta), 'Interpreter','latex');
    axis([0, x(end), -ymax, ymax]);

    % --- Animation pacing controls ---
    skip = 3;                       % show every 'skip' frames
    playback_speed = 1000;           % sim-seconds per real-second (larger = faster playback)

    % ---------- GIF saving setup ----------
    gif_filename = 'beam_dnc_two_obstacles.gif';
    frame_delay  = (skip*dt)/max(playback_speed, eps);
    firstFrame   = true;
    % NOTE: we use imwrite to write/append GIF frames (not the 'save' workspace function).

    for j = 1:skip:numel(t)
        set(plt, 'YData', U_anim(:,j));
        set(tipMarker, 'YData', U_anim(end,j));
        % Contact coloring at the tip
        if U_anim(end,j) <= y_minus || U_anim(end,j) >= y_plus
            set(tipMarker, 'MarkerEdgeColor','r','MarkerFaceColor','r');
        else
            set(tipMarker, 'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        drawnow;

        % ---------- Write current frame to GIF ----------
        frame = getframe(fig1);
        [imind, cmap] = rgb2ind(frame2im(frame), 256);
        if firstFrame
            imwrite(imind, cmap, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frame_delay);
            firstFrame = false;
        else
            imwrite(imind, cmap, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
        end

        % Pace the playback for on-screen animation (can be commented out if saving only)
        pause(frame_delay);
    end

    %================ nested solver for a single beta =======================
    % record_full = true returns the full space-time matrix U(x_i, t_j)
    function [u_tip, Ufull] = simulate_one_beta(beta, record_full)
        if nargin < 2, record_full = false; end

        % Initialize at j=0 and j=1  (u^1 = u^0 + dt*v0)
        u_prev = u0fun(x).';
        u_curr = u_prev + dt * v0fun(x).';
        u_prev(1) = phi(t(1));    % = 0
        u_curr(1) = phi(t(2));    % = 0

        % Tip history
        u_tip = zeros(N+1,1);
        u_tip(1) = u_prev(end);
        u_tip(2) = u_curr(end);

        % Optional full field storage for animation
        if record_full
            Ufull = zeros(M+1, N+1);
            Ufull(:,1) = u_prev;
            Ufull(:,2) = u_curr;
        else
            Ufull = [];
        end

        % Preallocate sparse matrix for unknowns (nodes 1..M)
        A = spalloc(M, M, 5*M);

        % Time stepping
        for j = 1:N-1
            tjp1 = t(j+2);  % time at level j+1
            fvec = arrayfun(@(xx) ffun(xx, tjp1), x).';   % size M+1

            % Reset A and b
            A(:) = 0;
            b = zeros(M,1);

            % ---- Row i=1 (ghost u_-1=u_1, and Dirichlet u_0=phi) ----
            A(1,1) = 1 + 7*gamma;
            if M >= 2, A(1,2) = -4*gamma; end
            if M >= 3, A(1,3) =  gamma;   end
            b(1) = dt^2 * fvec(2) + 2*u_curr(2) - u_prev(2) + 4*gamma * phi(tjp1); % phi=0

            % ---- Interior rows: i = 2..M-2 (push col==0 into RHS) ----
            for i = 2:M-2
                row = i;
                cols   = [i-2,  i-1,     i,   i+1,   i+2];
                coeffs = [gamma, -4*gamma, 1+6*gamma, -4*gamma, gamma];

                rhs = dt^2 * fvec(i+1) + 2*u_curr(i+1) - u_prev(i+1);

                for q = 1:5
                    col = cols(q);
                    c   = coeffs(q);
                    if col >= 1 && col <= M
                        A(row, col) = A(row, col) + c;
                    elseif col == 0
                        rhs = rhs + c * phi(tjp1); % phi=0 -> adds 0
                    end
                end
                b(row) = rhs;
            end

            % ---- Row i=M-1 ----
            if M >= 2
                row = M-1;
                A(row, M-3) = A(row, M-3) + gamma;
                A(row, M-2) = A(row, M-2) - 4*gamma;
                A(row, M-1) = A(row, M-1) + (1 + 5*gamma);
                A(row, M  ) = A(row, M  ) - 2*gamma;
                b(row) = dt^2 * fvec(M) + 2*u_curr(M) - u_prev(M);
            end

            % ---- Row i=M (tip law; decide by u_M^j) ----
            row = M;
            u_tip_prev = u_curr(end);  % u_M^j

            if u_tip_prev <= y_minus
                % Bottom contact
                y_wall = y_minus;
                A(row, M-2) = A(row, M-2) + 2*gamma;
                A(row, M-1) = A(row, M-1) - 4*gamma;
                A(row, M  ) = A(row, M  ) + (1 + 2*gamma*(s*(kappa + beta/dt) + 1));

                b(row) = dt^2 * fvec(M+1) ...
                       + 2*(1 + gamma*beta*s/dt)*u_curr(M+1) - u_prev(M+1) ...
                       + 2*gamma*s*kappa*y_wall;

            elseif u_tip_prev >= y_plus
                % Top contact
                y_wall = y_plus;
                A(row, M-2) = A(row, M-2) + 2*gamma;
                A(row, M-1) = A(row, M-1) - 4*gamma;
                A(row, M  ) = A(row, M  ) + (1 + 2*gamma*(s*(kappa + beta/dt) + 1));

                b(row) = dt^2 * fvec(M+1) ...
                       + 2*(1 + gamma*beta*s/dt)*u_curr(M+1) - u_prev(M+1) ...
                       + 2*gamma*s*kappa*y_wall;

            else
                % No contact
                A(row, M-2) = A(row, M-2) + 2*gamma;
                A(row, M-1) = A(row, M-1) - 4*gamma;
                A(row, M  ) = A(row, M  ) + (1 + 2*gamma);

                b(row) = dt^2 * fvec(M+1) + 2*u_curr(M+1) - u_prev(M+1);
            end

            % ---- Solve and roll time levels
            u_next_unknown = A \ b;

            u_next = zeros(M+1,1);
            u_next(1)     = phi(tjp1);      % = 0
            u_next(2:end) = u_next_unknown;

            % store
            if record_full
                Ufull(:, j+2) = u_next;
            end
            u_tip(j+2) = u_next(end);
            u_prev = u_curr;
            u_curr = u_next;
        end
    end
end
