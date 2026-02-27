function beam_contact_distributed_DNC_compare_betas()
% Compare the distributed DNC beam for multiple beta values on the SAME plots.
% Line plots (tip, metrics) are overlaid; image/3D plots are tiled by beta.
%
% Example:
%   beam_contact_distributed_DNC_compare_betas([0.02 0.05 0.1 0.2])

%% ---- Parameters (edit as needed) ----
L  = 1.0;   T  = 10.0;
dx = 0.01;  dt = 0.01;

alpha2  = 1.0;            % α^2
kappa   = 1e-2;           % κ
y_minus = -0.02;          % obstacle level (constant)
f0      = -0.2;           % constant forcing f(x,t)=f0
betaVec = [1e-3, 1e-2, 1e-1, 1, 5];

phi  = @(tt) 0.0;         % boundary excitation at x=0
u0f  = @(xx) 0.0;         % initial displacement
v0f  = @(xx) 0.0;         % initial velocity

% Solve for each beta and store results
nb   = numel(betaVec);
sol  = cell(nb,1);
for k = 1:nb
    sol{k} = solve_single_beta(betaVec(k), L,T,dx,dt, alpha2,kappa,y_minus,f0, phi,u0f,v0f);
end

% Shared grids (identical across beta runs)
x = sol{1}.x;  t = sol{1}.t;

%% ====== Overlaid line plots ======
cols = lines(nb);

% 1) Tip displacement vs time
fig1 = figure('Color','w','Name','Tip displacement vs time'); hold on;
for k = 1:nb
    plot(t, sol{k}.u(end,:), 'LineWidth', 1.8, 'Color', cols(k,:));
end
yline(y_minus,'--','LineWidth',1.4,'DisplayName','y_-');
xlabel('time','Interpreter','latex'); ylabel('$u(1,t)$','Interpreter','latex');
title('Right-end displacement for multiple $\beta$','Interpreter','latex');
legend([compose('$\\beta=%.3g$', betaVec), "y_-"], 'Interpreter','latex', 'Location','best');
prettify(gca);

% 5) Contact metrics (overlaid)
fig5 = figure('Color','w','Name','Contact metrics (multi-\beta)');
tl = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
for k = 1:nb, plot(t, sol{k}.frac_contact, 'LineWidth',1.8, 'Color', cols(k,:)); end
ylabel('fraction in contact','Interpreter','latex');
title('Contact fraction vs time','Interpreter','latex');
legend(compose('$\\beta=%.3g$', betaVec), 'Interpreter','latex', 'Location','best');
prettify(gca);

nexttile; hold on;
for k = 1:nb, plot(t, sol{k}.max_pen, 'LineWidth',1.8, 'Color', cols(k,:)); end
ylabel('max penetration','Interpreter','latex'); xlabel('time','Interpreter','latex');
title('Maximum penetration vs time','Interpreter','latex');
prettify(gca);

%% ====== Tiled image/surface plots (one panel per beta) ======

% 2) Contact maps
fig2 = figure('Color','w','Name','Contact map per \beta');
nrow = ceil(sqrt(nb)); ncol = ceil(nb/nrow);
tl2 = tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
for k = 1:nb
    nexttile;
    imagesc(t, x, sol{k}.contact); axis xy tight;
    colormap(gca, flipud(gray)); caxis([0 1]);
    title(sprintf('Contact: \\beta = %.3g', betaVec(k)), 'Interpreter','tex');
    xlabel('time'); ylabel('x'); prettify(gca);
end
cb = colorbar; cb.Layout.Tile = 'east'; cb.Ticks=[0 1]; cb.TickLabels={'No','Contact'};

% 3) Penetration heatmaps
fig3 = figure('Color','w','Name','Penetration heatmap per \beta');
tl3 = tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
for k = 1:nb
    nexttile;
    imagesc(t, x, sol{k}.penetration); axis xy tight;
    colormap(gca, parula);
    title(sprintf('Penetration: \\beta = %.3g', betaVec(k)), 'Interpreter','tex');
    xlabel('time'); ylabel('x'); prettify(gca);
end
cb = colorbar; cb.Layout.Tile = 'east';
ylabel(cb,'penetration $(y_- - u)_+$','Interpreter','latex');

% 4) 3D surfaces
fig4 = figure('Color','w','Name','Displacement surface per \beta');
tl4 = tiledlayout(nrow,ncol,'Padding','compact','TileSpacing','compact');
for k = 1:nb
    nexttile;
    surf(t, x, sol{k}.u, 'EdgeColor','none'); hold on;
    surf(t, x, y_minus*ones(size(sol{k}.u)), ...
         'FaceColor',[0 0 0], 'FaceAlpha',0.12, 'EdgeColor','none');
    shading interp; view(45,30); axis tight vis3d;
    title(sprintf('\\beta = %.3g', betaVec(k)), 'Interpreter','tex');
    xlabel('time'); ylabel('x'); zlabel('u'); prettify(gca);
end

end  % compare function

% ------------ Helper: solve for a single beta ------------
function S = solve_single_beta(beta, L,T,dx,dt, alpha2,kappa,y_minus,f0, phi,u0f,v0f)
M = round(L/dx);  N = round(T/dt);
x = linspace(0,L,M+1).';  t = linspace(0,T,N+1);
u = zeros(M+1, N+1);
gamma = alpha2*(dt^2)/(dx^4);

% Initial data (2nd-order start)
u(:,1) = u0f(x);
u(:,2) = u0f(x) + dt*v0f(x);
u(1,:) = phi(t);  % Dirichlet

for j = 1:(N-1)
    I = (u(:, j+1) <= y_minus);  I(1) = false;
    diag_add = dt^2 * (kappa * double(I(2:end)) + (beta/dt) * double(I(2:end)));
    rhs_add  = dt^2 * (kappa * y_minus * double(I(2:end)) ...
                      + (beta/dt) * double(I(2:end)) .* u(2:end, j+1));

    A = spalloc(M, M, 5*M);  b = zeros(M,1);

    % Row i=1
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));
    A(1,1) = A(1,1) + diag_add(1);  b(1) = b(1) + rhs_add(1);

    % Interior i=2..M-2
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            b(r) = b(r) - gamma*phi(t(j+1));
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);

        A(r, i) = A(r, i) + diag_add(i);
        b(r)    = b(r)    + rhs_add(i);
    end

    % Right end i=M-1
    if M >= 2
        A(M-1, M-3) = A(M-1, M-3) + gamma;
        A(M-1, M-2) = A(M-1, M-2) - 4*gamma;
        A(M-1, M-1) = A(M-1, M-1) + (1 + 5*gamma);
        A(M-1, M  ) = A(M-1, M  ) - 2*gamma;
        b(M-1) = b(M-1) + dt^2*f0 + 2*u(M, j+1) - u(M, j);
        A(M-1, M-1) = A(M-1, M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % Right end i=M
    A(M, M-2) = A(M, M-2) + 2*gamma;
    A(M, M-1) = A(M, M-1) - 4*gamma;
    A(M, M  ) = A(M, M  ) + (1 + 2*gamma);
    b(M) = b(M) + dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    A(M, M) = A(M, M) + diag_add(M);
    b(M)    = b(M)    + rhs_add(M);

    % Solve and write back
    w = A \ b;
    u(1,   j+2) = phi(t(j+2));
    u(2:end,j+2) = w;
end

% Derived fields
contact     = (u <= y_minus);
penetration = max(y_minus - u, 0);
frac_contact = mean(contact,1);
max_pen     = max(penetration,[],1);

% Pack
S.x=x; S.t=t; S.u=u;
S.contact=contact; S.penetration=penetration;
S.frac_contact=frac_contact; S.max_pen=max_pen;
S.beta = beta;
end

% ------------ Simple plot styling helper ------------
function prettify(ax)
set(ax, 'FontSize', 12, 'LineWidth', 1.0);
box(ax, 'on'); grid(ax, 'on');
end
