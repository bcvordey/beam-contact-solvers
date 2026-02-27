function [x,t,u] = beam_contact_distributed_DNC_FD()
% Implicit second-order FD for Euler–Bernoulli beam with DISTRIBUTED DNC contact (δ=1).
% PDE: u_tt + α^2 u_xxxx = f + r, with
% r^{j+1}_i ≈ χ^j_i [ κ (y_- - u^{j+1}_i) - β (u^{j+1}_i - u^j_i)/dt ],  χ^j_i = 1(u^j_i ≤ y_-).
%
% Left:  u(0,t)=φ(t), u_x(0,t)=0  (ghost u_{-1}=u_1)
% Right: u_xx(1,t)=0, u_xxx(1,t)=0 (free end closures)
%
% Returns: grid (x,t) and solution u(i,j) ≈ u(x_i, t_j)

%% ---------------- Parameters ----------------
L  = 1.0;   T  = 10.0;
dx = 0.01;  dt = 0.01;

alpha2  = 1.296;        % α^2
kappa   = 10;       % κ
beta    = 0.1;        % β
y_minus = -0.02;      % obstacle level (constant)
f0      = -0.2;       % constant forcing: f(x,t)=f0

phi = @(tt) 0.0;      % boundary excitation at x=0
u0f = @(xx) 0.0;      % initial displacement
v0f = @(xx) 0.0;      % initial velocity

%% ---------------- Grids & storage ----------------
M = round(L/dx);                   % nodes i=0..M
N = round(T/dt);                   % times j=0..N
x = linspace(0,L,M+1).';           % column (M+1)x1
t = linspace(0,T,N+1);             % row    1x(N+1)
u = zeros(M+1, N+1);

gamma = alpha2*(dt^2)/(dx^4);

%% ---------------- Initial conditions ----------------
u(:,1) = u0f(x);                   % u^0
u(:,2) = u0f(x) + dt*v0f(x);       % u^1 (2nd-order start)
u(1,:) = phi(t);                   % enforce Dirichlet at x=0 at all used times

%% ---------------- Time stepping ----------------
for j = 1:(N-1)
    % Contact indicator at time level j (boolean at all nodes 0..M)
    I = (u(:, j+1) <= y_minus);
    I(1) = false;  % node 0 is Dirichlet (not solved in the system)

    % Diagonal and RHS additions for DNC (δ=1) at unknown nodes 1..M
    diag_add = dt^2 * (kappa * double(I(2:end)) + (beta/dt) * double(I(2:end)));
    rhs_add  = dt^2 * (kappa * y_minus * double(I(2:end)) ...
                      + (beta/dt) * double(I(2:end)) .* u(2:end, j+1));

    % Assemble A (M x M) and b (M x 1) for [u_1..u_M]^{j+1}
    A = spalloc(M, M, 5*M);
    b = zeros(M,1);

    % ---- Row i=1 (uses u_{-1}=u_1 and u_0=phi)
    A(1,1) = 1 + 7*gamma;
    if M >= 2, A(1,2) = -4*gamma; end
    if M >= 3, A(1,3) =  gamma;   end
    b(1) = dt^2*f0 + 2*u(2, j+1) - u(2, j) + 4*gamma*phi(t(j+1));
    % DNC at i=1
    A(1,1) = A(1,1) + diag_add(1);
    b(1)   = b(1)   + rhs_add(1);

    % ---- Interior rows i=2..M-2 (careful with i-2 hitting the boundary)
    for i = 2:(M-2)
        r = i;
        if i-2 >= 1
            A(r, i-2) = A(r, i-2) + gamma;
        else
            % this would multiply u_0^{j+1} = phi(t_{j+1}) -> move to RHS
            b(r) = b(r) - gamma*phi(t(j+1));
        end
        A(r, i-1) = A(r, i-1) - 4*gamma;
        A(r, i  ) = A(r, i  ) + (1 + 6*gamma);
        A(r, i+1) = A(r, i+1) - 4*gamma;
        A(r, i+2) = A(r, i+2) + gamma;

        b(r) = b(r) + dt^2*f0 + 2*u(i+1, j+1) - u(i+1, j);

        % DNC diagonal/RHS at i
        A(r, i) = A(r, i) + diag_add(i);
        b(r)    = b(r)    + rhs_add(i);
    end

    % ---- Right-end rows (free-end closures)
    if M >= 2
        % i = M-1
        A(M-1, M-3) = A(M-1, M-3) + gamma;
        A(M-1, M-2) = A(M-1, M-2) - 4*gamma;
        A(M-1, M-1) = A(M-1, M-1) + (1 + 5*gamma);
        A(M-1, M  ) = A(M-1, M  ) - 2*gamma;
        b(M-1) = b(M-1) + dt^2*f0 + 2*u(M, j+1) - u(M, j);
        A(M-1, M-1) = A(M-1, M-1) + diag_add(M-1);
        b(M-1)      = b(M-1)      + rhs_add(M-1);
    end

    % i = M
    A(M, M-2) = A(M, M-2) + 2*gamma;
    A(M, M-1) = A(M, M-1) - 4*gamma;
    A(M, M  ) = A(M, M  ) + (1 + 2*gamma);
    b(M) = b(M) + dt^2*f0 + 2*u(M+1, j+1) - u(M+1, j);
    A(M, M) = A(M, M) + diag_add(M);
    b(M)    = b(M)    + rhs_add(M);

    % ---- Solve and write back
    w = A \ b;
    u(1,   j+2) = phi(t(j+2));  % node 0 (Dirichlet)
    u(2:end,j+2) = w;           % nodes 1..M
end

%% ===== Postprocessing & Plots =====
% Derived fields
contact     = (u <= y_minus);
penetration = max(y_minus - u, 0);
frac_contact = mean(contact,1);           % over x
max_pen     = max(penetration,[],1);      % over x

% A small helper for consistent figure styling
function prettify(ax)
    set(ax, 'FontSize', 12, 'LineWidth', 1.0);
    box(ax, 'on'); grid(ax, 'on');
end

% 1) Tip displacement vs time
fig1 = figure('Color','w','Name','Tip displacement'); 
plot(t, u(end,:), 'LineWidth',1.8); hold on;
yline(y_minus,'--','LineWidth',1.4,'DisplayName','y_-');
xlabel('time','Interpreter','latex'); 
ylabel('$u(1,t)$','Interpreter','latex');
title('Right-end displacement','Interpreter','latex');
legend('Location','best','Interpreter','latex');
prettify(gca);

% 2) Space–time contact indicator (white=no, black=yes)
fig2 = figure('Color','w','Name','Contact map');
imagesc(t, x, contact); axis xy tight;
colormap(flipud(gray)); caxis([0 1]); % 0 -> white, 1 -> black
cb = colorbar; set(cb,'Ticks',[0 1],'TickLabels',{'No','Contact'});
xlabel('time','Interpreter','latex'); 
ylabel('$x$','Interpreter','latex');
title('Contact indicator','Interpreter','latex');
prettify(gca);

% 3) Penetration heatmap
fig3 = figure('Color','w','Name','Penetration heatmap');
imagesc(t, x, penetration); axis xy tight;
colormap(parula); cb = colorbar; ylabel(cb,'penetration $(y_- - u)_+$','Interpreter','latex');
xlabel('time','Interpreter','latex'); 
ylabel('$x$','Interpreter','latex');
title('Penetration depth','Interpreter','latex');
prettify(gca);

% 4) 3D displacement surface with obstacle plane
fig4 = figure('Color','w','Name','Displacement surface');
surf(t, x, u, 'EdgeColor','none'); hold on;
surf(t, x, y_minus*ones(size(u)), ...
     'FaceColor','k', 'FaceAlpha',0.5, 'EdgeColor','r');
shading interp; view(45,30); axis tight vis3d;
xlabel('time','Interpreter','latex'); 
ylabel('$x$','Interpreter','latex'); 
zlabel('$u(x,t)$','Interpreter','latex');
title('Beam displacement surface with obstacle','Interpreter','latex');
prettify(gca);

% 5) Contact metrics: fraction of span in contact, and max penetration
fig5 = figure('Color','w','Name','Contact metrics');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(t, frac_contact,'LineWidth',1.8); 
ylabel('fraction in contact','Interpreter','latex');
title('Contact fraction vs time','Interpreter','latex');
prettify(gca);

nexttile;
plot(t, max_pen,'LineWidth',1.8); 
ylabel('max penetration','Interpreter','latex'); 
xlabel('time','Interpreter','latex');
title('Maximum penetration vs time','Interpreter','latex');
prettify(gca);

% exportgraphics(fig2,'contact_map.pdf','Resolution',300);

end
