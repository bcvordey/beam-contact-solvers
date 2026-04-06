function beam_dnc_fem
% -------------------------------------------------------------------------
% Euler–Bernoulli beam FEM (Hermite) with bilateral DNC at the free tip.
% - Cantilever: w=0, θ=0 at x=0
% - Sinusoidal tip force
% - DNC contact at tip (top & bottom obstacles)
% - Implicit Newmark
% - Robust animation (works even if figures are closed & recreated)
% -------------------------------------------------------------------------

%% ----- Physical & model parameters --------------------------------------
L   = 1.0;                 % [m]
b   = 0.02; h = 0.02;      % [m]
E   = 70e9;                % [Pa]
rho = 2700;                % [kg/m^3]
A   = b*h;                 % [m^2]
Izz = b*h^3/12;            % [m^4]

% Contact (DNC) at tip
y_top   = +0.0050;          % [m]
y_bot   = -0.0050;          % [m]
kappa   = 5e6;             % [N/m^(δ)]
beta_c  = 150.0;           % [N·s/m]
delta_c = 1.0;

% External tip force
Q0  = 5.0;
omega1 = (1.8751^2)*sqrt(E*Izz/(rho*A*L^4));
omega  = 0.9*omega1;
forceFun = @(t) Q0*sin(omega*t);

% Rayleigh damping
alphaR = 0.0;
betaR  = 1e-5;

%% ----- Discretization ----------------------------------------------------
Ne  = 40;
Le  = L/Ne;
nn  = Ne + 1;
nd  = 2*nn;                   % DOFs: [w1 θ1 w2 θ2 ... wN θN]
free = 3:nd;                  % clamp node 1
id_tip_w = nd-1;              % tip w DOF

% Assemble
M = sparse(nd,nd);
K = sparse(nd,nd);
for e = 1:Ne
    i1  = 2*(e-1)+1;
    loc = i1:(i1+3);
    [me, ke] = beamEB_hermite_mk(E, Izz, rho, A, Le);
    M(loc,loc) = M(loc,loc) + me;
    K(loc,loc) = K(loc,loc) + ke;
end
C = alphaR*M + betaR*K;

%% ----- Time integration (Newmark) ---------------------------------------
Tend = 5.0;
dt   = 2e-4;
Nt   = round(Tend/dt);

betaN  = 1/4;
gammaN = 1/2;

u = zeros(nd,1);
v = zeros(nd,1);
a = zeros(nd,1);

Aeff_full = M + gammaN*dt*C + betaN*dt^2*K;

tip_hist = zeros(Nt+1,1); tip_hist(1) = u(id_tip_w);
time_vec = (0:Nt)'*dt;

%% ----- Visualization setup ----------------------------------------------
x_nodes    = linspace(0,L,nn);
plot_scale = 1.0;
updateEveryN = max(1, round(0.002/dt));

[fig1, ax1, beam_line, tip_marker, fig2, ax2, tip_plot] = rebuild_figures();

%% ----- Time stepping -----------------------------------------------------
maxFP = 4;          % fixed-point iterations
tolFP = 1e-8;
relax = 0.8;

for n = 1:Nt
    tnp1 = n*dt;

    % Newmark predictors
    u_pred = u + dt*v + (0.5 - betaN)*dt^2*a;
    v_pred = v + (1 - gammaN)*dt*a;

    % External force vector
    Fext = sparse(nd,1);
    Fext(id_tip_w) = forceFun(tnp1);

    % Effective matrix on free set
    Aeff = Aeff_full(free,free);

    % Fixed-point for contact
    a_new = a; u_new = u; v_new = v; wtip_old = u(id_tip_w);
    for it = 1:maxFP
        if it == 1
            wtip = u_pred(id_tip_w);
            wdot = v_pred(id_tip_w);
        else
            wtip = u_new(id_tip_w);
            wdot = v_new(id_tip_w);
        end

        % DNC law at tip
        Fcontact = 0.0;
        if wtip > y_top
            pen = wtip - y_top;
            Fcontact = Fcontact - (kappa * pen^delta_c + beta_c * wdot);
        end
        if wtip < y_bot
            pen = y_bot - wtip;
            Fcontact = Fcontact + (kappa * pen^delta_c - beta_c * wdot);
        end

        Ftotal = Fext;
        Ftotal(id_tip_w) = Ftotal(id_tip_w) + Fcontact;

        rhs_full = Ftotal - C*v_pred - K*u_pred;
        rhs = rhs_full(free);

        a_free = Aeff \ rhs;
        a_new = zeros(nd,1); a_new(free) = a_free;

        u_next = u_pred + betaN*dt^2*a_new;
        v_next = v_pred + gammaN*dt*a_new;

        u_new = relax*u_next + (1-relax)*u_new;
        v_new = relax*v_next + (1-relax)*v_new;

        if abs(u_new(id_tip_w) - wtip_old) < tolFP
            break
        end
        wtip_old = u_new(id_tip_w);
    end

    % accept step
    u = u_new; v = v_new; a = a_new;

    % record tip
    tip_hist(n+1) = u(id_tip_w);

    % UI updates (rebuild if closed)
    if mod(n, updateEveryN)==0 || n==Nt
        if ~is_axes(ax1) || ~ishandle(beam_line) || ~ishandle(tip_marker)
            [fig1, ax1, beam_line, tip_marker] = rebuild_anim([], x_nodes, y_top, y_bot);
        end
        if ~is_axes(ax2) || ~ishandle(tip_plot)
            [fig2, ax2, tip_plot] = rebuild_tipfig(time_vec, tip_hist);
        end

        w_nodes = u(1:2:end);
        set(beam_line,'XData',x_nodes,'YData',plot_scale*w_nodes(:).');
        set(tip_marker,'XData',L,'YData',plot_scale*w_nodes(end));

        set(tip_plot,'XData',time_vec(1:n+1),'YData',tip_hist(1:n+1));
        drawnow limitrate
    end
end

disp('Simulation complete.');

% =================== nested helper functions =============================
    function tf = is_axes(ax)
        tf = ~isempty(ax) && ishandle(ax) && strcmp(get(ax,'Type'),'axes');
    end

    function [fig1, ax1, beam_line, tip_marker, fig2, ax2, tip_plot] = rebuild_figures()
        [fig1, ax1, beam_line, tip_marker] = rebuild_anim([], x_nodes, y_top, y_bot);
        [fig2, ax2, tip_plot]              = rebuild_tipfig(time_vec, tip_hist);
    end

    function [fig1, ax1, beam_line, tip_marker] = rebuild_anim(ax1_in, x_nodes, y_top, y_bot)
        if ~is_axes(ax1_in)
            fig1 = figure('Name','Beam with bilateral DNC: animation','Color','w');
            ax1  = axes('Parent',fig1);
        else
            ax1  = ax1_in;
            fig1 = ancestor(ax1,'figure');
            cla(ax1);
        end

        % Make this axes current, then use classic commands
        axes(ax1); %#ok<LAXES>
        hold on; grid on; box on;
        set(ax1,'XLim',[0 L*1.06]);
        set(ax1,'YLim',1.2*[min(y_bot, -0.012), max(y_top, 0.012)]);
        xlabel('x [m]'); ylabel('w(x,t) [m]');
        title('Beam deflection with bilateral DNC at tip');

        plot([0.95*L, L], [y_top y_top], 'r-', 'LineWidth',2, 'DisplayName','top obstacle');
        plot([0.95*L, L], [y_bot y_bot], 'b-', 'LineWidth',2, 'DisplayName','bottom obstacle');
        beam_line  = plot(x_nodes, zeros(size(x_nodes)), 'k-', 'LineWidth',1.6, 'DisplayName','beam');
        tip_marker = plot(L, 0, 'ko', 'MarkerFaceColor','k', 'DisplayName','tip');
        legend('Location','best');
    end

    function [fig2, ax2, tip_plot] = rebuild_tipfig(time_vec, tip_hist)
        fig2 = figure('Name','Tip displacement','Color','w');
        ax2  = axes('Parent',fig2);
        axes(ax2); %#ok<LAXES>
        hold on; grid on; box on;
        xlabel('time [s]'); ylabel('w(L,t) [m]');
        title('Tip displacement over time');
        tip_plot = plot(time_vec(1), tip_hist(1), 'k-');
    end
end

% -------- Element matrices for Euler–Bernoulli Hermite beam --------------
function [me, ke] = beamEB_hermite_mk(E, Izz, rho, A, Le)
% 4x4 Hermite element (w_i, θ_i, w_j, θ_j)
ke = (E*Izz/Le^3) * ...
    [ 12,    6*Le,   -12,    6*Le;
      6*Le,  4*Le^2, -6*Le,  2*Le^2;
     -12,   -6*Le,    12,   -6*Le;
      6*Le,  2*Le^2, -6*Le,  4*Le^2 ];

me = (rho*A*Le/420) * ...
    [ 156,     22*Le,    54,     -13*Le;
      22*Le,   4*Le^2,   13*Le,  -3*Le^2;
      54,      13*Le,    156,    -22*Le;
     -13*Le,  -3*Le^2,  -22*Le,   4*Le^2 ];
end
