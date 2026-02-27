%% rod_steady_plot.m
% Visualise the steady‑state displacement of a linear rod
% under constant body force f and DNC boundary at x = 0.

clear; clc; close all;

%% --- user parameters ----------------------------------------------------
l      = 1.0;          % natural length of the rod
c2     = 1.296;        % c^2 = E_Y / rho
f      = 0.8;          % body force   (>0)
kappa  = 2.0;          % obstacle stiffness
h      = 0.05;         % obstacle start position
Nx     = 400;          % # evaluation points along the rod
% -------------------------------------------------------------------------

x = linspace(-l,0,Nx).';                            % grid
u = -(f/(2*c2))*x.^2 - (f*l/c2)*x + h + f*l/(kappa*c2);   % formula (34)

figure('Color','w');
plot(x,u,'b-','LineWidth',1.6);  hold on;
yline(h,'k--','LineWidth',1.2);
xlabel('x'); ylabel('\bar u(x)');
title('Steady displacement profile of the rod');
legend('\bar u(x)','obstacle level h','Location','southoutside');
grid on;

%% --- derived quantities -------------------------------------------------
d_star = f*l/(kappa*c2);                 % penetration depth at x = 0
l_st   = l - f*l/(2*c2);                 % contracted length
E_rod  = f^2 * l^3 / (6*c2);             % elastic energy   (rod)
E_obs  = kappa * d_star^2;               % energy in obstacle
E_tot  = E_rod + E_obs;                  % total

fprintf('\nSteady penetration depth  d*  = %.4f\n', d_star);
fprintf('Contracted length         l_st = %.4f (natural l = %.2f)\n', l_st, l);
fprintf('Elastic energy (rod)      = %.4e\n', E_rod);
fprintf('Obstacle energy           = %.4e\n', E_obs);
fprintf('Total steady energy       = %.4e\n\n', E_tot);
