%% Clamped–Free Euler–Bernoulli beam: eigenvalues and characteristic function
% PDE:   rho*A*w_tt + E*I*w_xxxx = 0
% BCs:   X(0)=0, X'(0)=0 (clamped);  X''(L)=0, X'''(L)=0 (free)
% Char. equation:  f(lambda) = cos(lambda)*cosh(lambda) + 1 = 0

clear; clc;

%% (Optional) material/geometry (only used for omega)
E   = 1;      % Young's modulus
rho = 1;      % density
A   = 1;      % area
I   = 1;      % second moment
L   = 1;      % length

N = 6;        % number of eigenvalues wanted

% Characteristic function (vectorized)
fchar = @(lam) cos(lam).*cosh(lam) + 1;

%% Find roots by bracketing sign changes on a grid
lam_max = 20;                 % search up to this lambda
lam_grid = linspace(0, lam_max, 20001);
f_grid = fchar(lam_grid);

% Detect sign changes between consecutive points
sgn = sign(f_grid);
idx = find(sgn(1:end-1).*sgn(2:end) < 0);   % indices where a sign change occurs

lambdaVals = zeros(1,N);
k = 0;
for m = 1:length(idx)
    if k >= N, break; end
    a = lam_grid(idx(m));
    b = lam_grid(idx(m)+1);
    % Safety: ensure the bracket really changes sign
    if fchar(a)*fchar(b) > 0, continue; end
    k = k + 1;
    lambdaVals(k) = fzero(fchar, [a b]);   % root in [a,b]
end

if k < N
    warning('Only found %d roots up to lambda = %.1f. Increase lam_max.', k, lam_max);
    lambdaVals = lambdaVals(1:k);
    N = k;
end

%% Natural frequencies (optional)
omegaVals = (lambdaVals.^2 / L^2) * sqrt(E*I/(rho*A));

fprintf('  Mode   lambda_n         omega_n (rad/s)\n');
for n = 1:N
    fprintf('  %2d     %9.6f      %12.6f\n', n, lambdaVals(n), omegaVals(n));
end

%% Plot f(lambda) with roots marked
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; hold on; box on; grid on;
plot(lam_grid, f_grid, 'b-', 'LineWidth', 1.3);
yline(0,'k--','LineWidth',1.0);

% Mark roots
plot(lambdaVals, zeros(size(lambdaVals)), 'ro', 'MarkerSize', 7, ...
    'DisplayName', 'Roots');

% Annotate each root
for n = 1:N
    text(lambdaVals(n), 0, sprintf('  $\\lambda_{%d}=%.3f$', n, lambdaVals(n)), ...
        'VerticalAlignment','bottom','Interpreter','latex');
end

xlabel('$\lambda$');grid off;
ylabel('$f(\lambda)=\cos\lambda\,\cosh\lambda + 1$');
title('Characteristic Function and Clamped-Free Eigenvalues',Interpreter='latex');
legend({'$f(\lambda)$','$0$','$Roots$'},'Location','best');
hold off;

%% (Optional) stem plot of frequencies
% figure; box on; grid off;
% stem(1:N, omegaVals, 'filled','LineWidth',1.3);
% xlabel('Mode $n$');
% ylabel('$\omega_n$ (rad/s)');
% title('Natural frequencies for clamped--free beam (optional)');
