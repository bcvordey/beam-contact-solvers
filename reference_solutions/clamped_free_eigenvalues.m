%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Clamped-Free Beam Eigenvalue Problem
% and plot cos(lambda)*cosh(lambda) + 1 vs. lambda
%
% Euler-Bernoulli PDE:  rho A d^2 w/dt^2 + E I d^4 w/dx^4 = 0
%
% BCs: 
%   x=0 (clamped):    X(0)=0, X'(0)=0
%   x=L (free):       X''(L)=0, X'''(L)=0
%
% Characteristic eqn: cos(lambda)*cosh(lambda) + 1 = 0
%
% by [Your Name], [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% Beam and material parameters (example values)
E   = 1;     % Young's modulus (Pa)
rho = 1;      % density (kg/m^3)
A   = 1;    % cross-section area (m^2)
I   = 1; % second moment of area (m^4)
L   = 1.0;       % length of the beam (m)

% We want the first N modes
N   = 4;

%% Define the characteristic function: f(lambda) = cos(lambda)*cosh(lambda) + 1
fchar = @(lambda) cos(lambda).*cosh(lambda) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1: Numerically find the first N roots using fzero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Known approximate root locations for cos(lambda)*cosh(lambda) + 1 = 0
% (We can use these as initial guesses in fzero)
lambdaApprox = [1.8751, 4.6941, 7.8548, 10.9955, 14.1372, 17.2788];

% Preallocate
lambdaVals = zeros(1,N);
omegaVals  = zeros(1,N);

%% Solve for roots numerically using fzero
for n = 1:N
    guess = lambdaApprox(n);        % initial guess
    lambdaVals(n) = fzero(fchar, guess);
    
    % Once we have lambda, compute the natural frequency (rad/s):
    %   mu = lambda / L
    %   mu^4 = (rho*A/EI)*omega^2  -->  omega^2 = (EI/rho*A)*mu^4 = ...
    %   => omega = (lambda^2 / L^2) * sqrt(EI/(rho*A))
    lambda_n = lambdaVals(n);
    omegaVals(n)  = (lambda_n^2 / L^2) * sqrt(E*I/(rho*A));
end

%% Print out results
disp('  Mode   lambda_n       omega_n (rad/s)');
for n = 1:N
    fprintf('  %2d     %8.5f      %10.5f\n', n, lambdaVals(n), omegaVals(n));
end

%% Plot the frequencies vs. mode number
figure;
stem(1:N, omegaVals, 'LineWidth', 2);
xlabel('Mode number n');
ylabel('\omega_n (rad/s)');
title('Clamped-Free Beam Natural Frequencies');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2: Plot the characteristic function f(lambda) itself
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose a range of lambda values in which to see the zero crossings
lambdaRange = linspace(0, 1.5, 1000);

% Evaluate f(lambda) = cos(lambda)*cosh(lambda) + 1
fVals = fchar(lambdaRange);

% Plot
figure;
plot(lambdaRange, fVals,'b', 'LineWidth', 1.3); 
hold on;
yline(0, 'r*', LineWidth=1.3);  % reference line y=0

xlabel('\lambda');
ylabel('cos(\lambda) cosh(\lambda) + 1');
title('Characteristic function for the clamped-free beam');
grid on;

% Optionally, mark found roots on the plot
for n = 1:N
    plot(lambdaVals(n), 0, 'ro', 'MarkerSize', 8, ...
         'DisplayName', sprintf('\\lambda_{%d} \\approx %.3f', n, lambdaVals(n)));
end
legend('f(\lambda)','y=0','Roots','Location','best');
hold off;
