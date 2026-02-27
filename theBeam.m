clc;close all;clear;

% This code solves the free vibration of a 1D beam using Finite Difference 

%% User Data
L=1; % total length of the beam (m)
dx=0.025; %step size for the length (m)
x=0:dx:L; % locations of nodes
M=length(x); % number of nodes

T=2; % final computation time
dt=0.000004; % time increment
t=0:dt:T; % time points
N=length(t); % number of time steps used

animationSpeed = 500;
timeIntervalAdjust = 10;

% aluminum beam with rectangular cross section W=20cm, H=5cm
A=0.5*0.2; % m^2
rho=2700; % kg/m
E=69e10; % Pa
I=0.2*(0.05^3)/12; % m^4

alpha=((rho*A)/(E*I))^0.5;

disp('lambda must be less than or equal 1 for stability')
lambda = dt/((dx^2)*(alpha));

%% Boundary Conditions and Initial Condition

u=zeros(M,N); % intialization of displacements for all nodes at all time points

% Boundary conditions
u(1,:)=0;  
u(2,:)=u(1,:); 

% initial conditions
P=3*E*I*0.1/(L^3); % force for a 10-cm upward displacement at the end
u(:,1)=P*(x.^2).*(3*L-x)/(6*E*I); % initial displacements at all nodes, q(x)
u(:,2)=u(:,1);
ymax=max(abs(u(:,1)));

%% Solution

for j=3:N   % loop over time points
    for i=3:M-2  % loop over nodes
        u(i,j)=-(lambda^2)*(u(i-2,j-1)-4*u(i-1,j-1)-4*u(i+1,j-1)+u(i+2,j-1)+6*u(i,j-1))-u(i,j-2)+2*u(i,j-1);        
    end
    % Boundary conditions
    u(M-1,j)=2*u(M-2,j)-u(M-3,j);
    u(M,j)=3*u(M-2,j)-2*u(M-3,j);
end    


%% visualization

h=plot(x,u(:,1),'-','LineWidth',5);
xlabel('x');
ylabel('displacement u(x,t)');
title('free vibrations of a beam')
axis([0 1.1*L -1.1*ymax 1.1*ymax]);

for kk=2:animationSpeed:N
    set(h,'YData',u(:,kk));
    drawnow; 
end

% --- Sub-Sample Time and Data ---
idx        = 1:timeIntervalAdjust:length(t);
t_sub      = t(idx);
dist0M     = u(M,:);
dist0M_sub = dist0M(idx);

% --- Create Figure & Plot ---
fig1 = figure;
plot(t_sub, dist0M_sub, ...
     '-', ...                    % a line with markers at each sub-sample
     'LineWidth', 0.5, ...        % thicker line for clarity
     'MarkerSize', 1, ...         % small markers
     'DisplayName', 'Position of the tip of the beam');
hold on;                          % so we can add more plot elements easily

% --- Axis Labels & Title ---
xlabel('Time (s)');
ylabel('u(M,t)');
title({'Free end of the Beam vs. Time',sprintf('\\alpha=%.4f',alpha)});

% --- Cosmetics & Legend ---

legend('Location','northeast');
hold off;
% filename_fft = sprintf('ext_cantilever_tip_%.2g.jpg', alpha);
% saveas(fig1, filename_fft);  % Save the figure 
% 
% 
% %---------FFT-------%
% %-------------------%
% 
% uFFT = fft(dist0M_sub);
% fig_fft = figure;
% plot(abs(uFFT),'LineWidth', 1,'MarkerSize', 1.5)
% xlabel('Frequency (Hz)')
% ylabel('|FFT[u(M,t)]|')
% title({'FFT of Cantilever Beam without an Obstacle',sprintf('\\alpha=%.4f',alpha)});
% 
% filename_fft = sprintf('ext_cantilever_fft_%.2g.jpg', alpha);
% saveas(fig_fft, filename_fft);  % Save the figure  

% %--------SHIFT FFT -----%
% uFFTshift = fftshift(uFFT);
% 
% figure;
% 
% plot(abs(uFFTshift))
% title({'FFT SHIFT of Cantilever Beam'});
% xlabel('Frequency (Hz)')
% ylabel('|FFTshift[u(M,t)]|')