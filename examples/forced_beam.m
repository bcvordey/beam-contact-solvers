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

% % Special Parameters
beta = 1;
kappa = 0.1;
obstacle = -0.05;
Tol = 1e-16;
force = 1;


% aluminum beam with rectangular cross section W=20cm, H=5cm
A=1; % m^2
rho=1; % kg/m .   % Density of the beam (kg/m)
E=1; % Pa .   
I=1; % m^4

alpha=((rho*A)/(E*I))^0.5;

disp('lambda must be less than or equal 1 for stability')
lambda = dt/((dx^2)*(alpha));



tau = 1/alpha;
phi = dt*tau -2*beta*dx^3 +2*kappa*dx^3*dt;

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
        u(i,j)= -(lambda^2)*(u(i-2,j-1)-4*u(i-1,j-1)-4*u(i+1,j-1)+u(i+2,j-1)+6*u(i,j-1))-u(i,j-2)+2*u(i,j-1);   
    end

    % Boundary conditions
    u(M-1,j)=2*u(M-2,j)-u(M-3,j);
    u(M,j)=3*u(M-2,j)-2*u(M-3,j);
      
end  

%% visualization
% % Create a VideoWriter object
% video_filename = 'beam_vibration.mp4'; % Specify the video filename
% v = VideoWriter(video_filename); % Choose 'MPEG-4' for .mp4 format
% v.FrameRate = 10; % Set the frame rate (adjust as needed)
% open(v); % Open the video file for writing

h=plot(x,u(:,1),'-','LineWidth',5);
xlabel('x');
ylabel('displacements u(x,t)');
title('vibrations of a beam with Obstacle')
axis([0 1.1*L -1.1*ymax 1.1*ymax]);

for kk=2:500:N
    set(h,'YData',u(:,kk));
    drawnow; 
    % % Capture the current figure as a frame and write to the video
    % frame = getframe(gcf); % Capture the figure window
    % writeVideo(v, frame); % Write the frame to the video
    % 
end

% % Close the video file
% close(gca)
% close(v);
% 
% disp(['Video saved as ', video_filename]);