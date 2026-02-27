% obstacle_3D.m — Plot the distributed obstacle ψ(x) in 3D
% ψ(x) = -0.02 + 0.01 sin^2(π x), x ∈ [0,1]

clear; clc;

% --- Define obstacle
% psi = @(x) -0.02 + 0.01 * sin(pi*x).^2;
 % psi  = @(x) -(0.02 + 0.01*cos(50*pi*x));
% psi = @(x) -0.02 - 0.01*sin(pi*x).^2; 
psi = @(x) -0.02 + 0*x;  

% --- Grid along beam (x) and a thin span (y) just for 3D visualization
nx = 600;                 % x-resolution (increase for smoother surface)
ny = 10;                  % y-resolution
x  = linspace(0,1,nx);
y  = linspace(0,0.02,ny); % small width so it's a thin "sheet"
[X,Y] = meshgrid(x,y);

% --- Surface values
Z = psi(X);

% --- Plot obstacle surface
figure('Color','w');
hTop = surf(X,Y,Z, Z, 'EdgeColor','none'); % color by height
shading interp; colormap(parula); hold on;

% ===== Translucent fill under the obstacle =====
% Choose a reference plane to fill down to (adjust if you like)
z_floor = min(Z(:)) - 0.006;    % e.g., slightly below the lowest ψ(x)

% Bottom "floor" surface (optional, gives the fill a base)
hBottom = surf(X, Y, z_floor*ones(size(X)), ...
    'EdgeColor','none', 'FaceColor',[0.3 0.6 0.95], 'FaceAlpha',0.15);

% Side walls to close the shape, so it looks like a solid
% y = min
Xw = X(1,:);  Yw = Y(1,:);  ZwTop = psi(Xw);
surf([Xw; Xw], [Yw; Yw], [z_floor*ones(size(Xw)); ZwTop], ...
    'EdgeColor','none', 'FaceColor',[0.3 0.6 0.95], 'FaceAlpha',0.25);
% y = max
Xw = X(end,:); Yw = Y(end,:); ZwTop = psi(Xw);
surf([Xw; Xw], [Yw; Yw], [z_floor*ones(size(Xw)); ZwTop], ...
    'EdgeColor','none', 'FaceColor',[0.3 0.6 0.95], 'FaceAlpha',0.25);
% x = 0
Xw = X(:,1)'; Yw = Y(:,1)'; ZwTop = psi(Xw);
surf([Xw; Xw], [Yw; Yw], [z_floor*ones(size(Xw)); ZwTop], ...
    'EdgeColor','none', 'FaceColor',[0.3 0.6 0.95], 'FaceAlpha',0.25);
% x = 1
Xw = X(:,end)'; Yw = Y(:,end)'; ZwTop = psi(Xw);
surf([Xw; Xw], [Yw; Yw], [z_floor*ones(size(Xw)); ZwTop], ...
    'EdgeColor','none', 'FaceColor',[0.3 0.6 0.95], 'FaceAlpha',0.25);

% Make the top surface a bit opaque so the fill shows through
set(hTop, 'FaceAlpha', 0.95);

% Optional: centerline
yc = mean(y);
plot3(x, yc*ones(size(x)), psi(x), 'k-', 'LineWidth', 1.5);

% --- Beam: straight line above the obstacle along x ∈ [0,1]
gap    = 0.005;                   % vertical clearance above the highest ψ
z_beam = max(Z(:)) + gap;         % constant z-level for the beam
plot3(x, yc*ones(size(x)), z_beam*ones(size(x)), ...
      'k-', 'LineWidth', 3);      % straight beam line
% --- Clamp indicator at x = 0 (black dot on the beam)
plot3(0, yc, z_beam, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');



% --- Labels & view
xlabel('$x$','Interpreter','latex','FontSize',13);
% ylabel('$y$','Interpreter','latex','FontSize',13);
zlabel('$u(x)$','Interpreter','latex','FontSize',13);
title('Obstacle $\psi(x)=-0.02$ with Clamped Beam at $x=0$', ...
      'Interpreter','latex','FontSize',13);


axis tight; daspect([1 1 0.05]);
view(45,25); grid off; camlight headlight; lighting gouraud;
c = colorbar; c.Label.String = '\psi(x)';

% ===== Save at highest quality =====
outBase = 'obstacle_3D_clamped';

% Set a consistent on-screen size (exportgraphics uses this size)
set(gcf,'Units','inches');
set(gcf,'Position',[1 1 7 4]);   % width x height in inches
drawnow;
% 
% % 1) Vector PDF (best for print). Works if transparency is not critical.
% set(gcf,'Renderer','painters');  % vector renderer
% exportgraphics(gcf, [outBase '.pdf'], ...
%     'ContentType','vector', 'BackgroundColor','white');

% % 2) High-res PNG (best when you use transparency/alpha)
% set(gcf,'Renderer','opengl');    % correct alpha blending
% exportgraphics(gcf, [outBase '.png'], ...
%     'Resolution', 900, 'BackgroundColor','white');

% Optional: SVG (vector) and MATLAB figure for future edits
% exportgraphics(gcf, [outBase '.svg']);
% savefig(gcf, [outBase '.fig']);

% If the PDF shows artifacts due to transparency, rasterize the PDF instead:
exportgraphics(gcf, [outBase '_rasterized.pdf'], ...
    'ContentType','image', 'Resolution', 600, 'BackgroundColor','white');
