% Solves the 2D Heat equation, changed example wave2D/parabolic1D
% Thomas Keller, COMP670 HW3
clc
close all

addpath('../mole_MATLAB')

% Spatial discretization
k = 2;   % Order of accuracy
m = 3*k+1;  % Number of cells along the x-axis
n = m;   % Number of cells along the y-axis
a = 0;   % West
b = 1;   % East
c = 0;   % South
d = 1;   % North
dx = (b-a)/m;  % Step length along the x-axis
dy = (d-c)/n;  % Step length along the y-axis

% 2D Staggered grid
xgrid = [a a+dx/2 : dx : b-dx/2 b];
ygrid = [c c+dy/2 : dy : d-dy/2 d];

% Create 2D meshgrid
[X, Y] = meshgrid(xgrid, ygrid);

% Mimetic operator Laplacian 2D
L = lap2D(k, m, dx, n, dy);

% alpha
alpha = 1;

% Check neumann stability
% Neumann stability criterion
dt = dx^2/(5*alpha);

% BC
U = zeros(m+2, n+2);
U(1,:)=100;
U(:,1)=100;
U(end,:)=100;
U(:,end)=100;
U

U  = reshape(U, (m+2)*(n+2), 1);
U

% Simulation time
TIME = 0.3;

explicit = 1; % Explicit, for 0 -> implicit

if explicit
    % Explicit
    L = alpha*dt*L + speye(size(L)); % S = speye returns a sparse scalar 1 to L

    % Time integration loop
    for t = 1 : TIME/dt+1
        % Plot result explicit
        mesh(X, Y, reshape(U, m+2, n+2))
        title(['Parabolic 2D, Explicit Solution, Time = ' num2str(dt*t, '%1.2f')])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        caxis([0, 100]) % change colors of plot
        axis([0 1 0 1 0 101])
        drawnow
        U = L*U; % Apply the operator
    end
else
    % Implicit
    L = -alpha*dt*L + speye(size(L));
    
    % Time integration loop
    for t = 1 : TIME/dt+1
        % Plot result implicit
        mesh(X, Y, reshape(U, m+2, n+2))
        title(['Parabolic 2D, Implicit Solution, Time = ' num2str(dt*t, '%1.2f')])
        xlabel('x')
        ylabel('y')
        zlabel('z')
        colorbar
        caxis([0, 100]) % change colors of plot
        axis([0 1 0 1 0 101])
        drawnow
        U = L\U; % Solve a linear system of equations (unconditionally stable)
    end
end