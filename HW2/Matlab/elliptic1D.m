% Solves the 1D Poisson's equation with Robin boundary conditions

clc
close all

addpath('../mole_MATLAB')

% Domain's limits
west = 0; 
east = 1;

k = 4;  % Operator's order of accuracy %NEW 6->4
m = 2*k+2;  % Minimum number of cells to attain the desired accuracy %NEW 1->2
dx = (east-west)/m;  % Step length

L = lap(k, m, dx);  % 1D Mimetic laplacian operator

% Impose Robin BC on laplacian operator
a = 1; 
b = 1;
L = L + robinBC(k, m, dx, a, b);

% 1D Staggered grid
grid = [west west+dx/2 : dx : east-dx/2 east];
grid

% RHS
% U = exp(grid)'; 
U = 1+(grid)'; % NEW exp(grid)' -> 1+(grid)'
U(1) = 0;  % West BC
U(end) = 2*exp(1);  % East BC

U = L\U;  % Solve a linear system of equations
U

anal = ((12*exp(1)-13)/18).*(1+(grid))+(((grid).^2)./2)+(((grid).^3)./6); % analytical solution
anal

% Plot result
plot(grid, U, 'o')
hold on
plot(grid, anal)
legend('Approximated', 'Analytical', 'Location', 'NorthWest')
title('Poisson''s equation with Robin BC')
xlabel('x')
ylabel('u(x)')
