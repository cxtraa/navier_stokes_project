%Solving the Poisson equation using v. similar code to Laplace equation

clear;

%Initialising variables
nx = 50; xmax = 1; dx = xmax/(nx-1);
ny = 50; ymax = 1; dy = ymax/(nx-1);

x = linspace(0, xmax, nx);
y = linspace(0, 2, ny);
[X,Y] = meshgrid(x, y);

%Initialising p and initial conditions
p = zeros(ny, nx);
b = zeros(ny, nx);

b(int8(1/4 * ny), int8(1/4 * nx)) = 100;
b(int8(3/4 * ny), int8(3/4 * nx)) = -100;

%Assigning a dummy number to iter_dif
iter_dif = 1;

%We will stop once the fractional difference between p_current and p_next
%is less than 1e-4
iter_dif_target = 1e-4;

while iter_dif > iter_dif_target
    po = p; %Current p(j,i) which we will use to calculate the next p(j,i)

    %Iterative loop solving for p(j,i) from discretisation of Laplace
    %equation
    for i=2:nx-1
        for j=2:ny-1
            p(j, i) = ((dy^2)*(po(j,i+1) + po(j,i-1)) + (dx^2)*(po(j+1,i) + po(j-1,i)) - b(j,i)*(dx^2)*(dy^2))/(2 * (dx^2 + dy^2));            
        end
    end

    %Boundary conditions
    p(:, 1) = 0;
    p(1, :) = 0;
    p(ny, :) = 0;
    p(:, nx) = 0;
    
    %Fractional difference between p (new) and po (old)
    iter_dif = sum(abs(p(:)) - abs(po(:))) / sum(abs(po(:)));
end

%Plot the equilbrium solution
surf(X, Y, p); grid on;
xlabel('X'); ylabel('Y'); zlabel('p');
title('Solution to Poisson equation for given BCs and ICs');
