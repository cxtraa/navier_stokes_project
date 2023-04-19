%Solving the Laplace equation using an iterative method
%d^2p/dx^2 + d^2p/dy^2 = 0
%There is no time dependence - we will try to find the p(j, i) that
%satisifes this equation by letting it run until it reaches equilbrium
%So once the difference between the current iteration of p and the next is
%very small we know it has almost reached equilbrium so we will stop

clear;

%Initialising variables
nx = 31; xmax = 2; dx = xmax/(nx-1);
ny = 31; ymax = 2; dy = ymax/(nx-1);

x = linspace(0, xmax, nx);
y = linspace(0, 2, ny);
[X,Y] = meshgrid(x, y);

%Initialising p and initial conditions
p = zeros(ny, nx);
p(1, :) = x;
p(ny, :) = xmax - x;
p(:, 1) = y;
p(:, nx) = ymax - y;

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
            p(j, i) = ((dy^2)*(po(j,i+1) + po(j,i-1)) + (dx^2)*(po(j+1,i) + po(j-1,i)))/(2 * (dx^2 + dy^2));            
        end
    end

    %Boundary conditions
    p(1, :) = x;
    p(ny, :) = xmax - x;
    p(:, 1) = y;
    p(:, nx) = ymax - y;
    
    %Fractional difference between p (new) and po (old)
    iter_dif = sum(abs(p(:)) - abs(po(:))) / sum(abs(po(:)));
end

%Plot the equilbrium solution
surf(X, Y, p); grid on;
xlabel('X'); ylabel('Y'); zlabel('p');
title('Solution to Laplace equation for given BCs and ICs');
