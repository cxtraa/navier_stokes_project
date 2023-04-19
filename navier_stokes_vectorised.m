%Navier-stokes solver: cavity flow boundary conditions

clear;

%Setting up 2D space and time
nx = 41; xmax = 2; dx = xmax/(nx-1);
ny = 41; ymax = 2; dy = ymax/(ny-1);
nt = 500; dt = 0.001; 

x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
t = linspace(0, dt*nt, nt);
[X,Y] = meshgrid(x,y);

%r = density, nu = viscosity, iters = no. of iterations to calculate
%converged pressure field
r = 1;
nu = 0.6;
iters = 20;

u = zeros(ny, nx, nt);
v = zeros(ny, nx, nt);
p = zeros(ny, nx, nt);

%We iterate through time
for n = 1:nt+1

    %Calculate the pressure field by repeating this process iters times
    %before converging
    for k=1:iters        
        p(2:ny-1, 2:nx-1, n+1) = ((p(2:ny-1, 3:nx, n) + p(2:ny-1, 1:nx-2, n)).*(dy.^2)...
            + (p(3:ny, 2:nx-1, n) + p(1:ny-2, 2:nx-1, n)).*(dx.^2))/(2*((dx^2)+(dy.^2)))...
            - (r*dx.^2*dy.^2)/((2*((dx.^2)+(dy.^2)))) .* (1/dt .* ((u(2:ny-1, 3:nx, n)...
            - u(2:ny-1, 1:nx-2, n))/(2*dx) + (v(3:ny, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n))/(2*dy))...
            - ((u(2:ny-1, 3:nx, n) - u(2:ny-1, 1:nx-2, n))/(2*dx)).^2 - 2.*((u(3:ny, 2:nx-1, n)...
            - u(1:ny-2, 2:nx-1, n))/(2*dy)).*((v(2:ny-1, 3:nx, n) - v(2:ny-1, 1:nx-2, n))/(2*dx))...
            - ((v(3:ny, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n))/(2*dy)).^2);
        
        %Pressure field boundary conditions
        p(ny, :, n+1) = 0; %P = 0 @y=ymax
        p(1, :, n+1) = p(2, :, n+1); %dP/dy = 0 @y=0
        p(:, 1, n+1) = p(:, 2, n+1); %dP/dx = 0 @x=0
        p(:, nx, n+1) = p(:, nx-1, n+1); %dP/dx = 0 @x=xmax
    end          
    
    %Iterating through space to calculate the velocity field, using our
    %calculated pressure field which ensures continuity is satisfied
    u(2:ny-1, 2:nx-1, n+1) = u(2:ny-1, 2:nx-1, n) - u(2:ny-1, 2:nx-1, n).*(dt/dx).*(u(2:ny-1, 2:nx-1, n)...
        - u(2:ny-1, 1:nx-2, n)) - v(2:ny-1, 2:nx-1, n).*(dt/dy).*(u(2:ny-1, 2:nx-1, n) - u(1:ny-2, 2:nx-1, n))...
        - (dt/(2*r*dx)).*(p(2:ny-1, 3:nx, n) - p(2:ny-1, 1:nx-2, n)) + nu*(dt/(dx^2)*(u(2:ny-1, 3:nx, n)...
        - 2*u(2:ny-1, 2:nx-1, n) + u(2:ny-1, 1:nx-2, n)) + dt/(dy^2).*(u(3:ny, 2:nx-1, n) - 2*u(2:ny-1, 2:nx-1, n)...
        + u(1:ny-2, 2:nx-1, n)));

    v(2:ny-1, 2:nx-1, n+1) = v(2:ny-1, 2:nx-1, n) - u(2:ny-1, 2:nx-1, n).*(dt/dx).*(v(2:ny-1, 2:nx-1, n)...
        - v(2:ny-1, 1:nx-2, n)) - v(2:ny-1, 2:nx-1, n).*(dt/dy).*(v(2:ny-1, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n))...
        - (dt/(2*r*dy)).*(p(3:ny, 2:nx-1, n) - p(1:ny-2, 2:nx-1, n)) + nu*(dt/(dx^2)*(v(2:ny-1, 3:nx, n)...
        - 2*v(2:ny-1, 2:nx-1, n) + v(2:ny-1, 1:nx-2, n)) + dt/(dy^2).*(v(3:ny, 2:nx-1, n) - 2*v(2:ny-1, 2:nx-1, n)...
        + v(1:ny-2, 2:nx-1, n)));
    
    %Velocity field boundary conditions
    u(ny, :, n+1) = 5;
    u(:, 1, n+1) = 0;
    u(:, nx, n+1) = 0;
    u(1, :, n+1) = 0;    

    v(ny, :, n+1) = 0;
    v(:, 1, n+1) = 0;
    v(:, nx, n+1) = 0;
    v(1, :, n+1) = 0;
end

%{
%Display results as a vector field for the fluid velocity, showing the
%evolution of the velocity field in time
figure();
a = quiver(x, y, u(:,:,1), v(:,:,1)); axis tight manual; grid on;
xlabel('X'); ylabel('Y');
title('Solution of Navier-Stokes equations given ICs and BCs');
%b = contourf(X, Y, p(:,:,1));
for n=1:nt
    set(a, 'UData', u(:, :, n), 'VData', v(:, :, n));
    %set(b, 'ZData', p(:,:,n));
    pause(0.03);
end
%}

quiver(x, y, u(:, :, nt), v(:, :, nt), 'LineWidth', 0.5, 'AutoScaleFactor', 2); grid on; axis tight manual; axis equal;
xlabel('X'); ylabel('Y');
title('Solution of Navier-Stokes equations given ICs and BCs');
      



