%Navier-stokes solver: cavity flow boundary conditions

clear;

%Setting up 2D space and time
nx = 50; xmax = 2; dx = xmax/(nx-1);
ny = 50; ymax = 2; dy = ymax/(ny-1);
nt = 500; dt = 0.001; 

x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
t = linspace(0, dt*nt, nt);
[X,Y] = meshgrid(x,y);

%r = density, nu = viscosity, iters = no. of iterations to calculate
%converged pressure field
r = 1;
nu = 0.1;
iters = 50;

u = zeros(ny, nx, nt);
v = zeros(ny, nx, nt);
p = zeros(ny, nx, nt);

%We iterate through time
for n = 1:nt+1

    %Calculate the pressure field by repeating this process iters times
    %before converging
    for k=1:iters
        for i=2:nx-1
            for j=2:ny-1
                p(j,i,n+1) = ((p(j,i+1,n) + p(j,i-1,n))*(dy^2) + (p(j+1,i,n) + p(j-1,i,n))*(dx^2))/(2*((dx^2)+(dy^2))) - (r*dx^2*dy^2)/((2*((dx^2)+(dy^2)))) * (1/dt * ((u(j,i+1,n) - u(j,i-1,n))/(2*dx) + (v(j+1,i,n) - v(j-1,i,n))/(2*dy)) - ((u(j,i+1,n) - u(j,i-1,n))/(2*dx))^2 - 2*((u(j+1,i,n) - u(j-1,i,n))/(2*dy))*((v(j,i+1,n) - v(j,i-1,n))/(2*dx)) - ((v(j+1,i,n) - v(j-1,i,n))/(2*dy))^2);
            end
        end
    end
    
    %Pressure field boundary conditions
    p(:, 1, n+1) = 2;
    p(:, nx, n+1) = -2;
    p(1, :, n+1) = p(2, :, n+1); 
    p(ny, :, n+1) = p(ny-1, :, n+1);
    %p(:, 1) = p(:, 2); 
    %p(:, nx, n+1) = p(:, nx-1, n+1);
    
    %Iterating through space to calculate the velocity field, using our
    %calculated pressure field which ensures continuity is satisfied
    for j=2:ny-1
        for i=2:nx-1
            u(j,i,n+1) = u(j,i,n) - u(j,i,n)*(dt/dx)*(u(j,i,n) - u(j,i-1,n)) - v(j,i,n)*(dt/dy)*(u(j,i,n) - u(j-1,i,n)) - (dt/(2*r*dx))*(p(j,i+1,n) - p(j,i-1,n)) + nu*(dt/(dx^2)*(u(j,i+1,n) - 2*u(j,i,n) + u(j,i-1,n)) + dt/(dy^2)*(u(j+1,i,n) - 2*u(j,i,n) + u(j-1,i,n)));
            v(j,i,n+1) = v(j,i,n) - u(j,i,n)*(dt/dx)*(v(j,i,n) - v(j,i-1,n)) - v(j,i,n)*(dt/dy)*(v(j,i,n) - v(j-1,i,n)) - (dt/(2*r*dy))*(p(j+1,i,n) - p(j-1,i,n)) + nu*(dt/(dx^2)*(v(j,i+1,n) - 2*v(j,i,n) + v(j,i-1,n)) + dt/(dy^2)*(v(j+1,i,n) - 2*v(j,i,n) + v(j-1,i,n)));
        end
    end
    
    %Velocity field boundary conditions
    u(ny, :, n+1) = 0;
    u(:, 1, n+1) = 3;
    u(:, nx, n+1) = 0;
    u(1, :, n+1) = 0;

    v(ny, :, n+1) = 0;
    v(:, 1, n+1) = 0;
    v(:, nx, n+1) = 0;
    v(1, :, n+1) = 0;
end

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
    pause(0.01);
end

        



