%Navier-stokes solver: cavity flow boundary conditions

clear;

%Setting up 2D space and time
nx = 20; xmax = 2; dx = xmax/(nx-1);
ny = 20; ymax = 2; dy = ymax/(ny-1);
nt = 300; dt = 0.001; 

x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
t = linspace(0, dt*nt, nt);
[X,Y] = meshgrid(x,y);

%r = density, nu = viscosity, iters = no. of iterations to calculate
%converged pressure field
r = 1;
nu = 2;
iters = 50;
F = 1;

u = zeros(ny, nx, nt);
v = zeros(ny, nx, nt);
p = zeros(ny, nx, nt);

%We iterate through time
for n = 1:nt+1

    %Calculate the pressure field by repeating this process iters times
    %before converging
    for k=1:iters        
        p(2:ny-1, 2:nx-1, n+1) = ((p(2:ny-1, 3:nx, n) + p(2:ny-1, 1:nx-2, n)).*(dy.^2) + (p(3:ny, 2:nx-1, n) + p(1:ny-2, 2:nx-1, n)).*(dx.^2))/(2*((dx^2)+(dy.^2))) - (r*dx.^2*dy.^2)/((2*((dx.^2)+(dy.^2)))) .* (1/dt .* ((u(2:ny-1, 3:nx, n) - u(2:ny-1, 1:nx-2, n))/(2*dx) + (v(3:ny, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n))/(2*dy)) - ((u(2:ny-1, 3:nx, n) - u(2:ny-1, 1:nx-2, n))/(2*dx)).^2 - 2.*((u(3:ny, 2:nx-1, n) - u(1:ny-2, 2:nx-1, n))/(2*dy)).*((v(2:ny-1, 3:nx, n) - v(2:ny-1, 1:nx-2, n))/(2*dx)) - ((v(3:ny, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n))/(2*dy)).^2);
        
        %Periodic boundary conditions
        p(2:ny-1, nx, n+1) = ((p(2:ny-1, 1, n) + p(2:ny-1, nx-1, n)).*(dy.^2) + (p(3:ny, nx, n) + p(1:ny-2, nx, n)).*(dx.^2))/(2*((dx^2)+(dy.^2))) - (r*dx.^2*dy.^2)/((2*((dx.^2)+(dy.^2)))) .* (1/dt .* ((u(2:ny-1, nx, n) - u(2:ny-1, nx-1, n))/(2*dx) + (v(3:ny, nx, n) - v(1:ny-2, nx, n))/(2*dy)) - ((u(2:ny-1, 1, n) - u(2:ny-1, nx-1, n))/(2*dx)).^2 - 2.*((u(3:ny, nx, n) - u(1:ny-2, nx, n))/(2*dy)).*((v(2:ny-1, 1, n) - v(2:ny-1, nx-1, n))/(2*dx)) - ((v(3:ny, nx, n) - v(1:ny-2, nx, n))/(2*dy)).^2);
        p(2:ny-1, 1, n+1) = ((p(2:ny-1, 2, n) + p(2:ny-1, nx, n)).*(dy.^2) + (p(3:ny, 1, n) + p(1:ny-2, 1, n)).*(dx.^2))/(2*((dx^2)+(dy.^2))) - (r*dx.^2*dy.^2)/((2*((dx.^2)+(dy.^2)))) .* (1/dt .* ((u(2:ny-1, 1, n) - u(2:ny-1, nx, n))/(2*dx) + (v(3:ny, 1, n) - v(1:ny-2, 1, n))/(2*dy)) - ((u(2:ny-1, 2, n) - u(2:ny-1, nx, n))/(2*dx)).^2 - 2.*((u(3:ny, 1, n) - u(1:ny-2, 1, n))/(2*dy)).*((v(2:ny-1, 2, n) - v(2:ny-1, nx, n))/(2*dx)) - ((v(3:ny, 1, n) - v(1:ny-2, 1, n))/(2*dy)).^2);

        %Pressure field boundary conditions
        p(1, :, n+1) = p(2, :, n+1);
        p(ny, :, n+1) = p(ny-1, :, n+1);
    end          
    
    %Iterating through space to calculate the velocity field, using our
    %calculated pressure field which ensures continuity is satisfied
    u(2:ny-1, 2:nx-1, n+1) = u(2:ny-1, 2:nx-1, n) - u(2:ny-1, 2:nx-1, n).*(dt/dx).*(u(2:ny-1, 2:nx-1, n) - u(2:ny-1, 1:nx-2, n)) - v(2:ny-1, 2:nx-1, n).*(dt/dy).*(u(2:ny-1, 2:nx-1, n) - u(1:ny-2, 2:nx-1, n)) - (dt/(2*r*dx)).*(p(2:ny-1, 3:nx, n) - p(2:ny-1, 1:nx-2, n)) + nu*(dt/(dx^2)*(u(2:ny-1, 3:nx, n) - 2*u(2:ny-1, 2:nx-1, n) + u(2:ny-1, 1:nx-2, n)) + dt/(dy^2).*(u(3:ny, 2:nx-1, n) - 2*u(2:ny-1, 2:nx-1, n) + u(1:ny-2, 2:nx-1, n)) + F*dt);
    v(2:ny-1, 2:nx-1, n+1) = v(2:ny-1, 2:nx-1, n) - u(2:ny-1, 2:nx-1, n).*(dt/dx).*(v(2:ny-1, 2:nx-1, n) - v(2:ny-1, 1:nx-2, n)) - v(2:ny-1, 2:nx-1, n).*(dt/dy).*(v(2:ny-1, 2:nx-1, n) - v(1:ny-2, 2:nx-1, n)) - (dt/(2*r*dy)).*(p(3:ny, 2:nx-1, n) - p(1:ny-2, 2:nx-1, n)) + nu*(dt/(dx^2)*(v(2:ny-1, 3:nx, n) - 2*v(2:ny-1, 2:nx-1, n) + v(2:ny-1, 1:nx-2, n)) + dt/(dy^2).*(v(3:ny, 2:nx-1, n) - 2*v(2:ny-1, 2:nx-1, n) + v(1:ny-2, 2:nx-1, n)));

    %Periodic boundary conditions at x = 0, x = 2
    u(2:ny-1, 1, n+1) = u(2:ny-1, 1, n) - u(2:ny-1, 1, n).*(dt/dx).*(u(2:ny-1, 1, n) - u(2:ny-1, nx, n)) - v(2:ny-1, 1, n).*(dt/dy).*(u(2:ny-1, 1, n) - u(1:ny-2, 1, n)) - (dt/(2*r*dx)).*(p(2:ny-1, 2, n) - p(2:ny-1, nx, n)) + nu*(dt/(dx^2)*(u(2:ny-1, 2, n) - 2*u(2:ny-1, 1, n) + u(2:ny-1, nx, n)) + dt/(dy^2).*(u(3:ny, 1, n) - 2*u(2:ny-1, 1, n) + u(1:ny-2, 1, n)) + F*dt);
    u(2:ny-1, nx, n+1) = u(2:ny-1, nx, n) - u(2:ny-1, nx, n).*(dt/dx).*(u(2:ny-1, nx, n) - u(2:ny-1, nx-1, n)) - v(2:ny-1, nx, n).*(dt/dy).*(u(2:ny-1, nx, n) - u(1:ny-2, nx, n)) - (dt/(2*r*dx)).*(p(2:ny-1, 1, n) - p(2:ny-1, nx-1, n)) + nu*(dt/(dx^2)*(u(2:ny-1, 1, n) - 2*u(2:ny-1, nx, n) + u(2:ny-1, nx-1, n)) + dt/(dy^2).*(u(3:ny, nx, n) - 2*u(2:ny-1, nx, n) + u(1:ny-2, nx, n)) + F*dt);

    %Periodic boundary conditions at x = 0, x = 2
    v(2:ny-1, 1, n+1) = v(2:ny-1, 1, n) - u(2:ny-1, 1, n).*(dt/dx).*(v(2:ny-1, 1, n) - v(2:ny-1, nx, n)) - v(2:ny-1, 1, n).*(dt/dy).*(v(2:ny-1, 1, n) - v(1:ny-2, 1, n)) - (dt/(2*r*dy)).*(p(3:ny, 1, n) - p(1:ny-2, 1, n)) + nu*(dt/(dx^2)*(v(2:ny-1, 2, n) - 2*v(2:ny-1, 1, n) + v(2:ny-1, nx, n)) + dt/(dy^2).*(v(3:ny, 1, n) - 2*v(2:ny-1, 1, n) + v(1:ny-2, 1, n)));
    v(2:ny-1, nx, n+1) = v(2:ny-1, nx, n) - u(2:ny-1, nx, n).*(dt/dx).*(v(2:ny-1, nx, n) - v(2:ny-1, nx-1, n)) - v(2:ny-1, nx, n).*(dt/dy).*(v(2:ny-1, nx, n) - v(1:ny-2, nx, n)) - (dt/(2*r*dy)).*(p(3:ny, nx, n) - p(1:ny-2, nx, n)) + nu*(dt/(dx^2)*(v(2:ny-1, 1, n) - 2*v(2:ny-1, nx, n) + v(2:ny-1, nx-1, n)) + dt/(dy^2).*(v(3:ny, nx, n) - 2*v(2:ny-1, nx, n) + v(1:ny-2, nx, n)));    
    
    u(1, :, n+1) = 0;
    u(ny, :, n+1) = 0;
    v(1, :, n+1) = 0;
    v(ny, :, n+1) = 0;
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
    set(a, 'UData', u(:,:,n), 'VData', v(:,:,n));
    %set(b, 'ZData', p(:,:,n));
    pause(0.00001);
end
%}

quiver(x, y, u(:, :, nt), v(:, :, nt)); grid on; axis tight manual; axis equal;
xlabel('X'); ylabel('Y');
title('Solution of Navier-Stokes equations given ICs and BCs');




