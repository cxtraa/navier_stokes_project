%2D Burgers equation
%This combines nonlinear convection and diffusion in 2D

clear;

%Initialising variables
nx = 40; xmax = 2; dx = xmax./(nx-1);
ny = 40; ymax = 2; dy = ymax./(nx-1);
sigma = 0.05; mu = 0.05;
nt = 100; dt = sigma * dx * dy / mu;

x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
[X,Y] = meshgrid(x,y);
t = linspace(0, dt*nt, nt);

%Initial conditions
u = ones(ny, nx, nt);
v = ones(ny, nx, nt);

u(0.5/dy : 1/dy, 0.5/dx : 1/dx, 1) = 6;
v(0.5/dy : 1/dy, 0.5/dx : 1/dx, 1) = 3;

%Solving PDE via iterative method
for n=1:nt
    for j=2:ny-1
        for i=2:nx-1
            u(j,i,n+1) = u(j,i,n) - dt/dx*u(j,i,n)*(u(j,i,n)-u(j,i-1,n)) - dt/dy*v(j,i,n)*(u(j,i,n)-u(j-1,i,n)) + mu*(dt/dx^2)*(u(j,i+1,n)-2*u(j,i,n)+u(j,i-1,n)) + mu*(dt/dy^2)*(u(j+1,i,n)-2*u(j,i,n)+u(j-1,i,n));
            v(j,i,n+1) = v(j,i,n) - dt/dx*u(j,i,n)*(v(j,i,n)-v(j,i-1,n)) - dt/dy*v(j,i,n)*(v(j,i,n)-v(j-1,i,n)) + mu*(dt/dx^2)*(v(j,i+1,n)-2*v(j,i,n)+v(j,i-1,n)) + mu*(dt/dy^2)*(v(j+1,i,n)-2*v(j,i,n)+v(j-1,i,n));

            u(1, :) = 1;
            u(ny, :) = 1;
            u(:, 1) = 1;
            u(:, nx) = 1;

            v(1, :) = 1;
            v(ny, :) = 1;
            v(:, 1) = 1;
            v(:, nx) = 1;
        end
    end
end

figure();
a = quiver(x, y, u(:, :, 1), v(:, :, 1)); grid on; axis tight manual;
set(a, 'AutoScaleFactor', 1.5);
xlabel('X/m'); ylabel('Y/m');
title('Velocity field of fluid: viscous forces only');
for n=2:nt
    set(a, 'udata', u(:,:,n), 'vdata', v(:,:,n));
    pause(0.1);
end

%{
%Display the results
figure();
ax = axes;
hold(ax, 'on');
a = surface(ax, u(:, : , 1)); grid on; axis tight manual;
view(45, 45);
for n=2:nt
    delete(a);
    a = surface(ax, u(:, :, n));
    view(45, 45);
    pause(0.1);
end
%}