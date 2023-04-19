%Diffusion equation in 1D using CFL condition

clear all;

nx = 500; xmax = 5;
x = linspace(0, xmax, nx)
dx = xmax/(nx-1);
nt = 10000;
sigma = 0.1; nu  = 5;
dt = sigma * dx^2;
t = linspace(0, dt*nt, nt);

u = zeros(nt, nx);
u(1, :) = sin(3*pi/xmax * x).^2;
u(1, 1) = 2; u(1, nx) = -0.5;

for j = 1:1:nt-1
    for i = 2:1:nx-1
        u(j+1, i) = u(j, i) + nu * dt/(dx.^2) * (u(j, i+1) - 2*u(j, i) + u(j, i-1));
    end

    %Boundary conditions
    u(j+1, 1) = 2; 
    u(j+1, nx) = -0.5;
end

l = plot(x, u(1, :)); grid on; axis tight manual;
for n = 2:1:nt
    set(l, 'YData', u(n, :));
    pause(0.001);
end

[X, T] = meshgrid(x, t);
a = surface(X, T, u);
set(a, 'LineStyle', 'none');
