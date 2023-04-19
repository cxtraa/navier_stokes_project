%Burgers equation in 1D - this combines the diffusion and convection
%equations from earlier

clear;

nx = 101; xmax = 2*pi;
xx = linspace(0, xmax, nx);
dx = xmax/(nx-1);
nt = 1000;
sigma = 0.1; nu  = 0.05;
dt = dx * nu;
t = linspace(0, dt*nt, nt);

%Sawtooth function initial condition
%{
syms phi(x) u(x)
phi(x) = exp(-x^2 / (4*nu)) + exp(-(x-2*pi)^2 / (4*nu));
u(x) = -(2*nu)/phi(x) * diff(phi, x) + 4;
ufunc = matlabFunction(u(x));
%}

u = zeros(nt, nx);
u(1, :) = sin(5 * pi/xmax * xx) + 2;

for j = 1:1:nt-1
    for i = 2:1:nx-1
        u(j+1, i) = u(j, i) - (u(j, i) * dt/dx * (u(j, i) - u(j, i-1))) + (nu * dt/(dx.^2) * (u(j, i+1) - 2*u(j, i) + u(j, i-1)));
    end
    
    %Periodic boundary condition - the last point on the function u(x)
    %becomes the first point so the wave loops around on the plot
    u(j+1, 1) = u(j, 1) - (u(j, 1) * dt/dx * (u(j, 1) - u(j, nx-1))) + (nu * dt/(dx.^2) * (u(j, 2) - 2*u(j, 1) + u(j, nx-1)));
    u(j+1, nx) = u(j+1, 1);
end

l = plot(xx, u(1, :)); grid on; axis tight manual;
for n = 2:1:nt
    set(l, 'YData', u(n, :));
    pause(0.01);
end

figure();
[X, T] = meshgrid(xx, t);
a = surface(X, T, u);
set(a, 'LineStyle', 'none');
