%Solution to the heat equation via a finite difference method
%u is used to refer to temperature T

%Inputs
alpha = 1;

nx = 100; xmax = 10;
x = linspace(0, xmax, nx);
dx = xmax/(nx+1);

nt = 10000; tmax = 2;
t = linspace(0, tmax, nt);
dt = tmax/(nt+1);

s = dt/(dx^2);
u = zeros(nt, nx);

%Initial conditions
k = 5*pi/xmax;
for i=1:nx
    u(1, i) = cos(6*pi/xmax * x(i));
end

%Boundary conditions
u(:, 1) = -0.5;
u(:, nx) = 0.5;

%Finite difference method
for n=1:1:nt-1
    for i=2:1:nx-1
        %Derived from derivative approximations to heat equation
        u(n+1, i) = alpha * s * (u(n, i+1) - 2*u(n, i) + u(n, i-1)) + u(n, i);
    end   
end

%Plot surface in 3 dimensions showing x, t, and T
figure();
[X, T] = meshgrid(x, t);
a = surface(X, T, u); grid on;
set(a, 'LineStyle', 'none')
xlabel('x'); ylabel('t'); zlabel('T');
title('Surface showing solution to heat equation');

%Plot evolution of T(x) wrt. time
figure();
l = plot(x, u(1, :)); axis tight manual; grid on;
xlabel('x'); ylabel('T');
for j=1:1:nt-1
    set(l, 'XData', x, 'YData', u(j, :)); title('Evolution of T(x) wrt. time, Time = ', t(j));
    pause(0.0001);
end

