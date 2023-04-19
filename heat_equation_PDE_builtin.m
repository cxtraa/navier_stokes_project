%1D heat equation 2nd order PDE numerical solution
%Situation: 1D rod of length L with boundary conditions of T = T_L on left
%boundary and T = T_R on right boundary. Initial condition is some
%temperature distribution T(x, 0)

L = 2;
x = linspace(0, L, 500); %x coordinate along rod
t = linspace(0, 5, 1000);
[X, T] = meshgrid(x, t); %Creating mesh of x and t to use later
m = 0;
solution = pdepe(m, @PDEs, @IC, @BC, x, t);

%Producing a surface plot. x axis represents x, y axis represents time, z
%axis represents temperature. Shows evolution of T(x) wrt. time
figure;
a = surf(X, T, solution); %Producing a 3D surface
set(a, 'LineStyle', 'none')
title('1D heat equation solution over time')

%Producing a 2D plot that animates T(x) - it changes T(x) with each
%timestep to the new T(x) at that time
figure;
l = plot(x, solution(1, :)); hold on; axis manual; grid on;
xlabel('x/m'); ylabel('T/K');
for i = 1:length(t)
    title('Heat equation, Time = ', t(i))
    set(l, 'XData', x, 'YData', solution(i, :));
    pause(0.1);
end

%Specifying the PDE equation: du/dt = alpha * d^2u / dx^2
%c, f, and s are coefficients in order to 'fill in' MATLAB's PDE template
%found online
function [c, f, s] = PDEs(x, t, T, dTdx, alpha)
    alpha = 1e-2;
    c = 1;
    f = alpha * dTdx;
    s = 0;
end

%The initial condition. T0 is the temperature distribution T(x) at t = 0
function T0 = IC(x, L)
    L = 2;
    T0 = sin(5 * pi/L * x).^2;
end

%The boundary conditions. This sets T(0, t) = 0.5 and T(L, t) = 0.5
function [pl, ql, pr, qr] = BC(xl, Tl, xr, Tr, t)
    pl = Tl - 0.5;
    ql = 0;
    pr = Tr - 3;
    qr = 0;
end



