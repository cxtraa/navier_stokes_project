%Numerical solution to a PDE for 1D convection
%du/dt + c * du/dx = 0
%u = wave height, t = time, x = x coordinate, c = velocity

clear all;

nx = 100; %Number of x points
x = linspace(0, 5, nx); 
dx = 5./(nx-1); %xmax divided by number of points gives dx
nt = 1000; 
sigma = 0.1;
dt = sigma * dx;
t = linspace(0, dt*nt, nt);
[X, T] = meshgrid(x, t);

%Defining the initial condition for the wave
u = ones(1, nx); 
u(1, 0.5./dx : 1./dx+1) = 2; %Changing all values of u from x = 0.5 to x = 1 to u = 2 instead of 1
l = plot(x, u); axis tight manual; grid on; xlabel('x'); ylabel('u');
u_total = zeros(nt, nx);
u_total(1, :) = u;

%PDE solver - we need to iterate over n (index for t) and i (index for x)
for n=1:1:nt
    
    %This for loop is updating the whole wave function u using the
    %differences formula
    %Start from 2 as we reference the (i-1)th element
    %This whole loop is done nt times to produce the wave function at each
    %time t which we plot each time

    un = u; %Old u
    for i=2:1:nx        
        u(i) = un(i) - un(i) * dt/dx * (un(i) - un(i-1));
    end
    u_total(n, :) = u;
    set(l, 'XData', x, 'YData', u); title('1D Convection, Time: ', t(n))
    pause(0.0001);
end

figure();
a = surface(X, T, u_total); 
set(a, 'LineStyle', 'none');