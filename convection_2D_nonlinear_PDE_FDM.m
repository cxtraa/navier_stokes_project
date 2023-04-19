%Nonlinear convection in 2D
%Code is essentially the same as linear 2D convection, however, we
%introduce the 3D array v as there are two coupled PDEs now.

%Setting up all parameters
nx = 50; xmax = 3; dx = xmax./(nx-1);
ny = 50; ymax = 3; dy = ymax./(nx-1);
nt = 100; sigma = 0.1; dt = sigma*dx;
c = 1;

%Creating matrices
x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
[X, Y] = meshgrid(x, y);
u = ones(ny, nx, nt);
v = ones(ny, nx, nt);

%Setting the initial conditions
u((0.5./dy) : (1./dy), (0.5./dx) : (1./dx), 1) = 2;
v((0.5./dy) : (1./dy), (0.5./dx) : (1./dx), 1) = 5;

%Iterating through time, and both x and y coordinates
for n=1:1:nt
    for j=2:1:ny-1
        for i=2:1:nx-1
            %Iterative formula from the PDE
            u(j, i, n+1) = u(j, i, n) - u(j,i,n)*(dt/dx)*(u(j,i,n)-u(j,i-1,n)) - v(j,i,n)*(dt/dy)*(u(j,i,n)-u(j-1,i,n));
            v(j, i, n+1) = v(j, i, n) - u(j,i,n)*(dt/dx)*(v(j,i,n)-v(j,i-1,n)) - v(j,i,n)*(dt/dy)*(v(j,i,n)-v(j-1,i,n)); %Extra iteration for v

            %Boundary conditions (borders = 1 always)
            u(1, :, n+1) = 1;
            u(:, 1, n+1) = 1;
            u(ny, :, n+1) = 1;
            u(:, nx, n+1) = 1;

            v(1, :, n+1) = 1; %Extra BCs for v(x, y, t)
            v(:, 1, n+1) = 1;
            v(ny, :, n+1) = 1;
            v(:, nx, n+1) = 1;
        end
    end
end

%Velocity field - this code is the convective acceleration term in the
%Navier Stokes equation set to zero (i.e. fluid behaviour if net force = 0)
figure();
a = quiver(x, y, u(:, :, 1), v(:, :, 1)); grid on; axis tight manual;
set(a, 'AutoScaleFactor', 1.5);
xlabel('X/m'); ylabel('Y/m');
title('Velocity field of fluid: nonlinear 2D convection');
for n=2:nt
    set(a, 'udata', u(:,:,n), 'vdata', v(:,:,n));
    pause(0.2);
end


%Displaying the results
figure();
axh = axes;
hold(axh, 'on');
a = surface(axh, X, Y, u(:, :, 1)); grid on; axis tight manual;
xlabel('X/m'); ylabel('Y/m'); zlabel('u/ms^-1');
title('Evolution of velocity distribution of fluid in 2D space over time');
for n=2:nt
    delete(a);
    a = surface(axh, X, Y, u(:, :, n));
    pause(0.1);
end



