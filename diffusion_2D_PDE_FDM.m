%Diffusion equation in 2D

%Setting up all parameters
nu = 0.05;
nx = 100; xmax = 2; dx = xmax./(nx-1);
ny = 100; ymax = 2; dy = ymax./(nx-1);
nt = 200; sigma = 0.25; dt = sigma*dx*dy/nu;

%Creating matrices
x = linspace(0, xmax, nx);
y = linspace(0, ymax, ny);
[X, Y] = meshgrid(x, y);
u = ones(ny, nx, nt);

%Setting the initial conditions
u((0.5./dy) : (1./dy), (0.5./dx) : (1./dx), 1) = 2;

%Iterating through time, and both x and y coordinates
for n=1:1:nt
    for j=2:1:ny-1
        for i=2:1:nx-1
            %Iterative formula from the PDE
            u(j,i,n+1) = u(j,i,n) + nu*dt/dx^2*(u(j,i+1,n) - 2*u(j,i,n) + u(j,i-1,n)) + nu*dt/dy^2*(u(j+1,i,n) - 2*u(j,i,n) + u(j-1,i,n));

            %Boundary conditions (borders = 1 always)
            u(1, :, n+1) = 1;
            u(:, 1, n+1) = 1;
            u(ny, :, n+1) = 1;
            u(:, nx, n+1) = 1;
        end
    end
end

%Displaying the results
figure();
axh = axes;
hold(axh, 'on'); %Fix axis limits
a = surface(axh, X, Y, u(:, :, 1)); grid on; axis tight manual;
for n=2:nt
    %Delete old surface and plot new surface with 0.1s interval
    delete(a);
    a = surface(axh, X, Y, u(:, :, n));
    pause(0.1);
end



