%Velocity potential Laplace equation PDE solution
%Can be adjusted for other PDEs

clear all;

%Mesh variables
nx = 50; %No of x. grid points
ny = 50; %No of y. grid points
lx = 5;
ly = 5;

%Mesh grid
x = linspace(0, lx, nx);
y = linspace(0, ly, ny);
dx = x(2) - x(1); %The points are evenly spaced
dy = y(2) - y(1); 

%Matrix: M * phi = B. This contains all the difference equations. See notes
%taken on 'Solving Potential Flow'
N = nx * ny; %Number of unknowns
M = zeros(N, N);
B = zeros(N, 1);

%Interior mesh points
for i=2:nx-1 %Iterate over x direction, skip 1st and last x points
    for j=2:ny-1
        n = i + (j-1)*nx; %Turn i,j into nth item in phi matrix
        M(n, n) = -4; %All interior diagonals = -4
        M(n, n+1) = 1; %Neighbours = 1
        M(n, n-1) = 1;
        M(n, n+nx) = 1;
        M(n, n-nx) = 1;
        B(n, 1) = 0; %Source term (0 for Laplace equation)
    end
end

%Boundary conditions
%Left BC: phi = y
i = 1;
for j=1:ny
    n = i + (j-1)*nx;
    M(n, n) = 1;
    B(n, 1) = y(j);
end

%Right BC: phi = ly - y
i = nx;
for j=1:ny
    n = i + (j-1)*nx;
    M(n, n) = 1;
    B(n, 1) = ly - y(j);
end

%Top BC: phi = lx - x
j = ny;
for i=1:nx
    n = i + (j-1)*nx;
    M(n, n) = 1;
    B(n, 1) = lx - x(i);
end

%Bottom BC: phi = x
j = 1;
for i=1:nx
    n = i + (j-1)*nx;
    M(n, n) = 1;
    B(n, 1) = x(i);
end

phi_1D = M\B;

%Convert the list of phis back into the 2D mesh grid
for i = 1:nx
    for j = 1:ny
        n = i + (j-1)*nx;
        phi_2D(i, j) = phi_1D(n);
    end
end

%Compute the velocity field
%u = - dphi/dx and v = -dphi/dy
u = zeros(nx, ny); %Create empty matrices to store u and v components
v = zeros(nx, ny);
for i=2:nx-1 %Skip first and last points as these have missing neighbours
    for j=2:ny-1
        %Central difference definition of derivatives
        u(i, j) = -(phi_2D(i+1, j) - phi_2D(i-1, j))/(2*dx);
        v(i, j) = -(phi_2D(i, j+1) - phi_2D(i, j-1))/(2*dy);
    end
end

det(M);

figure();
surface(x, y, phi_2D'); xlabel('x'); ylabel('y'); zlabel('phi'); grid on;
title('Velocity potential in 2D space');

figure();
quiver(x, y, u', v'); xlabel('x'); ylabel('y'); grid on;
title('Velocity field in 2D space');
