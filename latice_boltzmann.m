%Lattice-Boltzmann method of fluid simulation over a cylinder

clear;

nx = 400;
ny = 100;
[X,Y] = meshgrid(1:nx, 1:ny);
tau = 0.53;
nt = 50;
p_every = 50;

cxs = [0 0 1 1  1  0 -1 -1 -1];
cys = [0 1 1 0 -1 -1 -1 0   1];
ws = [4/9 1/9 1/36 1/9 1/36 1/9 1/36 1/9 1/36];

F = ones(ny, nx, 9) + 0.01*randn(ny, nx, 9);
F(:, :, 4) = 2.3;

cylinder = (X - nx/4).^2 + (Y-ny/2).^2 < (ny/4).^2;

for n = 1:nt

    for i=1:9
        for cx = cxs
            for cy = cys
                F(:, :, i) = circshift(F(:, :, i), cx, 2);
                F(:, :, i) = circshift(F(:, :, i), cy, 1);
            end
        end
    end

    for i=1:9
        F_temp = F(:, :, i);
        bndryF_temp = F_temp(cylinder);
        bndryF(:, i) = bndryF_temp;
    end

    bndryF = bndryF(:, [1 6 7 8 9 2 3 4 5]);

    rho = sum(F, 3);
    vx = sum(F.*reshape(cxs,1,1,[]), 3) ./ rho;
    vy = sum(F.*reshape(cys,1,1,[]), 3) ./ rho;

    for i=1:9
        F_temp = F(:, :, i);
        bndryF_temp = bndryF(:, i);
        F_temp(cylinder) = bndryF_temp;
        F(:, :, i) = F_temp;
    end

    vx(cylinder) = 0;
    vy(cylinder) = 0;

    Feq = zeros(size(F));
    for i=1:9
        for cx = cxs
            for cy = cys
                for w = ws
                    Feq = w.*rho.*(1 + 3*(cx*vx + cy*vy) + 9/2 * (cx*vx + cy*vy).^2 + 3/2 * (vx.^2 + vy.^2).^2);
                end
            end
        end
    end
    
    F = F - (1/tau) * (F-Feq);
end

surf(sqrt(vx.^2 + vy.^2)); shading interp; colormap jet; axis equal;