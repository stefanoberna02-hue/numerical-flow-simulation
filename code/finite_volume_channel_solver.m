clear;
clc;
close all;

% Physical Parameters
rho = 1;  % Density kg/m^3
cp = 10;  % Specific heat capacity J/(K*kg)
k = 0.12;  % Thermal conductivity W/(m*K)
nu = 1e-2;  % Kinematic viscosity m^2/s
Gamma = k/cp;  % Heat diffusion coefficient kg/(m*s)
L = 10;  % Channel length m
H = 1;  % Channel height m
Pe = 16.5;  % Peclet number
u_mean = Pe*Gamma/(2*H*rho);  % Average speed m/s
Re = 2*H*u_mean/nu;  % Reynolds number

% Boundary Temperature Parameters
T_in = 50;  % Inlet temperature °C
T_wall = 100;  % Wall temperature °C

% Velocity Field Function
u_x = @(y) 6*u_mean.*(y/H).*(1-(y/H));

function [A, b, T_2D] = compute_sol(nx, ny, L, H, Gamma, rho, u_x, T_in, T_wall)
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny; % DOF number

% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Fill in the Equations for the internal CV
% Diffusion coefficient
D_e = Gamma*dy/dx;  % Eastern
D_w = Gamma*dy/dx;  % Western
D_n = Gamma*dx/dy;  % Northern
D_s = Gamma*dx/dy;  % Southern

for i = 2:(ny-1)
    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_E = D_e + max(-F_e,0);
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_E + a_W + a_N + a_S;

    for j = 2:(nx-1)
        n = (i-1)*nx + j;

        % Fill matrix A (The internal CV has no source term, so b(n) is zero
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        A(n,n-1) = -a_W;
        A(n,n-nx) = -a_N;
        A(n,n+nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for i = 2:ny-1
    n = (i-1)*nx + 1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann)
for i = 2:ny-1
    n = (i-1)*nx + nx;

    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_W + a_N + a_S;

    A(n,n) = a_P;
    A(n,n+1) = 0;
    A(n,n-1) = -a_W;
    A(n, n-nx) = -a_N;
    A(n, n+nx) = -a_S;
    b(n) = 0;
end

% Fill boundary conditions at the walls including 4 corner CVs (Dirichlet T=T_wall)
for j = 1:nx
    % Lower wall surface
    n = (ny-1)*nx + j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;

    % Upper wall surface
    n = j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
end

% Solve the linear equation system
T = A\b;

T_2D = reshape(T, nx, ny)';
end

nx = 50;
ny = 5;
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny;

[A, b, T_2D] = compute_sol(nx, ny, L, H, Gamma, rho, u_x, T_in, T_wall);
x_coords = linspace(0, L, nx);
y_coords = linspace(0, H, ny);
[X, Y] = meshgrid(x_coords, y_coords);

% Visualization and Validation of Temperature Fields
figure('Name', 'Coarse-mesh temperature field');
pcolor(X, Y, T_2D);
shading interp;
colorbar;
xlabel('x (m)', 'FontSize', 15);
ylabel('y (m)', 'FontSize', 15);
%title(['Channel temperature field（nx=' num2str(nx) ', ny=' num2str(ny) '）']);
axis equal tight;

%%% 12.
figure('Name', 'Boundary conditions validation');

% Verify inlet temperature (plotting the temperature curve at x=0)
subplot(2, 2, 1);
plot(y_coords, T_2D(:, 1), 'b-o', 'MarkerFaceColor', 'blue');
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Inlet');
grid on;

% Verify visualizartion of the outlet temperature
subplot(2, 2, 2);
plot(x_coords, T_2D(2,:), 'b');
hold on;
plot(x_coords, T_2D(ny-2,:), 'b');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Outlet');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=0)
subplot(2, 2, 3);
plot(x_coords, T_2D(ny, :), 'r-o');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = 0)');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=H)
subplot(2, 2, 4);
plot(x_coords, T_2D(1,:), 'r-o');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = H)');
grid on;

%%% 13.
figure();
plot(y_coords, T_2D(:, nx), 'b-o', 'MarkerFaceColor', 'blue', 'LineWidth', 2);
xlabel('y (m)',  'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Outlet');
grid on;

figure();
plot(x_coords, T_2D(ceil(ny/2),:), 'r-o', 'LineWidth', 2);
hold on;

T_mean_weighted_values = zeros(nx, 1);
int_u_x = zeros(nx, 1);

for j=1:nx
    for i=1:ny
        if i==1 || i==ny
            T_mean_weighted_values(j) = T_mean_weighted_values(j)+ u_x(y_coords(i))*(T_2D(i,j))*dy/2;
            int_u_x(j)=int_u_x(j)+u_x(y_coords(i))*dy/2;
        else
            T_mean_weighted_values(j) = T_mean_weighted_values(j)+ u_x(y_coords(i))*(T_2D(i,j))*dy;
            int_u_x(j)=int_u_x(j)+u_x(y_coords(i))*dy;
        end
    end
    T_mean_weighted_values(j)=T_mean_weighted_values(j)/int_u_x(j);
end

plot(x_coords, T_mean_weighted_values, 'b-o', 'LineWidth', 2);
xlabel('x (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Wall (y = H)');
grid on;

x_entry_index = 0;
for i=1:nx
    if T_2D(ceil(ny/2),i) >= 0.9 * T_wall
        x_entry_index = i;
        break
    end
end

xline(x_coords(x_entry_index), '-', 'LineWidth', 3);
legend('Centerline temperature', 'Weighted mean temperature', 'Entry length', 'FontSize', 20)
disp(x_coords(x_entry_index));

%%% 14.
function [error, residual, it, sol] = sor(A, b, w, guess, tol)
n = size(b,1);
sol = guess;
solnew = zeros(n,1);

error = [];
residual = [];

it = 0;
var1 = tol + 1;
var2 = tol + 1;
while var1 > tol || var2 > tol
    for i=1:n
        solnew(i) = (b(i) - A(i, i+1:end)*sol(i+1:end) - A(i, 1:i-1)*solnew(1:i-1)) / A(i, i);
    end

    solold = sol;
    sol = (1 - w) * sol + w * solnew;
    var1 = sqrt(sum((sol - solold).^2)) / sqrt(sum((solold).^2));
    var2 = sqrt(sum((A*sol - b).^2)) / sqrt(sum((diag(A)' * sol).^2));

    error(end + 1) = var1;
    residual(end + 1) = var2;
    it = it + 1;
end
end

guess = 50*ones(N, 1);

for j = 1:nx
    % Lower wall surface
    n = (ny-1)*nx + j;
    guess(n) = T_wall;

    % Upper wall surface
    n = j;
    guess(n) = T_wall;
end

[error1, residual1, it1, ~] = sor(A, b, 1, guess, 1e-5);
[error2, residual2, it2, sol] = sor(A, b, 1.5, guess, 1e-5);

figure();
iter1 = 1:1:it1;
iter2 = 1:1:it2;
plot(iter1, error1, 'LineWidth', 2);
hold on;
plot(iter2, error2, 'LineWidth', 2);
grid on;
xlabel('Iterations', 'FontSize', 20);
ylabel('Relative iteration error', 'FontSize', 20);
legend('Relative error (w=1)', 'Relative error (w=1.5)', 'FontSize', 20);

figure();
plot(iter1, residual1, 'LineWidth', 2);
hold on;
plot(iter2, residual2, 'LineWidth', 2);
grid on;
xlabel('Iterations', 'FontSize', 20);
ylabel('Normalized residuals', 'FontSize', 20);
legend('Normalized residual (w=1)', 'Normalized residual (w=1.5)', 'FontSize', 20);

disp(it1);
disp(it2);

%%% 15.
q = -((T_2D(ny-1, :) - T_2D(ny, :)) / dy)';
deltaT = T_wall - T_mean_weighted_values;
Nu = q ./ deltaT * 2 * H;

figure();
plot(x_coords, Nu, 'LineWidth', 2);
grid on;
xlabel('x (m)', 'FontSize', 20);
ylabel('Nu', 'FontSize', 20);

%%% 16.
nx_vec = [50, 100, 200, 400];
ny_vec = [5, 11, 21, 41];
x_coords_vec = zeros(4, 400);
y_coords_vec = zeros(4, 41);
outlet_temperatures = zeros(4, 41);
centerline_temperatures = zeros(4, 400);
vam_temperatures = zeros(4, 400);
nus = zeros(4, 400);

for i=1:4
    dx = L/(nx_vec(i)-1);
    dy = H/(ny_vec(i)-1);

    [A, b, T_2D] = compute_sol(nx_vec(i), ny_vec(i), L, H, Gamma, rho, u_x, T_in, T_wall);
    x_coords_vec(i, 1:nx_vec(i)) = linspace(0, L, nx_vec(i));
    y_coords_vec(i, 1:ny_vec(i)) = linspace(0, H, ny_vec(i));
    outlet_temperatures(i, 1:ny_vec(i)) = T_2D(:, nx_vec(i));
    centerline_temperatures(i, 1:nx_vec(i)) = T_2D(ceil(ny_vec(i)/2), :);

    int_u_x = zeros(nx_vec(i), 1);

    for j=1:nx_vec(i)
        for z=1:ny_vec(i)
            if z==1 || z==ny_vec(i)
                vam_temperatures(i, j) = vam_temperatures(i, j)+ u_x(y_coords_vec(i, z))*(T_2D(z,j))*dy/2;
                int_u_x(j)=int_u_x(j)+u_x(y_coords_vec(i, z))*dy/2;
            else
                vam_temperatures(i, j) = vam_temperatures(i, j) + u_x(y_coords_vec(i, z))*(T_2D(z,j))*dy;
                int_u_x(j)=int_u_x(j)+u_x(y_coords_vec(i, z))*dy;
            end
        end

        vam_temperatures(i, j) = vam_temperatures(i, j) / int_u_x(j);
    end

    q = -((T_2D(ny_vec(i)-1, :) - T_2D(ny_vec(i), :)) / dy);
    deltaT = T_wall - vam_temperatures(i,1:nx_vec(i));
    nus(i, 1:nx_vec(i)) = q ./ deltaT * 2 * H;

    x_entry_index = 0;
    for z=1:nx_vec(i)
        if T_2D(ceil(ny_vec(i)/2),z) >= 0.9 * T_wall
            x_entry_index = z;
            break
        end
    end

    disp('Entry length')
    disp(x_coords_vec(i, x_entry_index))
end

figure();
for i=1:4
    plot(y_coords_vec(i, 1:ny_vec(i)), outlet_temperatures(i,1:ny_vec(i)), 'LineWidth', 1);
    hold on;
end

grid on;
xlabel('y (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Outlet');
legend('Refinement level 1', 'Refinement level 2', 'Refinement level 3', 'Refinement level 4', 'FontSize', 20)

figure();
for i=1:4
    plot(x_coords_vec(i, 1:nx_vec(i)), centerline_temperatures(i,1:nx_vec(i)), 'LineWidth', 1);
    hold on;
end

grid on;
xlabel('x (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Centerline');
legend('Refinement level 1', 'Refinement level 2', 'Refinement level 3', 'Refinement level 4', 'FontSize', 20)

figure();
for i=1:4
    plot(x_coords_vec(i, 1:nx_vec(i)), vam_temperatures(i,1:nx_vec(i)), 'LineWidth', 1);
    hold on;
end

grid on;
xlabel('x (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Velocity-averaged');
legend('Refinement level 1', 'Refinement level 2', 'Refinement level 3', 'Refinement level 4', 'FontSize', 20)

figure();
for i=1:4
    plot(x_coords_vec(i, 1:nx_vec(i)), nus(i,1:nx_vec(i)), 'LineWidth', 1);
    hold on;
end

grid on;
xlabel('x (m)', 'FontSize', 20);
ylabel('Temperature (°C)', 'FontSize', 20);
%title('Nu');
ylim([0, 50])
legend('Refinement level 1', 'Refinement level 2', 'Refinement level 3', 'Refinement level 4', 'FontSize', 20)

ax_inset = axes('Position', [0.545, 0.3, 0.35, 0.25]);
for i=1:4
    plot(x_coords_vec(i, nx_vec(i)-2*2^i:nx_vec(i)), nus(i,nx_vec(i)-2*2^i:nx_vec(i)), 'LineWidth', 1);
    hold on;
end

grid on;
xlim([9.5, 10])
title(ax_inset, 'Zoomed-in view');
box on;

%%% 17.
function [A, b, T_2D] = compute_sol_quick(nx, ny, L, H, Gamma, rho, u_x, T_in, T_wall)
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny; % DOF number

% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Fill in the Equations for the internal CV
% Diffusion coefficient
D_e = Gamma*dy/dx;  % Eastern
D_w = Gamma*dy/dx;  % Western
D_n = Gamma*dx/dy;  % Northern
D_s = Gamma*dx/dy;  % Southern

for i = 2:(ny-1)
    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations (ud)
    a_E = D_e + max(-F_e,0);
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_E + a_W + a_N + a_S;

    n = (i-1)*nx + 2;

    A(n,n) = a_P;
    A(n,n+1) = -a_E;
    A(n,n-1) = -a_W;
    A(n,n-nx) = -a_N;
    A(n,n+nx) = -a_S;

    % Coefficients of Algebraic Equations (quick)
    a_EE = -max(0, -F_e/8);
    a_E = D_e - max(0,3*F_e/8) + max(0, -6*F_e/8) + max(0, -F_w/8);
    a_WW = -max(0, F_w/8);
    a_W = D_w + max(0,6*F_w/8) + max(0, F_e/8) - max(0, -3*F_w/8);
    a_N = D_n;
    a_S = D_s;
    a_P = a_E + a_EE + a_W + a_WW + a_N + a_S;

    for j = 3:(nx-1)
        n = (i-1)*nx + j;

        % Fill matrix A (The internal CV has no source term, so b(n) is zero
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        %A(n,n+2) = -a_EE;
        A(n,n-1) = -a_W;
        A(n,n-2) = -a_WW;
        A(n,n-nx) = -a_N;
        A(n,n+nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for i = 2:ny-1
    n = (i-1)*nx + 1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann)
for i = 2:ny-1
    n = (i-1)*nx + nx;

    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations (quick)
    a_WW = -max(0, F_w/8);
    a_W = D_w + max(0,6*F_w/8) + max(0, F_e/8) - max(0, -3*F_w/8);
    a_N = D_n;
    a_S = D_s;
    a_P = a_W + a_WW + a_N + a_S;

    A(n,n) = a_P;
    A(n,n-1) = -a_W;
    A(n,n-2) = -a_WW;
    A(n,n-nx) = -a_N;
    A(n,n+nx) = -a_S;
    b(n) = 0;
end

% Fill boundary conditions at the walls including 4 corner CVs (Dirichlet T=T_wall)
for j = 1:nx
    % Lower wall surface
    n = (ny-1)*nx + j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;

    % Upper wall surface
    n = j;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_wall;
end

% Solve the linear equation system
T = A\b;

T_2D = reshape(T, nx, ny)';
end

for i=1:4
    dx = L/(nx_vec(i)-1);
    dy = H/(ny_vec(i)-1);

    [A, b, T_2D] = compute_sol_quick(nx_vec(i), ny_vec(i), L, H, Gamma, rho, u_x, T_in, T_wall);
    x_coords_vec(i, 1:nx_vec(i)) = linspace(0, L, nx_vec(i));

    x_entry_index = 0;
    for z=1:nx_vec(i)
        if T_2D(ceil(ny_vec(i)/2),z) >= 0.9 * T_wall
            x_entry_index = z;
            break
        end
    end

    disp('Entry length')
    disp(x_coords_vec(i, x_entry_index))
end



figure();
plot(2.449, 1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none');
hold on;
plot(2.3232, 1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
plot(2.2111, 1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
plot(2.1805, 1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
plot(2.449, 2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none');
plot(2.2222, 2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
plot(2.1608, 2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
plot(2.1554, 2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
grid on
ylim([0.5, 2.5]);
ax = gca;
ax.YAxis.Visible = 'off';
xlim([2.1, 2.5]);
text(2.15, 0.75, 'Upwind', 'FontSize', 25)
text(2.15, 1.75, 'Quick', 'FontSize', 25)

figure();
base1 = 2.1805;
base2 = 2.1554;
plot(1, (2.449 - base1)/base1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none');
hold on;
plot(1, (2.3232 - base1)/base1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
plot(1, (2.2111 - base1)/base1, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
plot(1, 0, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
plot(2, (2.449 - base2)/base2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'none');
plot(2, (2.2222 - base2)/base2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
plot(2, (2.1608 - base2)/base2, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
plot(2, 0, '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
grid on
ylim([0, 0.25]);
ax = gca;
ax.XAxis.Visible = 'off';
xlim([0.75, 2.25]);
text(0.9, 0.2, 'Upwind', 'FontSize', 25)
text(1.9, 0.2, 'Quick', 'FontSize', 25)

%%% 18.

nx_vec = [400, 400, 400, 100, 200, 400];
ny_vec = [5, 11, 21, 41, 41, 41];
x_coords_vec = zeros(6, 400);
y_coords_vec = zeros(6, 41);
vam_temperatures = zeros(6, 400);
dxs = zeros(6, 1);
dys = zeros(6, 1);
nu_min = zeros(6, 1);

for i=1:6
    dx = L/(nx_vec(i)-1);
    dy = H/(ny_vec(i)-1);
    dxs(i) = dx;
    dys(i) = dy;

    [A, b, T_2D] = compute_sol(nx_vec(i), ny_vec(i), L, H, Gamma, rho, u_x, T_in, T_wall);
    x_coords_vec(i, 1:nx_vec(i)) = linspace(0, L, nx_vec(i));
    y_coords_vec(i, 1:ny_vec(i)) = linspace(0, H, ny_vec(i));

    int_u_x = zeros(nx_vec(i), 1);

    for j=1:nx_vec(i)
        for z=1:ny_vec(i)
            if z==1 || z==ny_vec(i)
                vam_temperatures(i, j) = vam_temperatures(i, j)+ u_x(y_coords_vec(i, z))*(T_2D(z,j))*dy/2;
                int_u_x(j)=int_u_x(j)+u_x(y_coords_vec(i, z))*dy/2;
            else
                vam_temperatures(i, j) = vam_temperatures(i, j) + u_x(y_coords_vec(i, z))*(T_2D(z,j))*dy;
                int_u_x(j)=int_u_x(j)+u_x(y_coords_vec(i, z))*dy;
            end
        end

        vam_temperatures(i, j) = vam_temperatures(i, j) / int_u_x(j);
    end

    q = -((T_2D(ny_vec(i)-1, :) - T_2D(ny_vec(i), :)) / dy);
    deltaT = T_wall - vam_temperatures(i,1:nx_vec(i));
    nu_min(i) = min(q ./ deltaT * 2 * H);
end



figure();
plot(dys(1:3), nu_min(1:3), 'LineWidth', 1);
hold on
plot(dxs(4:6), nu_min(4:6), 'LineWidth', 1)
grid on;
legend('Fixed dx, varying dy', 'Fixed dy, varying dx', 'FontSize', 20)
xlabel('dx / dy (m)', 'FontSize', 20)
ylabel('minimum Nu', 'FontSize', 20)

%%% 19.
% Just change Péclet number at the very top of the script.

%%% 20.+
function [A, b, T_2D] = compute_sol_neumann(nx, ny, L, H, Gamma, rho, cp, u_x, qwall, T_in)
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny; % DOF number

% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Fill in the Equations for the internal CV
% Diffusion coefficient
D_e = Gamma*dy/dx;  % Eastern
D_w = Gamma*dy/dx;  % Western
D_n = Gamma*dx/dy;  % Northern
D_s = Gamma*dx/dy;  % Southern

for i = 2:(ny-1)
    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_E = D_e + max(-F_e,0);
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_E + a_W + a_N + a_S;

    for j = 2:(nx-1)
        n = (i-1)*nx + j;

        % Fill matrix A (The internal CV has no source term, so b(n) is zero
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        A(n,n-1) = -a_W;
        A(n,n-nx) = -a_N;
        A(n,n+nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for i = 1:ny
    n = (i-1)*nx + 1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann)
for i = 1:ny
    n = (i-1)*nx + nx;

    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_W;

    if i ~= ny
        a_P = a_P + a_S;
        A(n, n+nx) = -a_S;
    end

    if i ~= 1
        a_P = a_P + a_N;
        A(n, n-nx) = -a_N;
    end

    if i == 1 || i == ny
        b(n) = dx*qwall/cp;
    end

    A(n,n) = a_P;
    A(n,n-1) = -a_W;
end

% Fill boundary conditions at the walls including 4 corner CV (Neumann)
% LOWER WALL
% Advection coefficient
y_e = 0;  % y_eastern
ux_e = u_x(y_e);  % Mean velocity
F_e = rho * ux_e * dy;  % Eastbound convection flux
F_w = F_e;  % Incompressible flow

% Coefficients of Algebraic Equations
a_E = D_e + max(-F_e,0);
a_W = D_w + max(F_w,0);
a_N = D_n;
a_P = a_E + a_W + a_N;

for j = 2:nx-1
    n = (ny-1)*nx + j;
    A(n,n) = a_P;
    A(n,n+1) = -a_E;
    A(n,n-1) = -a_W;
    A(n,n-nx) = -a_N;
    b(n) = dx*qwall/cp;
end

% UPPER WALL
% Advection coefficient
y_e = H;  % y_eastern
ux_e = u_x(y_e);  % Mean velocity
F_e = rho * ux_e * dy;  % Eastbound convection flux
F_w = F_e;  % Incompressible flow

% Coefficients of Algebraic Equations
a_E = D_e + max(-F_e,0);
a_W = D_w + max(F_w,0);
a_S = D_s;
a_P = a_E + a_W + a_S;

for j = 2:nx-1
    n = j;
    A(n,n) = a_P;
    A(n,n+1) = -a_E;
    A(n,n-1) = -a_W;
    A(n,n+nx) = -a_S;
    b(n) = dx*qwall/cp;
end

% Solve the linear equation system
T = A\b;
T_2D = reshape(T, nx, ny)';
end

nx = 50;
ny = 5;
dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny;

[A, b, T_2D] = compute_sol_neumann(nx, ny, L, H, Gamma, rho, cp, u_x, 10, T_in);
x_coords = linspace(0, L, nx);
y_coords = linspace(0, H, ny);
[X, Y] = meshgrid(x_coords, y_coords);

% Visualization and Validation of Temperature Fields
figure('Name', 'Coarse-mesh temperature field');
pcolor(X, Y, T_2D);
shading interp;
colorbar;
xlabel('x (m)', 'FontSize', 15);
ylabel('y (m)', 'FontSize', 15);
%title(['Channel temperature field（nx=' num2str(nx) ', ny=' num2str(ny) '）']);
axis equal tight;

figure('Name', 'Boundary conditions validation 2');
% Verify inlet temperature (plotting the temperature curve at x=0)
subplot(2, 2, 1);
plot(y_coords, T_2D(:, 1), 'b-o', 'MarkerFaceColor', 'blue');
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Inlet');
grid on;

% Verify visualizartion of the outlet temperature
subplot(2, 2, 2);
plot(y_coords, T_2D(:,nx), 'b');
hold on;
plot(y_coords, T_2D(:,nx-1), 'b');
xlabel('y (m)');
ylabel('Temperature (°C)');
title('Outlet');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=0)
subplot(2, 2, 3);
plot(x_coords, T_2D(ny, :), 'r-o');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = 0)');
grid on;

% Verify wall surface temperature (plotting the temperature curve at y=H)
subplot(2, 2, 4);
plot(x_coords, T_2D(1,:), 'r-o');
xlabel('x (m)');
ylabel('Temperature (°C)');
title('Wall (y = H)');
grid on;

T_mean_weighted_values = zeros(nx, 1);
int_u_x = zeros(nx, 1);

for j=1:nx
    for i=1:ny
        if i==1 || i==ny
            T_mean_weighted_values(j) = T_mean_weighted_values(j)+ u_x(y_coords(i))*(T_2D(i,j))*dy/2;
            int_u_x(j)=int_u_x(j)+u_x(y_coords(i))*dy/2;
        else
            T_mean_weighted_values(j) = T_mean_weighted_values(j)+ u_x(y_coords(i))*(T_2D(i,j))*dy;
            int_u_x(j)=int_u_x(j)+u_x(y_coords(i))*dy;
        end
    end
    T_mean_weighted_values(j)=T_mean_weighted_values(j)/int_u_x(j);
end

q = -((T_2D(ny-1, :) - T_2D(ny, :)) / dy)';
deltaT = T_2D(ny, :)' - T_mean_weighted_values;
Nu = q ./ deltaT * 2 * H;

figure();
plot(x_coords, Nu, 'LineWidth', 2);
grid on;
xlabel('x (m)', 'FontSize', 20);
ylabel('Nu', 'FontSize', 20);
ylim([6.8, 7.6])

%%% 21.
qwall = 10;
fprintf('\n--- Global Heat Balance (Watts) ---\n');

% 1. WALL FLUXES (Neumann) - These are specified, so they are correct.
% Positive when coming in
wall_area = L;
Q_upper_wall = -Gamma*dx*(sum(T_2D(2,2:nx-1) - T_2D(1, 2:nx-1)))/dy;
Q_upper_wall = Q_upper_wall - Gamma*(dx/2)*(T_2D(2,1) - T_2D(1, 1))/dy;
Q_upper_wall = Q_upper_wall - Gamma*(dx/2)*(T_2D(2,end) - T_2D(1, end))/dy;

Q_lower_wall = -Gamma*dx*(sum(T_2D(ny-1,2:nx-1) - T_2D(ny, 2:nx-1)))/dy;
Q_lower_wall = Q_lower_wall - Gamma*(dx/2)*(T_2D(ny-1,1) - T_2D(ny, 1))/dy;
Q_lower_wall = Q_lower_wall - Gamma*(dx/2)*(T_2D(ny-1,end) - T_2D(ny, end))/dy;

fprintf('  Heat IN from Upper Wall: %.4f W\n', Q_upper_wall);
fprintf('  Heat IN from Lower Wall: %.4f W\n', Q_lower_wall);

% 2. INLET FLUXES (Dirichlet)
Q_inlet_conv = 0;
Q_inlet_diff = 0;
for i = 1:ny
    y_e = H - (i-1)*dy;
    ux_e = u_x(y_e);

    % Determine face height: half-height for corners, full-height for interior
    if i == 1 || i == ny
        face_height = dy/2;
    else
        face_height = dy;
    end

    % Advective flux: m_dot * cp * T
    Q_inlet_conv = Q_inlet_conv + (ux_e * face_height) * rho * T_2D(i,2);%positive when entering

    % Diffusive flux: -k * A * (dT/dx)
    Q_inlet_diff = Q_inlet_diff - Gamma * face_height * (T_2D(i,2) - T_2D(i,1)) / dx;%poiìsiutuve when entering
end

fprintf('  Heat IN from Inlet (Conv): %.4f W\n', Q_inlet_conv);
fprintf('  Heat IN from Inlet (Diff): %.4f W\n', Q_inlet_diff);

% 3. OUTLET FLUXES (Neumann)
Q_outlet_conv = 0;
Q_outlet_diff = 0;
for i = 1:ny
    y_e = H - (i-1)*dy;
    ux_e = u_x(y_e);

    % Determine face height: half-height for corners, full-height for interior
    if i == 1 || i == ny
        face_height = dy/2;
    else
        face_height = dy;
    end

    % Advective flux: m_dot * cp * T
    Q_outlet_conv = Q_outlet_conv - (ux_e * face_height * rho)  * T_2D(i,nx-1);

    % Diffusive flux: -k * A * (dT/dx)
    Q_outlet_diff = Q_outlet_diff - Gamma * face_height * (T_2D(i,nx-1) - T_2D(i,nx-2)) / dx;
end

fprintf('  Heat OUT from Outlet (Conv): %.4f W\n', Q_outlet_conv);
fprintf('  Heat OUT from Outlet (Diff): %.4f W\n', Q_outlet_diff);

% 4. THE FINAL TALLY
total_flux = Q_upper_wall + Q_lower_wall + Q_inlet_conv + Q_inlet_diff + Q_outlet_conv + Q_outlet_diff;

fprintf('----------------------------------------\n');
fprintf('  TOTAL NET FLUX: %.4f W\n', total_flux);
fprintf('  (This value should be very close to 0.0)\n');


%%%%%%%%%%



%22. Non uniform wall flow

function [A, b, T_2D,T_outlet_mixedNeu] = compute_sol_mixed_neumann(nx, ny, L, H, Gamma, rho, cp, u_x, qwall, T_in, verbose)

dx = L/(nx-1);
dy = H/(ny-1);
N = nx*ny; % DOF number
x_coords = linspace(0, L, nx);
% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);

% Fill in the Equations for the internal CV
% Diffusion coefficient
D_e = Gamma*dy/dx;  % Eastern
D_w = Gamma*dy/dx;  % Western
D_n = Gamma*dx/dy;  % Northern
D_s = Gamma*dx/dy;  % Southern


% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N, N);
b = zeros(N, 1);


%normal filling of the matrix
for i = 2:(ny-1)
    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_E = D_e + max(-F_e,0);
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_E + a_W + a_N + a_S;

    for j = 2:(nx-1)
        n = (i-1)*nx + j;

        % Fill matrix A (The internal CV has no source term, so b(n) is zero
        A(n,n) = a_P;
        A(n,n+1) = -a_E;
        A(n,n-1) = -a_W;
        A(n,n-nx) = -a_N;
        A(n,n+nx) = -a_S;
    end
end

% Fill boundary conditions at the inlet (Dirichlet T=T_in)
for i = 1:ny
    n = (i-1)*nx + 1;
    A(n,:) = 0;
    A(n,n) = 1;
    b(n) = T_in;
end

% Fill the boundary conditions at the outlet (Neumann)
for i = 1:ny
    n = (i-1)*nx + nx;

    % Advection coefficient
    y_e = H - (i-1)*dy;  % y_eastern
    ux_e = u_x(y_e);  % Mean velocity
    F_e = rho * ux_e * dy;  % Eastbound convection flux
    F_w = F_e;  % Incompressible flow

    % Coefficients of Algebraic Equations
    a_W = D_w + max(F_w,0);
    a_N = D_n;
    a_S = D_s;
    a_P = a_W;

    if i ~= ny
        a_P = a_P + a_S;
        A(n, n+nx) = -a_S;
    end

    if i ~= 1
        a_P = a_P + a_N;
        A(n, n-nx) = -a_N;
    end

    if i == 1 || i == ny %impose the qwall flux on the south and north face, now it should be 0
        b(n) = 0;
    end

    A(n,n) = a_P;
    A(n,n-1) = -a_W;

end

% Fill boundary conditions at the walls including 4 corner CV (Neumann)



% LOWER WALL
% Advection coefficient
y_e = 0;  % y_eastern
ux_e = u_x(y_e);  % Mean velocity
F_e = rho * ux_e * dy;  % Eastbound convection flux
F_w = F_e;  % Incompressible flow

% Coefficients of Algebraic Equations
a_E = D_e + max(-F_e,0);
a_W = D_w + max(F_w,0);
a_N = D_n;
a_P = a_E + a_W + a_N;

for j = 2:nx-1
    n = (ny-1)*nx + j;
    A(n,n) = a_P;
    A(n,n+1) = -a_E;
    A(n,n-1) = -a_W;
    A(n,n-nx) = -a_N;
    if 2<=x_coords(j) && x_coords(j)<=5 
        b(n) = dx*qwall/cp;
    else
        b(n)=0;%otherwise b stays zero as per mixed conditions
    end
    
end

% UPPER WALL
% Advection coefficient
y_e = H;  % y_eastern
ux_e = u_x(y_e);  % Mean velocity
F_e = rho * ux_e * dy;  % Eastbound convection flux
F_w = F_e;  % Incompressible flow

% Coefficients of Algebraic Equations
a_E = D_e + max(-F_e,0);
a_W = D_w + max(F_w,0);
a_S = D_s;
a_P = a_E + a_W + a_S;

for j = 2:nx-1
    n = j;
    A(n,n) = a_P;
    A(n,n+1) = -a_E;
    A(n,n-1) = -a_W;
    A(n,n+nx) = -a_S;

    if 2<=x_coords(j) && x_coords(j)<=5 
        b(n) = dx*qwall/cp;
    else
        b(n)=0;%otherwise b stays zero as per mixed conditions
    end
end

% Solve the linear equation system
T = A\b;
T_2D = reshape(T, nx, ny)';
T_outlet_mixedNeu=mean(T_2D(:,nx));

if verbose
    fprintf("The average outlet temperature when using mixed Neumann with q=%f is T_Out_avg=%f\n", qwall, T_outlet_mixedNeu);
end

end


tested_q=0:1:10;
T_out_avg=zeros(length(tested_q),1);
i=1;
for q_w = tested_q
    [~ , ~, ~,T_outlet_mixedNeu] = compute_sol_mixed_neumann(nx, ny, L, H, Gamma, rho, cp, u_x, q_w, T_in, false);
    T_out_avg(i)=T_outlet_mixedNeu;
    i=i+1;
end


%%%PLOTS
plot(tested_q,T_out_avg')
xlabel('qwall (W/m^2)');
ylabel('Mean Outlet temperature (°C)');
title('Study of average Outlet Temperature in mixed Neumann conditions for different wall fluxes');
hold on
grid on;
%drowing the angular coefficient 
angular_coeffs= (T_out_avg(2:end)-T_out_avg(1:end-1)) ./ (tested_q(2:end)'-tested_q(1:end-1)');%column vector
m=mean(angular_coeffs');
L_tri = 0.5; % horizontal leg length
H_tri = m * L_tri; % vertical leg from slope
x0=mean(tested_q);
y0=mean(T_out_avg');
% Triangle coordinates
Xtri = [x0, x0 + L_tri, x0 + L_tri];
Ytri = [y0, y0, y0 + H_tri];
% Draw triangle
fill(Xtri, Ytri, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'k');
% Add label
text(x0 + L_tri*1.5, y0 , sprintf('m = %.2f', m), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
hold off

%%%Finding the required values through interpolation and check by rerunning
target_Ts=[60, 70 ,80];
q_interpolated=x0 + (target_Ts-y0)/m;

fprintf("\n--- Question 22 ---\n")
for i=1:length(q_interpolated)
    q_w = q_interpolated(i);
    [~ , ~, ~,T_outlet_mixedNeu] = compute_sol_mixed_neumann(nx, ny, L, H, Gamma, rho, cp, u_x, q_w, T_in, false);
    fprintf("qwall value = %f\nExpected temperature is: %f\tComputed temperature is: %f\n\n",q_w,target_Ts(i),T_outlet_mixedNeu)
end