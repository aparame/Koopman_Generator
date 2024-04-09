clear; close all; clc;

set(0,'DefaultLineLineWidth',2)
set(0,'defaultfigurecolor',[1 1 1])

addpath dynamics edmd eigenfunctions

%% Problem setup

clear; close all; clc;

% Dynamics
system.n = 2;

syms x1 x2 t

f = [((7.5*(x2^2)+5)*((x1^3)+x1+sin(x2))+(-x1+(x2^3)+2*x2)*cos(x2))/(9*(x1^2)*(x2^2)+6*(x1^2)+3*(x2^2)+cos(x2)+2);...
                (2.5*(x1^3)+2.5*x1-(3*(x1^2)+1)*(-x1+(x2^3)+2*x2)+2.5*sin(x2))/(9*(x1^2)*(x2^2)+6*(x1^2)+3*(x2^2)+cos(x2)+2)];

%%%%%%%%%%%%%%%%% Set up true eigenfunctions %%%%%%%%%%%%%%%%%%%%%
% complete code block to obtain true eig vals and eigfns
% true eigval (diagonal matrix)

% Calculate the Jacobian matrix
J = jacobian(f, [x1 x2]);

% Find the equilibrium point(s)
[x1_eq, x2_eq] = solve(f == [0; 0], [x1 x2]);

% Calculate the eigenvalues and eigenvectors at the equilibrium point(s)
A_cell = cell(length(x1_eq), 1);
for i = 1:length(x1_eq)
    J_eq = subs(J, [x1 x2], [x1_eq(i), x2_eq(i)]);
    A_cell{i} = double(subs(J_eq));
    [eig_vectors, lambda] = eig(J_eq);
    fprintf('Eigenvalues at equilibrium point (%f, %f):\n', double(x1_eq(i)), double(x2_eq(i)));
    disp(lambda);
    fprintf('Eigenvectors:\n');
    disp(eig_vectors);
end
lambda_val = diag(lambda);
% Compute the eigenvectors of the system matrix A
eigenfunctions = cell(length(A_cell), 1);
for i = 1:length(A_cell)
    [V, ~] = eig(A_cell{i});
    eigenfunctions{i} = V;
end

% % Display the eigenfunctions as functions of x
% for i = 1:length(eigenfunctions)
%     fprintf('Eigenfunction at equilibrium point (%f, %f):\n', double(x1_eq(i)), double(x2_eq(i)));
%     for j = 1:length(eigenfunctions{i})
%         fprintf('Eigenfunction %d:\n', j);
%         % Define the eigenfunction as a function of x
%         eigenfunc_x = matlabFunction(eigenfunctions{i}(:, j), 'vars', [x1, x2]);
%         disp(eigenfunc_x(x1, x2));
%     end
% end
% true eigenfunctions (inline function)
Phi = [x1 - 2*x2 - (x2^3); x1 + sin(x2) + x1^3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% corresponding system dynamics
% dynamics_fn = dynamics_custom(lambda, Phi);
f_func = matlabFunction(f,'vars',{x1,x2});

%%%%%%%%%%%%% Visualize the phase portrait %%%%%%%%%%%%%%%%%%%%
% complete code to plot phase portrait of the system
% Define range for plotting
x1_range = linspace(-5, 5, 20);
x2_range = linspace(-5, 5, 20);

% Create a grid of points
[x1_grid, x2_grid] = meshgrid(x1_range, x2_range);

% Evaluate the system function at each point in the grid
dx1_grid = zeros(size(x1_grid));
dx2_grid = zeros(size(x2_grid));
for i = 1:numel(x1_grid)
    dx = f_func(x1_grid(i), x2_grid(i));
    dx1_grid(i) = dx(1);
    dx2_grid(i) = dx(2);
end

% Plot the phase portrait using quiver plot
figure;
quiver(x1_grid, x2_grid, dx1_grid, dx2_grid);
xlabel('x_1');
ylabel('x_2');
title('Phase Portrait');
axis tight;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Linearize the system

%%%%%%%%%%%%% obtain linear dynamics around eqb point %%%%%%%%%%%
% complete code to obtain linearization of the system

A0 = J_eq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basis function setup

% define basis setup
    basis.dim = 2; % system dimension

    %%%%%%%%%%%%% Setup basis properties for EDMD %%%%%%%%%%%%
    % complete code
    % specify basis type as a string below ('monomials or rbf')
    % basis.type = ?? 

    switch basis.type
    case 'monomials'
        % specifiy degree of monomial sbelow
        % basis.deg = ??;
        basis.remove_linear = true;
    case 'rbf'
        % specifiy kernel width of rbfs
        % basis.gamma = 0.5;
    otherwise
        disp('Please specify type of basis')
    end

% obtain nonlinear part of basis functions and its jacobian
psi = get_basis(x, basis);
psi_fn = matlabFunction(psi,'Vars',{x});

psi_grad = jacobian(psi,x);
psi_grad_fn = matlabFunction(psi_grad,'Vars',{x});

basis.fn = psi_fn;
basis.fn_grad = psi_grad_fn;
basis.Nbs = length(psi);

%% Generate training data

%%%%%%%%%%%%% sample points in a uniform grid %%%%%%%%%%%%%%%%%%%%
% complete code to generate grid points
% X_samples = ??


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convex optimization of approximate principl eigfns of system
% use the provided convex_opt fun to obtain approximate eigfuns
[phi1_hat, phi2_hat] = convex_opt(f,A0,basis,X_samples);

%% Compare eigenfunction approximation with true eigenfunctions

%%%%%%%%%%%%%%%%%%%%% Plot true eigenfunctions %%%%%%%%%%%%%%%%%%%
% complete code to generate surface plot of the true eigfns


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Plot true eigenfunctions %%%%%%%%%%%%%%%%%%%%
% complete code to generate surface plot of the approximated eigfns


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% find avg approximation error %%%%%%%%%%%%%%%%
% complete code to find the avg pointwise error


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
