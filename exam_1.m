% clc; clear; close all;

%set default figure properties
set(0,'DefaultLineLineWidth',2)
set(0,'defaultfigurecolor',[1 1 1])

% import functions
addpath dynamics dmd edmd training prediction

%% Problem setup

dynamics_fn = @lorenz_system; % point to lorenz dynamics

% operator approximation method
% set to EDMD_flag = false for DMD.
% set to EDMD_flag = true for EDMD.
 

%% Load training dataset and visualize
% load the csv file for pendulum training
dt = 0.01;
x0 = [0; 1; 20];
xf = [0;0;0];
u = @(t) 0;
tspan = 0:dt:50;  % Adjusted to start from 0
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[t, data] = ode45(@(t,x) lorenz_system(t, x, u(t)), tspan, x0, options);  % Removed unnecessary input u


%% Plot zero/step control input
hold on
plot3(data(:,1), data(:,2), data(:,3));
plot3(x0(1),x0(2),x0(3), 'ro', 'MarkerSize', 20);
plot3(xf(1),xf(2),xf(3), 'ro', 'MarkerSize', 20);
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz Attractor');
view(3);
grid on;
hold off;
%% %%%%%%%%%%% Visualize the training dataset %%%%%%%%%%%%%%%%%% %%
% complete code to plot phase portraits of training data
n_traj = 100;
x0 = randi(20,3*n_traj,1);
u = @(t) sin(t);
tspan = 0:dt:10;  % Adjusted to start from 0
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
data_table = [];
for i = 1:3:3*n_traj
    [~, data] = ode45(@(t,x) lorenz_system(t, x, u(t)), tspan, x0(i:i+2,1), options);  % Removed unnecessary input u
    data_table = [data';data_table];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Generate time-shfited snapshots %%%%%%%%%%%%%%%%%%%%
% complete code to obtain snapshots X1 and X2 from dataset X

% STEP 0: data split
X1 = data_table(:,1:end-1);
X2 = data_table(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get Koopman operator
EDMD_flag = false;
if(EDMD_flag)
    % get Koopman operator using EDMD
    % define basis setup
    basis.dim = 3; % system dimension

    %%%%%%%%%%%%% Setup basis properties for EDMD %%%%%%%%%%%%
    % complete code
    % specify basis type as a string below ('monomials or rbf')

    basis.type = 'rbf';
    switch basis.type
    case 'monomials'
        % specifiy degree of monomial sbelow
        basis.deg = 5;
    case 'rbf'
        % specifiy kernel width of rbfs
        basis.centers = 50;
    otherwise
        disp('Please specify type of basis')
    end

    %%%%%%%%%% Obtain Koopman operator using EDMD %%%%%%%%%%%%
        
    % complete code for function get_EDMD()
    [A, C, D_sorted] = get_EDMD(X1,X2, basis);   
    operator.A = A;
    operator.C = C;
    D = D_sorted(abs(D_sorted) < 1);
    D = D(D~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    %%%%%%%%%% Obtain Koopman operator using DMD %%%%%%%%%%%%

    % complete code for function get_DMD()
    [A,Phi] = get_DMD(X1,X2);
    operator.A = A;
    operator.Phi = Phi;
    % D_sorted = eig(A);
    % D = D_sorted(abs(D_sorted) < 1);
    % D = D(D~=0);
    operator.D = eig(A);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Plot eigenvalues of Koopman operator
%%%%%%%%% Plot eigenvlaues of A in a unit circle %%%%%%%%%%%%%%%%
% complete code to plot discrete time eigenvalues of A

% Plot the eigenvalues in the complex plane (unit circle)
figure;
hold on;
plot(real(D), imag(D), 'bo'); % Plot eigenvalues as blue circles
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), 'r--'); % Plot unit circle in red dashed line
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues of Koopman Operator (DMD)');
legend('Eigenvalues', 'Unit Circle');
axis equal;
grid on;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% evaluate operator for n timesteps prediction
%% %%%%%%%%%%% Visualize the training dataset %%%%%%%%%%%%%%%%%% %%
% complete code to plot phase portraits of training data
%% evaluate operator for n timesteps prediction
% Use the following parameters for validation

prediction.n_steps = 20; % num timesteps to predict
prediction.dt = 0.1;
prediction.show_plot = true;


n_eval = 20;
x0_eval = randi(20,3*n_eval,1);
u = @(t) sin(t);
tspan = 0:dt:10;  % Adjusted to start from 0
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
X_eval = [];
for i = 1:3:3*n_eval
    [~, data_eval] = ode45(@(t,x) lorenz_system(t, x, u(t)), tspan, x0_eval(i:i+2,1), options);  % Removed unnecessary input u
    X_eval = [data_eval';X_eval];
end

%% Predict data using DMD/EDMD modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_pred = eval_prediction(data_table,operator,prediction);

%% evaluate the avg rmse and % error across for prediction
X_true = X_eval(1:prediction.n_steps,:);
RMSE = rmse(X_pred, X_true, prediction);






