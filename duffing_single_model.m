clc; clear; close all;

%set default figure properties
set(0,'DefaultLineLineWidth',2)
set(0,'defaultfigurecolor',[1 1 1])

% fix random seed
seed = 1;
rng(seed);

% import functions
addpath dynamics dmd edmd training prediction

%% Problem setup
% for duffing, states = [x1, x2].
dynamics_fn = @dynamics_duffing; % point to pendulum dynamics

% operator approximation method
% set to EDMD_flag = false for DMD.
% set to EDMD_flag = true for EDMD.
EDMD_flag = true; 


%% Load training trajectory data using full domain
data = load("training/duffing_training.csv");

%%%%%%%%%%%%% Visualize the training dataset %%%%%%%%%%%%%%%%%%%%
% complete code to plot phase portraits of training data
% Number of trajectories in the dataset
num_trajectories = 100;

% Sampling rate
sampling_rate = 0.1;

% Number of data points in each trajectory
num_points = 10;

% Time span for each trajectory
time_span = 1;

% Visualize a few sample trajectories by plotting a phase portrait (θ vs ẋ)
figure;
hold on;
data_new = [];
index = 1;
theta = data(index:num_points+1,1);
theta_dot = data(index:num_points+1,2);
plot(theta, theta_dot);
data_new = [theta; theta_dot];
for i = 2:num_trajectories % Plot up to 5 trajectories
    hold on;
    index = num_points*(i-1)+i;
    % Extract theta and theta_dot for the current trajectory
    theta = data(index:index+num_points,1);
    theta_dot = data(index:index+num_points,2);
    % Plot phase portrait (theta vs theta_dot)
    plot(theta, theta_dot);
    data_new = [data_new,[theta;theta_dot]];
end

% Set labels and title
xlabel('\theta');
ylabel('\theta-dot');
title('Phase Portrait of Duffing Training Data');

% Add legend
legend('Trajectory 1', 'Trajectory 2', 'Trajectory 3', 'Trajectory 4', 'Trajectory 5');

% Display grid
grid on;

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Generate time-shfited snapshots %%%%%%%%%%%%%%%%%%%%
% complete code to obtain snapshots X1 and X2 from dataset X
% STEP 0: data split
X1 = data_new(:,1:end-1);
X2 = data_new(:,2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get Koopman operator
if(EDMD_flag)
    % get Koopman operator using EDMD
    % define basis setup
    basis.dim = 2; % system dimension

    %%%%%%%%%%%%% Setup basis properties for EDMD %%%%%%%%%%%%
    % complete code
    % specify basis type as a string below ('monomials or rbf'
    basis.type = 'rbf';

    switch basis.type
    case 'monomials'
        % specifiy degree of monomial sbelow
        % basis.deg = ??;
    case 'rbf'
        % specifiy kernel width of rbfs
        basis.gamma = 50;
    otherwise
        disp('Please specify type of basis')
    end

    %%%%%%%%%% Obtain Koopman operator using EDMD %%%%%%%%%%%%
        
    % complete code for function get_EDMD()
    [A, C, D_sorted] = get_EDMD(X1, X2, basis);   
    operator.A = A;
    operator.C = C;
    D = D_sorted(abs(D_sorted) < 1);
    D = D(D~=0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    %%%%%%%%%% Obtain Koopman operator using DMD %%%%%%%%%%%%
    basis = [];

    % complete code for function get_DMD()
    A = get_DMD(X1, X2);
    operator.A = A;

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
% Use the following parameters for validation
prediction.n_traj = 9; %traj for evaluations
prediction.n_steps = 20; % num timesteps to predict
prediction.dt = 0.1;
prediction.show_plot = true;

% load the validation dataset
X_eval = load('prediction\duffing_validation.csv');
% Number of data points in each trajectory
prediction.num_points = size(X_eval, 1)/prediction.n_traj;

% Visualize a few sample trajectories by plotting a phase portrait (θ vs ẋ)
figure;
hold on;

index = 1;
theta_eval = X_eval(index:prediction.num_points,1);
theta_dot_eval = X_eval(index:prediction.num_points,2);
plot(theta_eval, theta_dot_eval);

for i = 2:prediction.n_traj % Plot up to 5 trajectories
    hold on;
    index = prediction.num_points *(i-1)+1;
    % Extract theta and theta_dot for the current trajectory
    theta_eval = X_eval(index:index+prediction.num_points -1,1);
    theta_dot_eval = X_eval(index:index+prediction.num_points -1,2);
    % Plot phase portrait (theta vs theta_dot)
    plot(theta_eval, theta_dot_eval);
end

% Set labels and title
xlabel('\theta_{Eval}');
ylabel('\theta-dot_{Eval}');
title('Phase Portrait of Pendulum Prediction Data');

% Add legend
legend('Trajectory 1', 'Trajectory 2', 'Trajectory 3', 'Trajectory 4', 'Trajectory 5', ...
    'Trajectory 6', 'Trajectory 7', 'Trajectory 8', 'Trajectory 9');

% Display grid
grid on;

hold off;

% Generate evaluation data
data_eval_new = zeros(size(X1, 1), prediction.n_traj  * 10);
data_eval_new(1:11,1) = X_eval(1:11,1); 
data_eval_new(12:22,1) = X_eval(1:11,2); 

%%%%%%%%%%%%% Get predictions using Koopman  %%%%%%%%%%%%%%%%%%%%
% complete code for function eval_prediction()
X_pred = eval_prediction(X_eval,operaotr,basis,prediction);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate the avg rmse and % error across for prediction
X_true = X_eval(1:prediction.n_steps,:);
RMSE = rmse(X_pred, X_true, prediction);
