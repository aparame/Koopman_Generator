clc; clear; close all;

%set default figure properties
set(0,'DefaultLineLineWidth',2)
set(0,'defaultfigurecolor',[1 1 1])

% import functions
addpath dynamics dmd edmd training prediction

%% Problem setup
% for pendulum, states = [theta, theta_dot].
dynamics_fn = @dynamics_pendulum; % point to pendulum dynamics

% operator approximation method
% set to EDMD_flag = false for DMD.
% set to EDMD_flag = true for EDMD.
EDMD_flag = false; 

%% Load training dataset and visualize
% load the csv file for pendulum training
data = load("training/pendulum_training.csv");


%%%%%%%%%%%%% Visualize the training dataset %%%%%%%%%%%%%%%%%%%%
% complete code to plot phase portraits of training data
% Load the training dataset from 'pendulum_training.csv'

% Number of trajectories in the dataset
num_trajectories = 126;

% Sampling rate
sampling_rate = 0.1;

% Number of data points in each trajectory
num_points = size(data, 1)/num_trajectories;

% Time span for each trajectory
time_span = 1;

% Visualize a few sample trajectories by plotting a phase portrait (θ vs ẋ)
figure;
hold on;
data_new = [];
index = 1;
theta = data(index:num_points,1);
theta_dot = data(index:num_points,2);
plot(theta, theta_dot);
for i = 2:5 % Plot up to 5 trajectories
    hold on;
    index = num_points*(i-1)+1;
    % Extract theta and theta_dot for the current trajectory
    theta = data(index:index+num_points-1,1);
    theta_dot = data(index:index+num_points-1,2);
    % Plot phase portrait (theta vs theta_dot)
    plot(theta, theta_dot);
end

% Set labels and title
xlabel('\theta');
ylabel('\theta-dot');
title('Phase Portrait of Pendulum Training Data');

% Add legend
legend('Trajectory 1', 'Trajectory 2', 'Trajectory 3', 'Trajectory 4', 'Trajectory 5');

% Display grid
grid on;

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% Generate time-shfited snapshots %%%%%%%%%%%%%%%%%%%%
% complete code to obtain snapshots X1 and X2 from dataset X
data_new = [data(1:11,1);data(1:11,2)];
for ii = 2:num_trajectories % Plot up to 5 trajectories
    index_1 = num_points*(ii-1)+1;
    % Extract theta and theta_dot for the current trajectory
    theta = data(index_1:index_1+num_points-1,1);
    theta_dot = data(index_1:index_1+num_points-1,2);
    data_new = [data_new,[theta;theta_dot]];
  
end
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
    % specify basis type as a string below ('monomials or rbf')

    % basis.type = ?? 
    switch basis.type
    case 'monomials'
        % specifiy degree of monomial sbelow
        % basis.deg = ??;
    case 'rbf'
        % specifiy kernel width of rbfs
        % basis.gamma = 0.5;
    otherwise
        disp('Please specify type of basis')
    end

    %%%%%%%%%% Obtain Koopman operator using EDMD %%%%%%%%%%%%
        
    % complete code for function get_EDMD()
    [A, C] = get_EDMD(X1,X2, basis);   
    operator.A = A;
    operator.C = C;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    %%%%%%%%%% Obtain Koopman operator using DMD %%%%%%%%%%%%
    basis = [];

    % complete code for function get_DMD()
    A = get_DMD(X1,X2);
    operator.A = A;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Plot eigenvalues of Koopman operator
%%%%%%%%% Plot eigenvlaues of A in a unit circle %%%%%%%%%%%%%%%%
% complete code to plot discrete time eigenvalues of A


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% evaluate operator for n timesteps prediction
% Use the following parameters for validation
prediction.n_traj = 9; %traj for evaluations
prediction.n_steps = 20; % num timesteps to predict
prediction.dt = 0.1;
prediction.show_plot = true;

% load the validation dataset
X_eval = load('prediction\pendulum_validation.csv');


%%%%%%%%%%%%% Get predictions using Koopman  %%%%%%%%%%%%%%%%%%%%
% complete code for function eval_prediction()
X_pred = eval_prediction(X_eval,operator,basis,prediction);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate the avg rmse and % error across for prediction
X_true = X_eval(1:prediction.n_steps,:);
RMSE = rmse(X_pred, X_true, prediction);
