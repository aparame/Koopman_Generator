clc
clear

%% Problem 2 Part A
pendulum_A = [0 1;-9.8 -2];
eig_A = eig(pendulum_A);

%% Problem 2 Part B
% Load the training dataset from 'pendulum_training.csv'
data = csvread('training/pendulum_training.csv');

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
for i = 1:5 % Plot up to 5 trajectories
    % Extract theta and theta_dot for the current trajectory
    index = num_points*(i-1)+2*i;
    theta = data(index:num_points,1);
    theta_dot = data(index:num_points,2);
    
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

%% Problem 2 Part 3

