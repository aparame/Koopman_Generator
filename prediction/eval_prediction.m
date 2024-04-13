function X_pred = eval_prediction(X_eval,operator,prediction)
% function to obtain trajectories predicted using Koopman operator
% Inputs
% X             : Validation dataset
% operator      : Koopman matrix approximated using DMD or EDMD
% basis         : structure with properites of basis functions
% prediction    : sturecure with properties of trajectoreis for prediction
% Outputs
% X_pred        : Dataset of trajectories predicted using Koopman


% Project evaluation data onto DMD modes (Higher space)
b = pinv(operator.Phi)*X_eval(:,1);
size(b)

% Evolve modal coefficients forward in time using DMD eigenvalues
omega = log(diag(operator.D)) / prediction.dt;
t = linspace(0,prediction.n_steps * prediction.dt,prediction.n_steps);  % Total time span of the data
r = rank(operator.A);

X_eval_pred_projected = zeros(r, length(t));
for i = 1:length(t)
    X_eval_pred_projected(:, i) = exp(omega * t(i))*b;
end

% Reconstruct predicted states using DMD modes
X_pred = operator.Phi * X_eval_pred_projected;
X_pred(:,1) = X_eval(:,1);


end

