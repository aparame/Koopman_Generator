function X_pred = eval_prediction(X_eval,X2,X1,operator,prediction)
% function to obtain trajectories predicted using Koopman operator
% Inputs
% X             : Validation dataset
% operator      : Koopman matrix approximated using DMD or EDMD
% basis         : structure with properites of basis functions
% prediction    : sturecure with properties of trajectoreis for prediction
% Outputs
% X_pred        : Dataset of trajectories predicted using Koopman

r = rank(operator.A);
[~, S, V] = svd(X1, 'econ');
Vr = V(:, 1:r);

Sr = S(1:r, 1:r);
[V_dmd, D] = eig(operator.A);
Phi = X2 * Vr / Sr * V_dmd;
size(Phi)
% Project evaluation data onto DMD modes
X_eval_projected = Phi' .* X_eval;

% Evolve modal coefficients forward in time using DMD eigenvalues
omega = log(diag(D)) / prediction.dt;
t = linspace(0,prediction.n_steps * prediction.dt,prediction.n_steps);  % Total time span of the data

X_eval_pred_projected = zeros(r, size(X_eval, 2));
for i = 1:prediction.n_steps
    X_eval_pred_projected(:, i) = exp(omega * t(i)) .* X_eval_projected(:, i);
end

% Reconstruct predicted states using DMD modes
X_pred = Phi * X_eval_pred_projected;





end

