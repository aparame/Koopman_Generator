function DMD = get_DMD(X1,X2)
% Load your data matrix X (each column represents a snapshot of the system)
% Assuming X is already loaded with appropriate data

% Perform Dynamic Mode Decomposition (DMD)
[U, S, V] = svd(X1, 'econ');
r = rank(X1); % Choose the rank of the approximation (usually determined by energy content)

U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

A_tilde = U_r' *X2 * V_r / S_r;

% Compute the DMD modes
[eigenvectors, ~] = eig(A_tilde);
Phi = X2 * V_r / S_r * eigenvectors;

% Construct the Koopman operator
DMD = Phi / X1;

   
end