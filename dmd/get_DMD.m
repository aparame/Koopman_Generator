function [A_tilde,Phi] = get_DMD(X1,X2,r)
% Load your data matrix X (each column represents a snapshot of the system)

% Assuming X is already loaded with appropriate data

% Perform Dynamic Mode Decomposition (DMD)
[U, E, V] = svd(X1, 'econ');


U_r = U(:, 1:r);
E_r = E(1:r, 1:r);
V_r = V(:, 1:r);

% STEP 2: low-rank subspace matrix
%         (similarity transform, least-square fit matrix, low-rank subspace matrix)
A_tilde = U_r' * X2 * V_r / E_r;

% Compute the DMD modes
[eigenvectors, ~] = eig(A_tilde);


Phi = X2 * V_r *E_r^(-1) * eigenvectors;

 
end