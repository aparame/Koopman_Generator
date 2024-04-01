function EDMD = get_EDMD(X1, X2, basis)
% function to generate A B matrices using EDMD
% Inputs
% X1            : state trajecotries of n traj of length 1:m-1
% X2            : state trajecotries of n traj of length 2:m
% basis         : structure with properites of basis functions
% Outputs
% A             : n x N lifted Koopman operator
% C             : matrix to map lifted states z to original states x

psi = get_basis(X1,basis);
% Compute the Koopman operator using least squares
A = psi\Y;

% Compute eigenvalues and eigenvectors of the Koopman operator
[V, D] = eig(A);

% Sort eigenvalues in descending order
[eigenvalues_sorted, idx] = sort(diag(D), 'descend');
V_sorted = V(:, idx);

% Extract eigenvalues and eigenvectors corresponding to stable modes
Lambda = diag(eigenvalues_sorted);
psi = psi * V_sorted;
