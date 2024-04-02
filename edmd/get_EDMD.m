function [A,C,lambda] = get_EDMD(X1, X2, basis)
% function to generate A B matrices using EDMD
% Inputs
% X1            : state trajecotries of n traj of length 1:m-1
% X2            : state trajecotries of n traj of length 2:m
% basis         : structure with properites of basis functions
% Outputs
% A             : n x N lifted Koopman operator
% C             : matrix to map lifted states z to original states x

psiX1 = get_basis(X1,basis);
psiX2 = get_basis(X2,basis);

% Compute the Koopman operator using least squares
A = pinv(psiX1)*psiX2;
C = [X2;psiX2']*pinv(X1);
% Compute eigenvalues and eigenvectors of the Koopman operator
[V, D] = eig(A);
mu = diag(D);
% Sort eigenvalues in descending order
[eigenvalues_sorted, idx] = sort(mu, 'descend');
V_sorted = V(:, idx);

% Extract eigenvalues and eigenvectors corresponding to stable modes
lambda = diag(eigenvalues_sorted);
psi = psiX2 * V_sorted;
