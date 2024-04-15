function [operator] = get_EDMD(X1, X2, U1,basis)
% function to generate A B matrices using EDMD
% Inputs
% X1            : state trajecotries of n traj of length 1:m-1
% X2            : state trajecotries of n traj of length 2:m
% basis         : structure with properites of basis functions
% Outputs
% A             : n x N lifted Koopman operator
% C             : matrix to map lifted states z to original states x
Z2 = [];
Z1 = [];
for i = 1:size(X1,2)
    psiX1 = get_basis(X1(:,i),basis);
    psiX2 = get_basis(X2(:,i),basis);
    Z1 = [Z1,psiX1];
    Z2 = [Z2,psiX2];
end

Z1_aug = [Z1;U1];


% % Perform Dynamic Mode Decomposition (DMD)
% [U1, E1, V1] = svd(Z1_aug, 'econ');
% [U11, E11, V11] = svd(Z1, 'econ');



% STEP 2: low-rank subspace matrix
%         (similarity transform, least-square fit matrix, low-rank subspace matrix)
A_edmd_aug = Z2*pinv(Z1);



C = X1*pinv(Z1);


A = A_edmd_aug(:,1:height(Z1));
B = A_edmd_aug(:, end);

operator.A = A;
operator.C = C;
operator.B = B;