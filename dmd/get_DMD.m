function DMD = get_DMD(X1,X2)
% Load your data matrix X (each column represents a snapshot of the system)

u0 = X1(:, 1);
% Assuming X is already loaded with appropriate data
dt = 0.1;
% Perform Dynamic Mode Decomposition (DMD)
[U, S, V] = svd(X1, 'econ');

r = 10;
U_r = U(:, 1:r);
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);

% STEP 2: low-rank subspace matrix
%         (similarity transform, least-square fit matrix, low-rank subspace matrix)
A_tilde = U_r' * X2 * V_r / S_r;

% Compute the DMD modes
[eigenvectors, D] = eig(A_tilde);
lambda = diag(D);       % eigen value
omega = log(lambda)/dt; % log of eigen value

Phi = X2 * V_r * S_r^(-1) * eigenvectors;

% STEP 5: reconstruct the signal
b = pinv(Phi)*u0;  % pseudo-inverse initial conditions
t = 11;
u_modes = zeros(r,length(t));  % DMD reconstruction for every time point
for i = 1:length(t)
    u_modes(:,i) =(b.*exp(omega*(t(i))));
end
% Construct the Koopman operator
DMD = Phi * u_modes;

% Now you have the finite-dimensional approximation for the Koopman operator in K

   
end