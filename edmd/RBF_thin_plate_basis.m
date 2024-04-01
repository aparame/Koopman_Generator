function Psi = RBF_thin_plate_basis(X, centers)
%  [Psi, DPsi] = GaussianRBF_basis(deg, dim) returns a gaussian RBF function and its derivative 
%   There is not a 1 included in Psi, only the linear and higher order terms
%   Centers = matrix of center points (size = (dim of x, K))
%   dim = number of states

% x=sym('x',[dim,1]);
% assume(x,'real')

% Calculate distances between each data point and the centers
distances = pdist2(X', centers', 'euclidean');

% Compute Thin Plate RBF values
Psi = distances .^ 2 .* log(distances + eps); % Adding eps to avoid log(0)

end