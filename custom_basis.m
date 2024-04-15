function Psi = custom_basis(x)
%  [Psi, DPsi] = GaussianRBF_basis(deg, dim) returns a gaussian RBF function and its derivative
%   There is not a 1 included in Psi, only the linear and higher order terms
%   Centers = matrix of center points (size = (dim of x, K))
%   dim = number of states

% x=sym('x',[dim,1]);
% assume(x,'real')

Psi = [];
G = [];

Psi = [x(1),x(2),x(3),x(1)*x(2),x(1)*x(3)]';

end