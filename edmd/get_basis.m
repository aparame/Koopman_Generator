function psi = get_basis(X, basis)
% Input
% X         : states
% n         : number of basis
% output
% basis     : basis functions

%% build observables [x,basis]
psi = []; %  basis function [linear, nonlienar]
% psi = custom_basis(X);
if(strcmp(basis.type,'monomials'))
    deg = basis.deg; % degeree of monomial
    dim = basis.dim; % state dimension
    if(~isfield(basis,'remove_lin'))
        remove_lin = false;
    else
        remove_lin = basis.remove_lin;
    end
    psi = monomial_basis(X,deg,dim,remove_lin);

elseif(strcmp(basis.type,'rbf'))
    gamma = basis.gamma;
    % 
    % % random sampling
    dom = [-2,2];
    c1 = dom(1) + (dom(1)-dom(2))*rand(150,1);
    c2 = dom(1) + (dom(1)-dom(2))*rand(150,1);
    c3 = dom(1) + (dom(1)-dom(2))*rand(150,1);
    centers = [c1(:),c2(:),c3(:)];

    % % unifom sampling
    % % grid_x= -2:0.2:2;
    % % grid_y= -2:0.2:2;
    % % grid_z= -2:0.2:2;
    % % [c1,c2,c3] = meshgrid(grid_x,grid_y,grid_z);
    % % centers = [c1(:),c2(:),c3(:)];

    % gaussian rbf
    psi =RBF_gaussian_basis(X, centers, gamma);
    % thin plate rbf
    % psi = RBF_thin_plate_basis(X, centers);

elseif(strcmp(basis.type,'custom'))
     psi = custom_basis(X);
end