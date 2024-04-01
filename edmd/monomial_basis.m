function monomials = monomial_basis(x, deg, dim, remove_lin)
%  [Psi, DPsi] = monomial_basis(deg, dim) returns a monomial basis function and its derivative 
%   There is not a 1 included in Psi, only the linear and higher order terms
%   deg = degree of monomial
%   dim = number of states
k = linspace(2, deg, deg-1);
d = dim;

if (deg == 0)
    Psi = zeros(dim, 1);
elseif (deg == 1)
    Psi = x.'
else
    Psi = x.';
    for i = 1:size(k,2)
        m = nchoosek(k(i)+d-1,d-1); 
        dividers = [zeros(m,1),nchoosek((1:(k(i)+d-1))',d-1),ones(m,1)*(k(i)+d)]; 
        a = diff(dividers,1,2)-1;
        size(a)
        size(x.')
        for j = 1:size(a,1)
            Psi = [Psi prod(x.' .^ a(j,:))];
        end
    end
    if(remove_lin)
        Psi = Psi(dim+1:end);
    end
end
monomials = Psi';

end