function dxdt = dynamics_duffing(t,x,u)
    delta = 0.5; 
    scaling = 1;
    dxdt = zeros(2, 1);
    dxdt(1) = x(2);
    dxdt(2) = u + x(1) - delta*x(2) - scaling.*x(1)^3;
end
