function dxdt = dynamics_duffing_MPC(t,x,U_log,T_log)
    u = interp1(T_log(2:end), U_log, t);
    delta = 0.5; 
    scaling = 1;
    dxdt = zeros(2, 1);
    dxdt(1) = x(2);
    dxdt(2) = u + x(1) - delta*x(2) - scaling.*x(1)^3;
end
