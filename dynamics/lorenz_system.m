function dxdt = lorenz_system(t,x,u)
    % lorenz_system: Computes the derivative of the Lorenz system
    % Inputs:
    %   x: State vector [x, y, z]
    %   t: Time (not used in the Lorenz system, included for compatibility)
    % Output:
    %   dxdt: Derivative of the state vector
    
    % Lorenz system parameters
    sigma = 10;
    rho = 28;
    beta = 8/3;
    
    % Compute derivatives
    dxdt = zeros(3, 1);
    dxdt(1) = sigma * (x(2) - x(1));
    dxdt(2) = x(1) * (rho - x(3)) - x(2) + u;
    dxdt(3) = x(1) * x(2) - beta * x(3);
end