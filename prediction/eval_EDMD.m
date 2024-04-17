function X_pred = eval_EDMD(X_eval,U,basis,operator,prediction)
% function to obtain trajectories predicted using Koopman operator
% Inputs
% X             : Validation dataset
% operator      : Koopman matrix approximated using DMD or EDMD
% basis         : structure with properites of basis functions
% prediction    : sturecure with properties of trajectoreis for prediction
% Outputs
% X_pred        : Dataset of trajectories predicted using Koopman


t = linspace(0,prediction.n_steps * prediction.dt,prediction.n_steps);  % Total time span of the data
X0 = X_eval(:,1);
z0 = get_basis(X0,basis);
z = z0;
Z_pred = [];

for i = 1:length(t)
    z_next = operator.A*z + operator.B*U(:,i);
    z = z_next;
    Z_pred = [Z_pred,z];
end
X_pred = operator.C*Z_pred;
figure(2)
hold on
p1 = plot(X_pred(1,:),X_pred(2,:),'DisplayName','X_{pred}');
p2 = plot(X_eval(1,1:length(t)),X_eval(2,1:length(t)),'DisplayName','X_{eval}');
legend([p1,p2])
hold off

end

