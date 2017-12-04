% function of the stage cost
% We consider the following stage cost:
%   l(k, x, u) = eps * (x - x_ref(k))^2  +  u^2
%
% The function evaluates the stage cost at the current time k with the
% current state x and control u.
% The parameter epsilon defines the weighting of the state penalization
% from the reference trajectory x_ref.
% The function also returns the gradient of the stage cost w.r.t. the state
% and the control (this is used in the optimization)
function [l, gradl] = stagecost(k,y,u)
    global epsilon;
    global n_y;
    global n_u;
    y_ref = 0.5 * ones(n_y,1);
    u_ref = 0.5 * ones(n_u,1);
    Q = eye(n_y);
    R = eye(n_u);
	l = 1/2 .* (epsilon .* (y - y_ref)' * Q * (y - y_ref) + (u - u_ref)' * R * (u - u_ref));
    
    gradl = [epsilon*Q*(y-y_ref); R*(u-u_ref)];
end 