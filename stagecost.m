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
function [l, gradl] = stagecost(y,u,w, params)
    n_y = params.n_y;
    n_u = params.n_u;
    n_w = params.n_w;
    
    epsilon_y = params.epsilon_y;
    epsilon_u = params.epsilon_u;
    epsilon_w = params.epsilon_w;
    
    y_ref = params.y_ref;
    u_ref = params.u_ref;
    w_ref = params.w_ref;
    
    Q = epsilon_y * eye(n_y);
    R = epsilon_u * eye(n_u);
    W = epsilon_w * eye(n_w);
    
	l = 1/2 .* ((y - y_ref)' * Q * (y - y_ref) ... 
        + (u - u_ref)' * R * (u - u_ref) ...
        + (w - w_ref)' * W * (w - w_ref));
    
    gradl = [Q*(y-y_ref); R*(u-u_ref); W*(w-w_ref)];
end 