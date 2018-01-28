% function of the terminal cost
function [l, gradl] = terminalcost(y, params)
    n_y = params.n_y;
    epsilon_y = params.epsilon_y;
    
    y_ref = params.y_ref;
    
    Q = epsilon_y * eye(n_y);
    
	l = 1/2 .* ((y - y_ref)' * Q * (y - y_ref));
    gradl = [Q*(y-y_ref)];
end 