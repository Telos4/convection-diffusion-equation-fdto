% function of the terminal cost
function [l, gradl] = terminalcost(k,y)
    global epsilon;
    global n_y;
    y_ref = 0.5 * ones(n_y,1);
    Q = eye(n_y);
	l = 1/2 .* (epsilon .* (y - y_ref)' * Q * (y - y_ref));
    gradl = [epsilon*Q*(y-y_ref)];
end 