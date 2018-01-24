
% solve the optimal control problem using a full discretization scheme
%   k   - current time index
%   N   - horizon length
%   x0  - initial value at current time
%
% The full discretization scheme treats both states and controls as
% optimization variables. In total there are 2*N optimization variables
% that correspond to the following states/controls:
%   z = [y0 y1 ... yN u0 u1 ... uN-1 w0 w1 ... wN-1 v0 v1 ... vN-1]
%
% The system dynamic is implemented as a linear constraint to the 
% optimization problem. 
% The time-varying constraints for the state and the constraints for the
% control are implemented as simple lower and upper bounds (box
% constraints) of the optimization variable z.
function [u_ol, w_ol, y_ol, fval,exitflag] = ocp_full_discretization(params)
% extract parameters
n_y = params.n_y;
n_u = params.n_u;
n_w = params.n_w;
n_v = params.n_v;
N = params.N;

n_z = (N+1) * n_y + N * n_u + N * n_w + N * n_v;

% initial value
z = zeros(n_z,1); 

% system dynamic as equality constraints
[Aeq, beq] = assemble_eq_constraints(params);

% no linear inequality constraints
A = [];
b = [];

% define state and control constraints
[lb, ub] = assemble_state_control_constraints(params);

% set option to use derivative information
options = optimset('Algorithm', 'interior-point', ...
    'GradObj', 'on', ...
    'GradConstr', 'on', ...
    'DerivativeCheck', 'off', ...
    'Hessian', 'user-supplied', 'HessFcn', @(z,lambda)hessian_objective(z,lambda,params), ...
    'Display', 'off');

%tic
% run optimization
[z,fval,exitflag, output] = fmincon(@(z) objective(z, params), z, A, b, Aeq, beq, lb, ub, @(z) nonlinear_constraints(z, params), options);
%toc

% if exitflag == 0
%     save('test_data', 'z')
% end

% extract optimal control sequence from result of the optimization
y_ol = zeros(N+1, n_y);
u_ol = zeros(N, n_u);
w_ol = zeros(N, n_w);
for j = 0:N
    y_ol(j+1,:) = z(j*n_y+1:(j+1)*n_y)';
end
for j = 0:N-1
    u_ol(j+1,:) = z((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u)';
end
for j = 0:N-1
    w_ol(j+1,:) = z((N+1)*n_y+N*n_u+j*n_w+1:(N+1)*n_y+N*n_u+(j+1)*n_w)';
end

% checks
if exitflag ~= 1 
    disp(sprintf('WARNING: exitflag = %d', exitflag))
end
if norm(Aeq * z - beq) > 1.0e-3
    disp(sprintf('WARNING: |Aeq * z - beq| = %f', norm(Aeq * z - beq)))
end
end

% objective function
function [J, gradJ, hessJ] = objective(z, params)
    n_y = params.n_y;
    n_u = params.n_u;
    n_w = params.n_w;
    N = params.N;
    J = 0;
    gradJ = zeros(size(z));

    for j = 0:N-1
        yj = z(j*n_y+1:(j+1)*n_y);
        uj = z((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u);
        wj = z((N+1)*n_y+N*n_u+j*n_w+1:(N+1)*n_y+N*n_u+(j+1)*n_w);
        [l, gradl] = stagecost(yj,uj,wj, params);
        J = J + l;
        gradJ(j*n_y+1:(j+1)*n_y) = gradl(1:n_y);
        gradJ((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u) = gradl(n_y+1:n_y+n_u);
        gradJ((N+1)*n_y+N*n_u+j*n_w+1:(N+1)*n_y+N*n_u+(j+1)*n_w) = gradl(n_y+n_u+1:n_y+n_u+n_w);
        % v doesn't appear in cost function so we don't need to compute the
        % gradient
    end
    
    % cost at final time
    yN = z(N*n_y+1:(N+1)*n_y);
    [l, gradl] = terminalcost(yN, params);
    J = J + l;
    gradJ(N*n_y+1:(N+1)*n_y) = gradl(1:n_y);
end

% The system dynamic is defined as an equality constraint for the 
% optimization problem.
% The matrix looks like this:
% A_eq =   [ -B_y    A                     -b_u               0             I       ]
%          [      -B_y    A                    -b_u             0             I     ]
%          [                 ...                     ...          ...         ...   ]
%          [                   -B_y    A                -b_u          0           I ]
% The rhs vector contains the following entries:
% b_eq = [    b_y_out * y_out(0)   ]
%        [    b_y_out * y_out(2)   ]
%        [           ...           ]
%        [    b_y_out * y_out(N-1) ]
function [Aeq, beq] = assemble_eq_constraints(params)
    n_y = params.n_y;
    n_u = params.n_u;
    n_w = params.n_w;
    n_v = params.n_v;
    N = params.N;
    A = params.A;
    B_y = params.B_y;
    b_u = params.b_u;
    b_y_out = params.b_y_out;
    y_out = params.y_out;

    n_z = (N+1) * n_y + N * n_u + N * n_w + N * n_v;
    
    Aeq = zeros(N * n_y, n_z);
    for i = 1:N
        Aeq((i-1)*n_y+1:i*n_y,i*n_y+1:(i+1)*n_y) = A;
        Aeq((i-1)*n_y+1:i*n_y,(i-1)*n_y+1:i*n_y) = -B_y;        
        Aeq((i-1)*n_y+1:i*n_y, (N+1)*n_y+(i-1)*n_u+1:(N+1)*n_y+i*n_u) = -b_u;
        Aeq((i-1)*n_y+1:i*n_y, (N+1)*n_y+N*n_u+N*n_w+(i-1)*n_y+1:(N+1)*n_y+N*n_u+N*n_w+i*n_v) = eye(n_v);
    end
    Aeq = sparse(Aeq);  % convert to sparse matrix

    % rhs is time varying data and initial value of x0
    beq = zeros(N * n_y, 1);
    for i = 1:N
        beq((i-1)*n_y+1:i*n_y) = b_y_out * y_out(i);
    end
end

function [c, ceq, grad_c, grad_ceq] = nonlinear_constraints(z, params)
    n_y = params.n_y;
    n_u = params.n_u;
    n_w = params.n_w;
    n_v = params.n_v;
    N = params.N;
    B_w = params.B_w;
    
    n_z = (N+1) * n_y + N * n_u + N * n_w + N * n_v;
    
    c = [];
    ceq = zeros(N * n_v, 1);
    for i = 1:N
        wk = z((N+1)*n_y+N*n_u+(i-1)*n_w+1:(N+1)*n_y+N*n_u+i*n_w);
        vk = z((N+1)*n_y+N*n_u+N*n_w+(i-1)*n_v+1:(N+1)*n_y+N*n_u+N*n_w+i*n_v);
        yk1 = z(i*n_v+1:(i+1)*n_v);
        ceq((i-1)*n_v+1:i*n_v) = vk - wk' * B_w * yk1;
    end
    
    grad_c = [];
    
    grad_ceq = zeros(N * n_v, n_z);
    for i = 1:N
        wk = z((N+1)*n_y+N*n_u+(i-1)*n_w+1:(N+1)*n_y+N*n_u+i*n_w);
        yk1 = z(i*n_v+1:(i+1)*n_v);
        
        grad_ceq((i-1)*n_v+1:i*n_v,i*n_v+1:(i+1)*n_v) = -wk' * B_w;
        grad_ceq((i-1)*n_v+1:i*n_v,(N+1)*n_y+N*n_u+(i-1)*n_w+1:(N+1)*n_y+N*n_u+i*n_w) = -B_w * yk1;
        grad_ceq((i-1)*n_v+1:i*n_v,(N+1)*n_y+N*n_u+N*n_w+(i-1)*n_v+1:(N+1)*n_y+N*n_u+N*n_w+i*n_v) = eye(n_v);
    end
    grad_ceq = sparse(grad_ceq');
end

function [lb, ub] = assemble_state_control_constraints(params)
    n_y = params.n_y;
    n_u = params.n_u;
    n_w = params.n_w;
    n_v = params.n_v;
    N = params.N;
    y0 = params.y0;
    k = params.k;
    
    n_z = (N+1) * n_y + N * n_u + N * n_w + N * n_v;
    
    lb = -Inf * ones(n_z, 1);
    ub =  Inf * ones(n_z, 1);
    
    lb_y = -1.0;
    ub_y =  1.0;
    
    % set initial condition
    lb(1:n_y) = y0;
    ub(1:n_y) = y0;
    
    % constraints for state
    for i = 2:N+1
        [lb_y, ub_y] = bounds(k, n_y);
        lb(n_y * (i-1)+1:n_y*i) = lb_y;
        ub(n_y * (i-1)+1:n_y*i) = ub_y;
    end

    % constraints for the control
    lb((N+1)*n_y+1:(N+1)*n_y+N*n_u) = -0.25;
    ub((N+1)*n_y+1:(N+1)*n_y+N*n_u) = 0.25;
end