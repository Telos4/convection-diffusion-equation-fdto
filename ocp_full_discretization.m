
% solve the optimal control problem using a full discretization scheme
%   k   - current time index
%   N   - horizon length
%   x0  - initial value at current time
%
% The full discretization scheme treats both states and controls as
% optimization variables. In total there are 2*N optimization variables
% that correspond to the following states/controls:
%   z = [x0 u0 x1 u1 ... xN-1 uN-1]
%
% The system dynamic is implemented as a linear constraint to the 
% optimization problem. 
% The time-varying constraints for the state and the constraints for the
% control are implemented as simple lower and upper bounds (box
% constraints) of the optimization variable z.
function [u_ol, y_ol, fval,exitflag] = ocp_full_discretization(k, N, y0, A, B_y, b_u, b_y_out, y_outs)
global free_initial_value;
global n_y;
global n_u;

n_z = (N+1) * n_y + N * n_u;

% initial value
z = zeros(n_z,1);

% system dynamic as equality constraints
[Aeq, beq] = assemble_eq_constraints(N, A, B_y, b_u, b_y_out, y_outs);

% no linear inequality constraints
A = [];
b = [];

% define state and control constraints
[lb, ub] = assemble_state_control_constraints(k, N, y0);

% set option to use derivative information
options = optimset('Algorithm', 'interior-point', 'GradObj', 'on', 'Hessian', ...
   'user-supplied', 'HessFcn', @(z,lambda)hessian_objective(z,lambda,N), ...
   'Display', 'off');

tic
% run optimization
[z,fval,exitflag, output] = fmincon(@(z) objective(k,N,z), z, A, b, Aeq, beq, lb, ub, [], options);
toc

% extract optimal control sequence from result of the optimization
y_ol = zeros(N+1, n_y);
u_ol = zeros(N, n_u);
for j = 0:N
    y_ol(j+1,:) = z(j*n_y+1:(j+1)*n_y)';
end
for j = 0:N-1
    u_ol(j+1,:) = z((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u)';
end
end

% objective function
function [J, gradJ, hessJ] = objective(k, N, z)
    global n_y;
    global n_u;
    J = 0;
    gradJ = zeros(size(z));

    for j = 0:N-1
        yj = z(j*n_y+1:(j+1)*n_y);
        uj = z((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u);
        [l, gradl] = stagecost(k+j,yj,uj);
        J = J + l;
        gradJ(j*n_y+1:(j+1)*n_y) = gradl(1:n_y);
        gradJ((N+1)*n_y+j*n_u+1:(N+1)*n_y+(j+1)*n_u) = gradl(n_y+1:n_y+n_u);
    end
    
    % cost at final time
    yN = z(N*n_y+1:(N+1)*n_y);
    [l, gradl] = terminalcost(k, yN);
    J = J + l;
    gradJ(N*n_y+1:(N+1)*n_y) = gradl(1:n_y);
end

% The system dynamic is defined as an equality constraint for the 
% optimization problem.
% The matrix looks like this:
% A_eq =   [ -B_y    A                     b_u              ]
%          [      -B_y    A                    b_u          ]
%          [                 ...                    ...     ]
%          [                   -B_y    A                b_u ]
% The rhs vector contains the following entries:
% b_eq = [    b_y_out * y_out(0)   ]
%        [    b_y_out * y_out(2)   ]
%        [           ...           ]
%        [    b_y_out * y_out(N-1) ]
function [Aeq, beq] = assemble_eq_constraints(N, A, B_y, b_u, b_y_out, y_out)
    global n_y;
    global n_u;

    Aeq = zeros(N * n_y, (N+1) * n_y + N * n_u);
    for i = 1:N
        Aeq((i-1)*n_y+1:i*n_y,i*n_y+1:(i+1)*n_y) = A;
        Aeq((i-1)*n_y+1:i*n_y,(i-1)*n_y+1:i*n_y) = -B_y;
        
        Aeq((i-1)*n_y+1:i*n_y, (N+1)*n_y+(i-1)*n_u+1:(N+1)*n_y+i*n_u) = -b_u;
    end
    Aeq = sparse(Aeq);  % convert to sparse matrix

    % rhs is time varying data and initial value of x0
    beq = zeros(N * n_y, 1);
    for i = 1:N
        beq((i-1)*n_y+1:i*n_y) = b_y_out * y_out(i);
    end
end

function [lb, ub] = assemble_state_control_constraints(k, N, y0)
    global n_y;
    global n_u;
    n_z = (N+1) * n_y + N * n_u;
    
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
    lb((N+1)*n_y+1:(N+1)*n_y+N*n_u) = 0.25;
    ub((N+1)*n_y+1:(N+1)*n_y+N*n_u) = 0.75;
end