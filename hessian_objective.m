% Returns the hessian matrix of the objective function
% The objective is
%   J(k,x,u) = \sum_{j=1}^{N} epsilon * (x(j+1) - x_ref(k+j+1))^2 + u(j)^2
% We consider the optimization variables z = [x2 u1 x3 u2 ... xN+1 uN],
% Thus the hessian has the following structure:
% H =    [ e  0  ...                        ]
%        [ 0  1  0 ...                      ]
%        [ 0  0  e  0  ...                  ]
%        [ 0  0  0  1  0  ...               ]
%        [                   ...            ]
%        [                    ... 0  e  0   ]
%        [                      ...  0  1   ]
% The matrix is stored as a sparse matrix.
function [ H ] = hessian_objective( z, lambda, N )
global epsilon;
global n_y;
global n_u;
n_z = (N+1) * n_y + N * n_u;

inds = [1:(N+1)*n_y (N+1)*n_y+1:(N+1)*n_y+N*n_u];
vals = [epsilon*ones(1,(N+1)*n_y) ones(1,N*n_u)];
H = sparse(inds, inds, vals);
%H = speye(n_z);
end

