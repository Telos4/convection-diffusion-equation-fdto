% Returns the hessian matrix of the objective function
% The objective is
%   J(k,x,u) = \sum_{j=0}^{N-1} epsilon * (y(j) - y_ref(j))^T * Q * (y(j) - y_ref(j)) 
%                                       + u(j)^T * R * u(j)
%                 + epsilon * (y(N) - y_ref(N))^T * Q * (y(N) - y_ref(N))
% We consider the optimization variables z = [y0 y1 ... yN u0 u1 ... uN-1],
% Thus the hessian has the following structure:
% H =    [ e*Q                        ]
%        [      ...                   ]
%        [           e*Q              ]
%        [                R           ]
%        [                   ...      ]
%        [                        R   ]
% The matrix is stored as a sparse matrix.
function [ H ] = hessian_objective( z, lambda, N )
global epsilon;
global n_y;
global n_u;
n_z = (N+1) * n_y + N * n_u;

inds = [1:(N+1)*n_y (N+1)*n_y+1:(N+1)*n_y+N*n_u];
vals = [epsilon*ones(1,(N+1)*n_y) ones(1,N*n_u)];
H = sparse(inds, inds, vals);
end

