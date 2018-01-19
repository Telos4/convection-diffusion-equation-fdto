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
%        [                        R              ]
%        [                           W           ]
%        [                              ...      ]
%        [                                  W    ]
% The matrix is stored as a sparse matrix.
function [ H ] = hessian_objective( z, lambda, params )
n_y = params.n_y;
n_u = params.n_u;
n_w = params.n_w;
N = params.N;
    
epsilon_y = params.epsilon_y;
epsilon_u = params.epsilon_u;
epsilon_w = params.epsilon_w;

n_z = (N+1) * n_y + N * n_u + N * n_w + N * n_y;

inds = [1:(N+1)*n_y (N+1)*n_y+1:(N+1)*n_y+N*n_u (N+1)*n_y+N*n_u+1:(N+1)*n_y+N*n_u+N*n_w];
vals = [epsilon_y*ones(1,(N+1)*n_y) epsilon_u*ones(1,N*n_u) epsilon_w*ones(1,N*n_w)];
H = sparse(inds, inds, vals, n_z, n_z);
end

