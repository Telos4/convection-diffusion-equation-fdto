load('sys.mat')

A
B_y
b_u
b_y_out
u
y_out

n_y = size(A,1);

M = zeros(n_y*N);
b = zeros(n_y*N,1);

% assemble system matrix
for i = 1:N
    M((i-1)*n_y+1:i*n_y,(i-1)*n_y+1:i*n_y) = A;
end
for i = 2:N
    M((i-1)*n_y+1:i*n_y,(i-2)*n_y+1:(i-1)*n_y) = -B_y;
end

% assemble rhs
for i = 1:N
    b((i-1)*n_y+1:i*n_y) = b_u * u(i) + b_y_out * y_out(i);
end

% solve system
y = M \ b;

ys = [];
for i = 1:N
    ys = [ys; y((i-1)*n_y+1:i*n_y)'];
end
ys