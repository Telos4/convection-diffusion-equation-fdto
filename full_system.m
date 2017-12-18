clc;
close all;

global n_y;
global n_u;
global epsilon;

load('sys.mat')

A;
B_y;
b_u = b_u(:);
b_y_out = b_y_out(:);

n_y = size(A,1);
n_u = 1;
epsilon = 0.001;

k = 0;
y0 = 0.5 * ones(1, n_y);

N = 10;
L = 40;

y_out = zeros(N+L,1);
for i = 1:N+L
    y_out(i) = 0.5 + 0.3 * sin(0.1 * i);
end

y_cl = [y0];
for j = 1:L
    [u_ol, y_ol] = ocp_full_discretization(k+j, N, y0, A, B_y, b_u, b_y_out, y_out(j:j+N-1));
    
    u_ol(1)
    
    y0 = y_ol(2,:);
    
    y_cl = [y_cl; y0];
end

v = VideoWriter('heat.avi');
v.Quality = 100;
%v.LosslessCompression = true
v.FrameRate = 10;
open(v);
figure(1);
%F(L) = struct('cdata',[],'colormap',[]);
for i = 1:L
    [lb,ub] = bounds(k+i,n_y);
    plot((1:n_y)/n_y, lb, 'r'); hold on;
    plot((1:n_y)/n_y, ub, 'r');
    plot([0.95 1.0], [0.25 0.25], 'b');
    plot([0.95 1.0], [0.75 0.75], 'b');
    plot(1/(2*n_y):1/n_y:1, y_cl(i,:), 'k');
    axis([0 1 0 1])
    
    xlabel('$x$ (Pos.)','interpreter','latex'); ylabel(['$y$ (Temp.)'],'interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    text(0.8, 0.9, ['t = ' num2str(i/100)])
    
    frame = getframe(gcf);
    writeVideo(v,frame);
     hold off;
end

close(v);