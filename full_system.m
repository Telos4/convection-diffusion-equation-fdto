clc;
close all;
clear all;

load('sys.mat')

params = {};

% system matrices
params.A = A;
params.B_y = B_y;
params.B_w = B_w;
params.b_u = b_u(:);
params.b_y_out = b_y_out(:);

% dimensions
params.n_y = size(A,1);
params.n_u = 1;
params.n_w = 1;

% parameters for optimization
params.y_ref = 0.5 * ones(params.n_y,1);
params.u_ref = 0.5 * ones(params.n_u,1);
params.w_ref = 0.0 * ones(params.n_w,1);

% weights for objective function
params.epsilon_y = 0;%1.0e-6;
params.epsilon_u = 1.0;
params.epsilon_w = 1.0;

% horizon length
params.N = 10;

% initial time
params.k = 0;

% initial state
params.y0 = 0.5 * ones(1, params.n_y);




L = 100;
y_out = zeros(params.N+L,1);
for i = 0:params.N+L
    y_out(i+1) = 0.5 + 0.5 * sin(0.1 * i);
end



u_guess = zeros(params.N,1);
w_guess = zeros(params.N,1);

tic
y_cl = [params.y0];
u_cl = [];
w_cl = [];
for j = 0:L-1
    params.y_out = y_out(params.k+1:params.k+params.N);
    [u_ol, w_ol, y_ol] = ocp_full_discretization(params, u_guess, w_guess);
    
    disp(sprintf('j = %d, u = %f, w = %f', [j u_ol(1) w_ol(1)]));
    
    params.y0 = y_ol(2,:);
    params.k = params.k + 1;
    
    y_cl = [y_cl; fliplr(params.y0)];
    u_cl = [u_cl; u_ol(1)];
    w_cl = [w_cl; w_ol(1)];
    
    u_guess = [u_ol(2:end); u_ol(end)];
    w_guess = [w_ol(2:end); w_ol(end)];
end
toc
approx_cl_cost = sum((u_cl - params.u_ref).*(u_cl - params.u_ref)) + sum((w_cl - params.w_ref) .* (w_cl - params.w_ref));
disp(sprintf('approximate closed loop cost: %f', approx_cl_cost));
    

v = VideoWriter('heat.avi');
v.Quality = 100;
%v.LosslessCompression = true
v.FrameRate = 10;
open(v);
figure(1);

ax = gca();

for i = 1:L
    [lb,ub] = bounds(0+i,params.n_y);
    
    plot(ax,(1:params.n_y)/params.n_y, lb, 'r'); hold on;
    plot(ax,(1:params.n_y)/params.n_y, ub, 'r');
    plot(ax,[0.95 1.0], [0.25 0.25], 'b');
    plot(ax,[0.95 1.0], [0.75 0.75], 'b');
    plot(ax,1/(2*params.n_y):1/params.n_y:1, y_cl(i,:), 'k');
    plot(ax,0,y_out(i),'r*');
    plot(ax,1,u_cl(i),'b*');
    % draw arrows for convection velocity
    for j=1:2:params.n_y
        vstart = [j/params.n_y y_cl(i,j)];
        vend   = [j/params.n_y-0.001*sign(w_cl(i)) y_cl(i,j)];
%         arrow(vstart, vend)
    end
    axis([0 1 0 1])

    xlabel('$x$ (Pos.)','interpreter','latex'); ylabel(['$y$ (Temp.)'],'interpreter','latex');
    set(ax,'TickLabelInterpreter','latex');
    text(0.8, 0.9, ['t = ' num2str(i/100)])
    text(0.8, 0.7, ['w = ' num2str(w_cl(i))])
    hold off;
    
    drawnow();
    frame = getframe(gcf);
    size(frame.cdata);
    writeVideo(v,frame);
     
end

close(v);