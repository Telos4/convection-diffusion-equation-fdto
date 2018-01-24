clc;
close all;
clear all;

load('sys.mat')

params = {};

% system matrices
params.A = A;
params.B_y = B_y;
params.b_u = b_u(:);
params.b_y_out = b_y_out(:);

% dimensions
params.n_y = size(A,1);
params.n_u = 1;
params.n_w = 0;
params.n_v = 0;%params.n_y;

% parameters for optimization
params.y_ref = 0.5 * ones(params.n_y,1);
params.u_ref = 0.5 * ones(params.n_u,1);

% weights for objective function
params.epsilon_y = 0.0;
params.epsilon_u = 1.0;

% additional parameters for control by convection term
params.convection = 1;
if params.convection == 1
    params.B_w = B_w;
    params.n_w = 1;
    params.n_v = params.n_y;
    params.w_ref = 0.0 * ones(params.n_w,1);
    params.epsilon_w = 1.0;
else
    params.n_w = 0;
    params.n_v = 0;
    params.B_w = zeros(params.n_v, params.n_v); % empty matrix
    params.w_ref = 0.0 * ones(params.n_w,1);    % empty matrix
    params.epsilon_w = 0.0;
end

% horizon length
params.N = 2;

% initial time
params.k = 0;

% initial state
params.y0 = 0.5 * ones(1, params.n_y);


L = 50;
y_out = zeros(params.N+L,1);
for i = 0:params.N+L
    y_out(i+1) = 0.5 + 0.3 * sin(0.1 * i);
end

u_guess = zeros(params.N,1);
w_guess = zeros(params.N,1);

tic
y_cl = [params.y0];
u_cl = zeros(L, params.n_u);
w_cl = zeros(L, params.n_w);
for j = 0:L-1
    params.y_out = y_out(params.k+1:params.k+params.N);
    [u_ol, w_ol, y_ol] = ocp_full_discretization(params);
    
    if params.convection == 1
        disp(sprintf('j = %d, u = %f, w = %f', [j u_ol(1) w_ol(1)]));
    else
        disp(sprintf('j = %d, u = %f', [j u_ol(1)]));
    end
    
    params.y0 = y_ol(2,:);
    params.k = params.k + 1;
    
    y_cl = [y_cl; fliplr(params.y0)];
    u_cl(j+1) = u_ol(1);
    if params.convection == 1
       w_cl(j+1) = w_ol(1);
    end
end
toc

if params.convection == 1
    approx_cl_cost = sum((u_cl - params.u_ref).*(u_cl - params.u_ref)) + sum((w_cl - params.w_ref) .* (w_cl - params.w_ref));
else
    approx_cl_cost = sum((u_cl - params.u_ref).*(u_cl - params.u_ref));
end
disp(sprintf('approximate closed loop cost: %f', approx_cl_cost));
    

delete('heat.avi'); % delete video file to avoid errors
v = VideoWriter('heat.avi');
v.Quality = 100;
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
%     for j=1:2:params.n_y
%         vstart = [j/params.n_y y_cl(i,j)];
%         vend   = [j/params.n_y-0.001*sign(w_cl(i)) y_cl(i,j)];
% %         arrow(vstart, vend)
%     end
    axis([0 1 0 1])

    xlabel('$x$ (Pos.)','interpreter','latex'); ylabel(['$y$ (Temp.)'],'interpreter','latex');
    set(ax,'TickLabelInterpreter','latex');
    text(0.8, 0.9, ['t = ' num2str(i/100)])
    
    if params.convection == 1
        text(0.8, 0.7, ['w = ' num2str(w_cl(i))])
    end
    hold off;
    
    drawnow();
    frame = getframe(gca);
    size(frame.cdata);
    writeVideo(v,frame);
     
end

close(v);