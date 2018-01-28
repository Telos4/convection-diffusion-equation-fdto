function [ approx_cl_cost ] = closed_loop_simulation( params, N, L, plot_closed_loop )

% horizon length
params.N = N;

y_out = zeros(params.N+L,1);
for i = 0:params.N+L
    y_out(i+1) = 0.3 * sin(0.1 * i);
end

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

if plot_closed_loop == 1
    delete('heat.avi'); % delete video file to avoid errors
    v = VideoWriter('heat.avi');
    v.Quality = 100;
    v.FrameRate = 10;
    open(v);
    figure(1);

    ax = gca();

    for i = 1:L
        [lb,ub] = bounds(0+i,params.n_y);
        lb = params.lb_y * lb;
        ub = params.ub_y * ub;

        plot(ax,(1:params.n_y)/params.n_y, lb, 'r'); hold on;
        plot(ax,(1:params.n_y)/params.n_y, ub, 'r');
        plot(ax,[0.95 1.0], [0.25 0.25], 'b');
        plot(ax,[0.95 1.0], [-0.25 -0.25], 'b');
        plot(ax,1/(2*params.n_y):1/params.n_y:1, y_cl(i,:), 'k');
        plot(ax,0,y_out(i),'r*');
        plot(ax,1,u_cl(i),'b*');
        % draw arrows for convection velocity
    %     for j=1:2:params.n_y
    %         vstart = [j/params.n_y y_cl(i,j)];
    %         vend   = [j/params.n_y-0.001*sign(w_cl(i)) y_cl(i,j)];
    % %         arrow(vstart, vend)
    %     end
        axis([0 1 -0.5 0.5])

        xlabel('$x$ (Pos.)','interpreter','latex'); ylabel(['$y$ (Temp.)'],'interpreter','latex');
        set(ax,'TickLabelInterpreter','latex');
        text(0.8, 0.4, ['t = ' num2str(i/100)])

        if params.convection == 1
            text(0.8, 0.3, ['w = ' num2str(w_cl(i))])
        end
        hold off;

        drawnow();
        frame = getframe(gca);
        size(frame.cdata);
        writeVideo(v,frame);

    end

    close(v);
end

end

