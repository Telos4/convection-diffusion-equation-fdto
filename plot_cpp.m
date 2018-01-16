clc;
close all;

n_y = 101;
L = 200;
y = load('solution_y.txt');
y = fliplr(y);
u = load('solution_u.txt');


v = VideoWriter('heat.avi');
v.Quality = 100;
%v.LosslessCompression = true
v.FrameRate = 10;
open(v);
figure(1);
%F(L) = struct('cdata',[],'colormap',[]);
for i = 1:L
    plot([0.25 0.75], [0.35 0.35], 'r');     hold on;
    plot([0.25 0.75], [0.65 0.65], 'r');
    plot([0.95 1.0], [0.25 0.25], 'b');
    plot([0.95 1.0], [0.75 0.75], 'b');
    %plot(1/(2*n_y):1/n_y:1, y(:,i), 'k');
    plot(linspace(0, 1, n_y), y(i,:));
    plot(1, u(i), '-o');
    axis([0 1 0 1])
    xlabel('$x$ (Pos.)','interpreter','latex'); ylabel(['$y$ (Temp.)'],'interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    text(0.8, 0.9, ['t = ' num2str(i/100)])
    
    frame = getframe(gcf);
    writeVideo(v,frame);
     hold off;
end

close(v);