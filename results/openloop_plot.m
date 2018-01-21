clc;
close all;


filepart1 = 'results50/openloop_y';
file = strcat(filepart1, num2str(1), '.txt');
Z = load(file);

L = 200;
n_y = size(Z,2);
N = size(Z,1)-1;

X = (0:0.01:1);

x = [0 0 200 200];
y = [0.25 0.75 0.75 0.25];
z1 = [0.35 0.35 0.35 0.35];
z2 = [0.65 0.65 0.65 0.65];

yu = [0 0.01 0.01 0];
zu1 = [0.75 0.75 0.75 0.75];
zu2 = [0.25 0.25 0.25 0.25];



v = VideoWriter('heat.avi');
v.Quality = 100;
%v.LosslessCompression = true
v.FrameRate = 10;
open(v);


figure(1);
for i = 1:L


    file = strcat(filepart1, num2str(i), '.txt');     
    Z = load(file);
    Z = Z';
    Y = (i-1:i-1+N);

    surf(Y,X,Z)      
    
    hold on;
    axis([0 200 0 1 0 1])
  

    patch(x, y, z1, 'red', 'FaceAlpha',.3);
    patch(x, y, z2, 'red', 'FaceAlpha',.3);
    
    patch(x, yu, zu1, 'black', 'FaceAlpha', 0.5);
    patch(x, yu, zu2, 'black', 'FaceAlpha', 0.5);
    
    view(340, 10)
   
    xlabel('$x$ (timestep)','interpreter','latex');
    ylabel(['$y$ (Pos.)'],'interpreter','latex');
    zlabel('$z$ (Temp.)','interpreter','latex');
    
    
    set(gca,'TickLabelInterpreter','latex');
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    hold off;
end

close(v);
