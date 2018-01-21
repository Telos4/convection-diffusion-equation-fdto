x = load('closed_cost_vs_steps_20.txt');

figure(2);

hold on;
plot(x)
xlabel('steps','interpreter','latex');
ylabel(['total closed loop cost'],'interpreter','latex');
title('closed loop cost vs steps')
hold off;