x = load('closed_cost_vs_horizon_length.txt');
figure(1);

hold on;
plot(x)
xlabel('$N$ (MPC horizon)','interpreter','latex');
ylabel(['$F$ (closed loop cost)'],'interpreter','latex');

hold off;