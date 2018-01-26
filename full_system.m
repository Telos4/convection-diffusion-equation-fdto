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

% parameters for optimization
params.y_ref = 0.0 * ones(params.n_y,1);
params.u_ref = 0.0 * ones(params.n_u,1);

% weights for objective function
params.epsilon_y = 0.0;
params.epsilon_u = 1.0;

params.lb_u = -0.25;
params.ub_u =  0.25;

params.lb_y = -0.15;
params.ub_y = 0.15;

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

% initial time
params.k = 0;

% initial state
y0 = 0.0 * ones(1, params.n_y);

Ns = [1:1:25];
L = 100;
J_cls = [];
for i = 1:length(Ns)
    N = Ns(i)
    params.y0 = y0;
    cl = closed_loop_simulation(params, N, L);
    J_cls = [J_cls; cl];
end

figure;
plot(Ns, J_cls)
xlim([Ns(1) Ns(end)])
xlabel('N','interpreter','latex'); ylabel(['$J^{cl}_{' num2str(L) '}(0, 0, \mu_N)$'],'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gcf, 'PaperPosition', [0 0 15 7]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 7]);
saveas(gcf, 'closed-loop-cost', 'pdf')
