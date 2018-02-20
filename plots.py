"""
===========
MovieWriter
===========

This example uses a MovieWriter directly to grab individual frames and write
them to a file. This avoids any event loop integration, but has the advantage
of working with even the Agg backend. This is not recommended for use in an
interactive setting.

"""
# -*- noplot -*-

import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

from subprocess import call
from pathlib import Path

lb_y = -0.15
ub_y =  0.15

lb_u = -0.25
ub_u =  0.25

delta_t = 1.0e-2
linewidth_global = 2.0

class SimulationResult:
    def __init__(self, results_folder):
        p = Path(results_folder)

        q = p / 'closedloop_y.txt'
        if q.exists():
            self.y_cl = np.fliplr(np.loadtxt(str(q), ndmin=2))

        q = p / 'closedloop_u.txt'
        if q.exists():
            self.u_cl = np.loadtxt(str(q), ndmin=1)

        q = p / 'closedloop_w.txt'
        if q.exists():
            self.w_cl = np.loadtxt(str(q), ndmin=1)

        q = p / 'closedloop_cost.txt'
        if q.exists():
            self.closed_loop_cost = np.loadtxt(str(q), ndmin=1)
            self.final_closed_loop_cost = self.closed_loop_cost[-1]

        # load open loops
        open_loops_y = list(p.glob('openloop_y*'))
        if len(open_loops_y) > 0:
            self.y_ol = []
        for i in range(0, len(open_loops_y)):
            q = p / ('openloop_y'+str(i)+'.txt')
            if q.exists():
                self.y_ol.append(np.fliplr(np.loadtxt(str(q), ndmin=1)))

        open_loops_u = list(p.glob('openloop_u*'))
        if len(open_loops_u) > 0:
            self.u_ol = []
        for i in range(0, len(open_loops_u)):
            q = p / ('openloop_u'+str(i)+'.txt')
            if q.exists():
                self.u_ol.append(np.loadtxt(str(q), ndmin=1))

        open_loops_w = list(p.glob('openloop_w*'))
        if len(open_loops_w) > 0:
            self.w_ol = []
        for i in range(0, len(open_loops_w)):
            q = p / ('openloop_w'+str(i)+'.txt')
            if q.exists():
                self.w_ol.append(np.loadtxt(str(q), ndmin=1))

        # find parameters
        try:
            self.L = len(self.closed_loop_cost)
        except AttributeError:
            print("closed loop not available")

        try:
            self.N = len(self.u_ol[0])      # MPC horizon
            self.n_y = len(self.y_ol[0][0]) # number of finite element nodes
        except AttributeError:
            print("open loop not available")



    def plot_closed_loop(self, output_file, reference=None, L=None):
        if L is None:
            L = self.L
        else:
            L = min(L, self.L)
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='closed_loop_N='+str(self.N), artist='Matplotlib',
                    comment='Movie support!')
        writer = FFMpegWriter(fps=15, metadata=metadata)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        h = 1.0/(self.n_y - 1)
        domain = np.arange(0.0, 1.0+h, h)
        subdomain = np.arange(0.25, 0.75+h, h)
        fig, ax = plt.subplots()
        #ax.hold(False)

        with writer.saving(fig, output_file, 100):
            for i in range(0, L):
                ax.plot(domain, self.y_cl[i], color='k', linewidth=linewidth_global, label='$y\;(N=' + str(self.N) + ')$')

                subdivisions = 20
                for j in range(0, subdivisions):
                    k = (self.n_y-1)/subdivisions * j
                    k1 = (self.n_y-1)/subdivisions * (j+1)
                    lx = (self.n_y-1)/subdivisions * h
                    ly = self.y_cl[i][k1] - self.y_cl[i][k]
                    if self.w_cl[i] <= 0.0:
                        px1 = h * k
                        py1 = self.y_cl[i][k]
                        px2 = lx
                        py2 = ly
                    else:
                        px1 = h * k1
                        py1 = self.y_cl[i][k1]
                        px2 = -lx
                        py2 = -ly
                    head_width = ly * min(self.w_cl[i], 0.5) * 5
                    head_length = lx/2 * min(self.w_cl[i], 0.5) * 5
                    overhang = 0.9
                    head_lr = False

                    ax.arrow(px1, py1, px2, py2, head_width=head_width, head_length=head_length, fc='k',
                             ec='b', overhang=overhang, length_includes_head=True, head_starts_at_zero=head_lr)

                ax.hold(True)
                # plot constraints
                ax.plot(subdomain, lb_y * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)
                ax.plot(subdomain, ub_y * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)

                ax.plot([0.95, 1.0], [ub_u, ub_u], color='b', linewidth=linewidth_global)
                ax.plot([0.95, 1.0], [lb_u, lb_u], color='b', linewidth=linewidth_global)

                # plot text with time
                ax.text(0.8, 0.4, 't = {:5.3f}'.format((i)*delta_t), bbox={'facecolor':'white', 'pad':10})
                ax.hold(False)

                if reference:
                    ax.hold(True)
                    ax.plot(domain, reference.y_ol[0][i], color='g', linewidth=linewidth_global, label='$y^*$')
                    ax.hold(False)

                ax.set_xlim([0.0, 1.0])
                ax.set_ylim([-0.5, 0.5])

                ax.set(xlabel='$\Omega$', ylabel='$y$')
                ax.xaxis.label.set_size(20)
                ax.yaxis.label.set_size(20)

                ax.legend(loc='upper left')
                plt.show()
                writer.grab_frame()

    def plot_open_loop(self,output_file,reference=None,k=0):
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='open_loop_N='+str(self.N)+'_k='+str(k), artist='Matplotlib',
                    comment='Movie support!')
        writer = FFMpegWriter(fps=15, metadata=metadata)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        h = 1.0/(self.n_y - 1)
        domain = np.arange(0.0, 1.0+h, h)
        subdomain = np.arange(0.25, 0.75+h, h)
        fig, ax = plt.subplots()

        with writer.saving(fig, output_file, 100):
            for i in range(0, self.N):
                ax.plot(domain, self.y_ol[k][i], color='k', linewidth=linewidth_global, label='$y_{ol} \; (N=' + str(self.N) +')$')

                subdivisions = 20
                for j in range(0, subdivisions):
                    g = (self.n_y - 1) / subdivisions * j
                    g1 = (self.n_y - 1) / subdivisions * (j + 1)
                    lx = (self.n_y - 1) / subdivisions * h
                    ly = self.y_ol[k][i][g1] - self.y_ol[k][i][g]
                    if self.w_ol[k][i] <= 0.0:
                        px1 = h * g
                        py1 = self.y_ol[k][i][g]
                        px2 = lx
                        py2 = ly
                    else:
                        px1 = h * g1
                        py1 = self.y_ol[k][i][g1]
                        px2 = -lx
                        py2 = -ly
                    head_width = ly * min(self.w_ol[k][i], 0.5) * 5
                    head_length = lx / 2 * min(self.w_ol[k][i], 0.5) * 5
                    overhang = 0.9
                    head_lr = False

                    ax.arrow(px1, py1, px2, py2, head_width=head_width, head_length=head_length, fc='k',
                             ec='b', overhang=overhang, length_includes_head=True, head_starts_at_zero=head_lr)

                ax.hold(True)
                # plot constraints
                ax.plot(subdomain, lb_y * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)
                ax.plot(subdomain, ub_y * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)

                ax.plot([0.95, 1.0], [ub_u, ub_u], color='b', linewidth=linewidth_global)
                ax.plot([0.95, 1.0], [lb_u, lb_u], color='b', linewidth=linewidth_global)

                # plot text with time
                ax.text(0.8, 0.4, 't = {:5.3f}'.format((i+1)*delta_t), bbox={'facecolor':'white', 'pad':10})
                ax.hold(False)

                if reference:
                    ax.hold(True)
                    ax.plot(domain, reference.y_ol[k][1+i], color='g', linewidth=linewidth_global, label='$y^*$')
                    ax.hold(False)

                ax.set_xlim([0.0, 1.0])
                ax.set_ylim([-0.5, 0.5])

                ax.set(xlabel='$\Omega$', ylabel='$y$')
                ax.xaxis.label.set_size(20)
                ax.yaxis.label.set_size(20)

                ax.legend(loc='upper left')
                plt.show()
                writer.grab_frame()

    def plot_closed_loop_cost(self):
        domain = np.arange(0.0, self.L)
        fig, ax = plt.subplots()
        ax.hold(False)

        ax.plot(domain, self.closed_loop_cost)
        plt.show()

    # def plot_open_loop(self):
    #     domain = np.arange(0.0, 1.0, 1.0/self.n_y)
    #     fig, ax = plt.subplots()
    #     ax.hold(False)
    #
    #     for i in range(0, self.N):
    #         ax.plot(domain, self.y_ol[i])
    #         ax.set_xlim([0.0, 1.0])
    #         ax.set_ylim([-0.5, 0.5])
    #         plt.show()
    #     pass

def plot_closed_loop_convergence(results, reference, output_file='test.pdf'):
    fig, ax = plt.subplots()
    ax.hold(True)

    for r in results:
        y = r.y_cl
        y_ref = reference.y_ol[0]

        norm_diffs = []
        for j in range(0,r.L):
            norm_diffs.append(np.linalg.norm(y[j] - y_ref[j]))

        domain = np.arange(0, 0+r.L)
        ax.plot(domain, norm_diffs, label='N='+str(r.N), linewidth=linewidth_global)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$k$', fontsize=20)
    plt.ylabel('$\|y_{\mu_{N}}(k,x) - y^*(k)\|$', fontsize=20)
    plt.title('Convergence of MPC trajectories')
    plt.savefig(output_file)
    plt.show()

def plot_cumulative_closed_loop_cost(results, reference, output_file='test.pdf'):
    fig, ax = plt.subplots()
    ax.hold(True)

    for r in results:
        domain = np.arange(0.0, r.L)
        ax.plot(domain, r.closed_loop_cost, label='N='+str(r.N), linewidth=linewidth_global)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$L$', fontsize=20)
    plt.ylabel('$J^{cl}_{L}(y,\mu_N)$', fontsize=20)
    plt.title('Cumulative closed-loop cost')
    plt.savefig(output_file)
    plt.show()

def plot_cost_convergence(results, output_file='test.pdf'):
    L = 0
    Ns = []
    Js = []
    for r in results:
        Ns.append(r.N)
        Js.append(r.final_closed_loop_cost)

        L = max(L, r.L)
    Ns, Js = zip(*sorted(zip(Ns, Js)))

    fig, ax = plt.subplots()
    ax.hold(False)
    ax.plot(Ns, Js, linewidth=linewidth_global)
    ax.set_xlim([np.min(Ns), np.max(Ns)])

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$N$', fontsize=20)
    plt.ylabel('$J^{cl}_{' + str(L) + '}(y,\mu_N)$', fontsize=20)
    plt.title('Convergence of closed-loop cost')
    plt.savefig(output_file)
    plt.show()

def plot_turnpike(results, reference, output_file='test.pdf'):
    fig, ax = plt.subplots()
    ax.hold(True)

    y_ref = reference.y_ol[0]
    for r in results:
        y = r.y_ol[0]

        norm_diffs = []
        for j in range(0,r.N+1):
            norm_diffs.append(np.linalg.norm(y[j] - y_ref[0+j]))

        domain = np.arange(0, 0+r.N+1)
        ax.plot(domain, norm_diffs, label='N=' + str(r.N), linewidth=linewidth_global)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$k$', fontsize=20)
    plt.ylabel('$\|y_{u^*_{N}}(k,x) - y^*(k)\|$', fontsize=20)
    plt.title('Turnpike behaviour')
    plt.savefig(output_file)
    plt.show()

def run_simulations(Ns, L, exec_folder, result_folder, prefix="", ref=False):
    # run simulations
    for N in Ns:
        folder_prefix = prefix + "N=" + str(N) + "_L=" + str(L) + "_"
        if ref:
            call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                  "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                  "--result_folder_prefix=" + folder_prefix, "--fi"])
        else:
            call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                  "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                  "--result_folder_prefix=" + folder_prefix])



if __name__ == "__main__":
    exec_folder = 'cpp/cmake-build-debug/'  # folder with executable
    result_folder = 'results/'              # folder where results are stored

    sim = False
    if sim == True:
        # generate results
        min_N = 40
        max_N = 40
        #Ns = range(min_N,max_N+1)
        Ns = [10, 20, 30, 40]
        L = 75
        run_simulations(Ns, L, exec_folder, result_folder, prefix="mpc_")

        # reference solution (= open loop simulation with long horizon)
        run_simulations([max_N + L], 1, exec_folder, result_folder, prefix="ref_", ref=True)

        # overtaking solution (= open loop simulation with long horizon)
        run_simulations([max_N + L], 1, exec_folder, result_folder, prefix="opt_")

    #run_simulations([1], 250, exec_folder, result_folder, prefix="unc_")

    # handle results
    # MPC simulations
    p = Path(result_folder)
    #result_folder_list = map(str, list(p.glob('mpc_N=5*')) + list(p.glob('mpc_N=10*')) + list(p.glob('mpc_N=20*')) + list(p.glob('mpc_N=30*'))
    #                        + list(p.glob('mpc_N=40*')) + list(p.glob('mpc_N=49*')))    # find all folders with results
    result_folder_list = map(str, list(p.glob('mpc_N=*')))
    #result_folder_list = map(str, list(p.glob('mpc_N=5_*')) + list(p.glob('mpc_N=10*')) + list(p.glob('mpc_N=20_*'))
    #                            + list(p.glob('mpc_N=30_*')) + list(p.glob('mpc_N=40_*')) + list(p.glob('mpc_N=50_*')))
    mpc_list = []
    # save results as objects
    for r in result_folder_list:
        mpc_list.append(SimulationResult(r))
    mpc_list = sorted(mpc_list, key=lambda r: r.N)  # sort by horizon length

    # reference trajectory
    ref_result_folder = map(str,list(p.glob('ref_*')))[0]
    ref_result = SimulationResult(ref_result_folder)

    # overtaking optimal trajectory
    opt_result_folder = map(str, list(p.glob('opt_*')))[0]
    opt_result = SimulationResult(opt_result_folder)

    # uncontrolled
    #unc_folder = map(str, list(p.glob('unc_*')))[0]
    #unc_result = SimulationResult(unc_folder)

    # create turnpike plot
    #plot_turnpike(mpc_list, ref_result, output_file='figures/turnpike.pdf')


    # create plot for convergence of mpc cost
    #mpc_list_ = filter(lambda x: x.N <= 30, mpc_list)
    #plot_cost_convergence(mpc_list_, output_file='figures/cost_convergence.pdf')

    # create plot for convergence of mpc trajectories
    #plot_closed_loop_convergence(mpc_list, ref_result, output_file='figures/trajectory_convergence.pdf')

    # create plot for cumulative closed loop cost
    #plot_cumulative_closed_loop_cost(mpc_list, ref_result, output_file='figures/cumulative_cost.pdf')

    #unc_result.plot_closed_loop('figures/uncontrolled.mp4', L=100)

    # create videos of mpc closed loop simulations
    for r in mpc_list:
        r.plot_closed_loop('figures/heat_N={}.mp4'.format(r.N), L=100, reference=opt_result)
        #r.plot_open_loop('figures/ol_heat_N={}.mp4'.format(r.N), k=0, reference=ref_result)


    pass



    #res1.plot_open_loop()
    #res5.plot_closed_loop(ref_ol_2)
