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

class SimulationResult:
    def __init__(self, results_folder):
        p = Path(results_folder)

        q = p / 'closedloop_y.txt'
        if q.exists():
            self.y_cl = np.loadtxt(str(q), ndmin=2)

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
                self.y_ol.append(np.loadtxt(str(q), ndmin=1))

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
            self.N = len(self.u_ol[0])
            self.n_y = len(self.y_ol[0][0])
        except AttributeError:
            print("open loop not available")



    def plot_closed_loop(self, reference=None):
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='closed_loop_L='+str(self.L), artist='Matplotlib',
                    comment='Movie support!')
        writer = FFMpegWriter(fps=15, metadata=metadata)

        domain = np.arange(0.0, 1.0, 1.0/self.n_y)
        fig, ax = plt.subplots()
        ax.hold(False)

        with writer.saving(fig, "results/sim_N=" + str(self.N) + ".mp4", 100):
            for i in range(0, self.L):
                ax.plot(domain, self.y_cl[i])

                if reference:
                    ax.hold(True)
                    ax.plot(domain, reference.y_ol[0][1+i])
                    ax.hold(False)

                ax.set_xlim([0.0, 1.0])
                ax.set_ylim([-0.5, 0.5])
                plt.show()
                writer.grab_frame()

    def plot_closed_loop_cost(self):
        domain = np.arange(0.0, self.L)
        fig, ax = plt.subplots()
        ax.hold(False)

        ax.plot(domain, self.closed_loop_cost)
        plt.show()

    def plot_open_loop(self):
        domain = np.arange(0.0, 1.0, 1.0/self.n_y)
        fig, ax = plt.subplots()
        ax.hold(False)

        for i in range(0, self.N):
            ax.plot(domain, self.y_ol[i])
            ax.set_xlim([0.0, 1.0])
            ax.set_ylim([-0.5, 0.5])
            plt.show()
        pass

def movie_example():
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib',
            comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    fig = plt.figure()
    l, = plt.plot([], [], 'k-o')

    plt.xlim(-5, 5)
    plt.ylim(-5, 5)

    x0, y0 = 0, 0

    with writer.saving(fig, "writer_test.mp4", 100):
        for i in range(100):
            x0 += 0.1 * np.random.randn()
            y0 += 0.1 * np.random.randn()
        l.set_data(x0, y0)
        writer.grab_frame()

def plot_closed_loop_convergence(results, reference, output_file='test.pdf'):
    fig, ax = plt.subplots()
    ax.hold(True)

    for r in results:
        y = r.y_cl
        y_ref = reference.y_ol[0]

        norm_diffs = []
        for j in range(0,r.L):
            norm_diffs.append(np.linalg.norm(y[j] - y_ref[1+j]))

        domain = np.arange(0, 0+r.L)
        ax.plot(domain, norm_diffs, label='N='+str(r.N))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$k$')
    plt.ylabel('$\|x_{\mu_{N}}(k,x) - x^*(k)\|$')
    plt.title('Convergence of MPC trajectories')
    plt.savefig(output_file)
    plt.show()

def plot_cumulative_closed_loop_cost(results, reference, output_file='test.pdf'):
    fig, ax = plt.subplots()
    ax.hold(True)

    for r in results:
        domain = np.arange(0.0, r.L)
        ax.plot(domain, r.closed_loop_cost, label='N='+str(r.N))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$L$')
    plt.ylabel('$J^{cl}_{L}(x,\mu_N)$')
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
    ax.plot(Ns, Js)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$N$')
    plt.ylabel('$J^{cl}_{' + str(L) + '}(x,\mu_N)$')
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
        ax.plot(domain, norm_diffs, label='N=' + str(r.N))

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.legend(loc='best')
    plt.xlabel('$k$')
    plt.ylabel('$\|x_{u^*_{N}}(k,x) - x^*(k)\|$')
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
    result_folder = 'test_results/'              # folder where results are stored

    sim = False
    if sim == True:
        # generate results
        max_N = 5
        L = 50
        run_simulations(range(1,max_N), L, exec_folder, result_folder, prefix="mpc_")

        # reference solution (= open loop simulation with long horizon)
        run_simulations([max_N + L], 1, exec_folder, result_folder, prefix="ref_", ref=False)

    # handle results
    # MPC simulations
    p = Path(result_folder)
    result_folder_list = map(str, list(p.glob('mpc_N=5*')) + list(p.glob('mpc_N=10*')) + list(p.glob('mpc_N=20*')) + list(p.glob('mpc_N=30*'))
                            + list(p.glob('mpc_N=40*')) + list(p.glob('mpc_N=49*')))    # find all folders with results
    mpc_list = []
    # save results as objects
    for r in result_folder_list:
        mpc_list.append(SimulationResult(r))

    # reference trajectory
    ref_result_folder = map(str,list(p.glob('ref_*')))[0]
    ref_result = SimulationResult(ref_result_folder)

    # create turnpike plot
    plot_turnpike(mpc_list, ref_result, output_file='figures/turnpike.pdf')

    # create plot for convergence of mpc cost
    plot_cost_convergence(mpc_list, output_file='figures/cost_convergence.pdf')

    # create plot for convergence of mpc trajectories
    plot_closed_loop_convergence(mpc_list, ref_result, output_file='figures/trajectory_convergence.pdf')

    # create plot for cumulative closed loop cost
    plot_cumulative_closed_loop_cost(mpc_list, ref_result, output_file='figures/cumulative_cost.pdf')

    # create videos of mpc closed loop simulations
    #for r in mpc_list:
    #    r.plot_closed_loop(ref_result)

    pass



    #res1.plot_open_loop()
    #res5.plot_closed_loop(ref_ol_2)