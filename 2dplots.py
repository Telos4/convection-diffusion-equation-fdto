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

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

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
            self.N = len(self.u_ol[0])      # MPC horizon
            self.n_y = len(self.y_cl[0]) # number of finite element nodes
            self.n_disc = int(np.sqrt(self.n_y))
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


        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #ax.azim = 90
        #ax.elev = 0

        h = 1./(self.n_disc - 1)
        # Make data.
        X = np.arange(0, 1+0.5*h, h)
        Y = np.arange(0, 1+0.5*h, h)
        X, Y = np.meshgrid(X, Y)

        distance_of_arrows =  5
        Xa = np.arange(0, 1+0.5*h, distance_of_arrows * h)
        Ya = np.arange(0, 1+0.5*h, distance_of_arrows * h)
        Xa, Ya = np.meshgrid(Xa, Ya)

        with writer.saving(fig, output_file, 100):
            for i in range(0, L):
                Z = np.reshape(self.y_cl[i], (self.n_disc, self.n_disc))
                Za = np.zeros((len(Xa), len(Ya)))
                for j in range(0, len(Xa)):
                    for k in range(0, len(Ya)):
                        Za[j][k] = Z[distance_of_arrows * j][distance_of_arrows * k]

                Za_left = np.zeros((len(Xa), len(Ya)))
                for j in range(0, len(Xa)):
                    for k in range(0, len(Ya)):
                        if k == 0:
                            Za_left[j][k] = Z[distance_of_arrows * j][distance_of_arrows * k]
                        else:
                            Za_left[j][k] = Z[distance_of_arrows * j][distance_of_arrows * k - 1]

                Za_right = np.zeros((len(Xa), len(Ya)))
                for j in range(0, len(Xa)):
                    for k in range(0, len(Ya)):
                        if distance_of_arrows * k >= len(X) - 1:
                            Za_right[j][k] = Z[distance_of_arrows * j][distance_of_arrows * k]
                        else:
                            Za_right[j][k] = Z[distance_of_arrows * j][distance_of_arrows * k + 1]
                #print(Za)

                #print(Z)
                #print(X)
                # Plot the surface.
                #ax.hold(False)
                surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, vmin = -0.3, vmax = 0.3, linewidth=0.0, antialiased=False, alpha = 0.2)
                #surf = ax.plot_wireframe(X, Y, Z)

                #save difference between Za and Za shifted left/right for z-component of vector, y=0, x depends on sign of w (either +h or -h), multiply x,z by w for length
                factor = 100
                if self.w_cl[i] >= 0:
                    ax.quiver(Xa, Ya, Za, -h * self.w_cl[i] * factor, 0, (Za_left - Za) * self.w_cl[i] * factor)
                else:
                    ax.quiver(Xa, Ya, Za, -h * self.w_cl[i] * factor, 0, (Za_right - Za) * self.w_cl[i] * factor)


                # Customize the z axis.
                ax.set_xlim([0.0, 1.0])
                ax.set_ylim([0.0, 1.0])
                ax.set_zlim([-0.5, 0.5])

                # plot text with time
                ax.text(0.8, 0.4, 1.0, 't = {:5.3f}'.format((i)*delta_t), bbox={'facecolor':'white', 'pad':10})

                ax.set(xlabel='$x$', ylabel='$y$', zlabel = '$value$')
                ax.xaxis.label.set_size(20)
                ax.yaxis.label.set_size(20)
                ax.zaxis.label.set_size(20)

                # Add a color bar which maps values to colors.
                #if(i == 0):
                    #fig.colorbar(surf, shrink=0.5, aspect=5)
                #plt.show()
                writer.grab_frame()
                plt.cla()

    def plot2d(self):



        fig = plt.figure()
        ax = fig.gca(projection='3d')

        h = 1./(self.n_disc - 1)
        # Make data.
        X = np.arange(0, 1+h, h)
        Y = np.arange(0, 1+h, h)
        X, Y = np.meshgrid(X, Y)
        Z = np.reshape(self.y_cl[40], (self.n_disc, self.n_disc))

        #print(Z)
        #print(X)
        # Plot the surface.
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        # Customize the z axis.
        ax.set_zlim(-1.01, 1.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()




def run_simulations(Ns, L, exec_folder, result_folder, prefix="", ref=False):
    # run simulations
    for N in Ns:
        folder_prefix = prefix + "N=" + str(N) + "_L=" + str(L) + "_"
        if ref:
            call([exec_folder + "heateq_opt", "-c", "-d", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                  "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                  "--result_folder_prefix=" + folder_prefix, "--fi", "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=0"])
        else:
            call([exec_folder + "heateq_opt", "-c", "-d", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                  "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                  "--result_folder_prefix=" + folder_prefix, "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=0"])



if __name__ == "__main__":
    exec_folder = 'cpp/'  # folder with executable
    result_folder = 'results3/'              # folder where results are stored

    sim = False
    if sim == True:
        # generate results
        min_N = 40
        max_N = 40
        #Ns = range(min_N,max_N+1)
        Ns = [5]
        L = 100
        run_simulations(Ns, L, exec_folder, result_folder, prefix="mpc_")

        # reference solution (= open loop simulation with long horizon)
        #run_simulations([max_N + L], 1, exec_folder, result_folder, prefix="ref_", ref=True)

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
    #ref_result_folder = map(str,list(p.glob('ref_*')))[0]
    #ref_result = SimulationResult(ref_result_folder)

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
        r.plot_closed_loop('figures/heat_N={}.mp4'.format(r.N), L=100)
        #r.plot_closed_loop('figures/heat_N={}.mp4'.format(r.N), L=100, reference=ref_result)
        #r.plot_open_loop('figures/ol_heat_N={}.mp4'.format(r.N), k=0, reference=ref_result)


    pass



    #res1.plot_open_loop()
    #res5.plot_closed_loop(ref_ol_2)
