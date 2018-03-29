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
import ConfigParser as cp
import sys
import mpl_toolkits.mplot3d.art3d as art3d

from subprocess import call
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.patches import Rectangle

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
        config = cp.ConfigParser()
        q = p / 'parameters.txt'
        if q.exists():
            config.read(str(q))

            self.n_disc = config.getint('param', 'discretization_parameter')
            self.L = config.getint('param', 'steps')
            self.N = config.getint('param', 'MPC_horizon')

            self.alpha = config.getfloat('param', 'alpha')
            self.beta = config.getfloat('param', 'beta')
            self.y_ref = config.getfloat('param', 'y_ref')
            self.gamma = config.getfloat('param', 'gamma')
            self.u_ref = config.getfloat('param', 'u_ref')
            self.u_upper = config.getfloat('param', 'u_upper')
            self.u_lower = config.getfloat('param', 'u_lower')
            self.y_upper = config.getfloat('param', 'y_upper')
            self.y_lower = config.getfloat('param', 'y_lower')
            self.w_upper = config.getfloat('param', 'w_upper')
            self.w_lower = config.getfloat('param', 'w_lower')
            self.boundary_right = config.getfloat('param', 'boundary_right')
            self.boundary_left = config.getfloat('param', 'boundary_left')
            self.boundary_top = config.getfloat('param', 'boundary_top')
            self.boundary_bot = config.getfloat('param', 'boundary_bot')
            self.std_dev = config.getfloat('param', 'std_dev')

            self.dim2 = config.getboolean('param', 'dim2')
            self.convection = config.getboolean('param', 'convection')
            self.inexact_data = config.getboolean('param', 'inexact_data')
            self.function = config.get('param', 'function')
        else:
            print("can't find parameters.txt")
            sys.exit()


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

        #2 dimensions
        if self.dim2:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            #ax.azim = 90
            #ax.elev = 0

            h = 1./(self.n_disc)
            # Make data.
            X = np.arange(0, 1+0.5*h, h)
            Y = np.arange(0, 1+0.5*h, h)
            X, Y = np.meshgrid(X, Y)

            distance_of_arrows = 2
            Xa = np.arange(0, 1+0.5*h, distance_of_arrows * h)
            Ya = np.arange(0, 1+0.5*h, distance_of_arrows * h)
            Xa, Ya = np.meshgrid(Xa, Ya)

            with writer.saving(fig, output_file, 100):
                for i in range(0, L):
                    Z = np.reshape(self.y_cl[i], (self.n_disc + 1, self.n_disc + 1))
                    Z = np.fliplr(Z)
                    # Plot the surface.
                    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, vmin = -0.3, vmax = 0.3, linewidth=0.0, antialiased=False, alpha = 0.3)


                    #visualization of convection
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

                    #save difference between Za and Za shifted left/right for z-component of vector, y=0, x depends on sign of w (either +h or -h)
                    #multiply x, z by factor, so all arrows have length 100*w_cl[i]

                    factor = 0.1
                    if self.w_cl[i] >= 0:
                        ax.quiver(Xa, Ya, Za, -h * self.w_cl[i] * factor, 0, (Za_left - Za) * self.w_cl[i] * factor)
                    else:
                        ax.quiver(Xa, Ya, Za, -h * self.w_cl[i] * factor, 0, (Za_right - Za) * self.w_cl[i] * factor)

                    #upper boundary y
                    p1 = Rectangle((self.boundary_bot, self.boundary_left), self.boundary_right - self.boundary_left, self.boundary_top - self. boundary_bot,
                        alpha = 0.9, color = 'green')
                    ax.add_patch(p1)
                    art3d.pathpatch_2d_to_3d(p1, z = self.y_upper, zdir="z")

                    #lower boundary y
                    p2 = Rectangle((self.boundary_bot, self.boundary_left), self.boundary_right - self.boundary_left, self.boundary_top - self. boundary_bot,
                        alpha = 0.9, color = 'green')
                    ax.add_patch(p2)
                    art3d.pathpatch_2d_to_3d(p2, z = self.y_lower, zdir="z")

                    #upper boundary u
                    p3 = Rectangle((0.98, 0), 0.02, 1, alpha = 0.9, color = 'green')
                    ax.add_patch(p3)
                    art3d.pathpatch_2d_to_3d(p3, z = self.u_upper, zdir="z")

                    #lower boundary u
                    p3 = Rectangle((0.98, 0), 0.02, 1, alpha = 0.9, color = 'green')
                    ax.add_patch(p3)
                    art3d.pathpatch_2d_to_3d(p3, z = self.u_lower, zdir="z")

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
                    if(i == 0):
                        fig.colorbar(surf, shrink=0.5, aspect=5)

                    #plt.show()

                    writer.grab_frame()
                    plt.cla()


        #1 dimension
        else:
            h = 1.0/(self.n_disc)
            domain = np.arange(0.0, 1.0+0.5*h, h)
            subdomain = np.arange(0.25, 0.75+0.5*h, h)
            fig, ax = plt.subplots()
            #ax.hold(False)

            with writer.saving(fig, output_file, 100):
                for i in range(0, L):
                    ax.plot(domain, self.y_cl[i], color='k', linewidth=linewidth_global, label='$y\;(N=' + str(self.N) + ')$')

                    subdivisions = 20
                    for j in range(0, subdivisions):
                        k = (self.n_disc-1)/subdivisions * j
                        k1 = (self.n_disc-1)/subdivisions * (j+1)
                        lx = (self.n_disc-1)/subdivisions * h
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
                    ax.plot(subdomain, self.y_lower * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)
                    ax.plot(subdomain, self.y_upper * np.ones(subdomain.shape), color='r', linewidth=linewidth_global)

                    ax.plot([0.95, 1.0], [self.u_upper, self.u_upper], color='b', linewidth=linewidth_global)
                    ax.plot([0.95, 1.0], [self.u_lower, self.u_lower], color='b', linewidth=linewidth_global)

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
                    #plt.show()
                    writer.grab_frame()

    def plot_open_loop(self,output_file,reference=None,k=0):
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='open_loop_N='+str(self.N)+'_k='+str(k), artist='Matplotlib',
                    comment='Movie support!')
        writer = FFMpegWriter(fps=15, metadata=metadata)

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        h = 1.0/(self.n_disc - 1)
        domain = np.arange(0.0, 1.0+0.5*h, h)
        subdomain = np.arange(0.25, 0.75+0.5*h, h)
        fig, ax = plt.subplots()

        with writer.saving(fig, output_file, 100):
            for i in range(0, self.N):
                ax.plot(domain, self.y_ol[k][i], color='k', linewidth=linewidth_global, label='$y_{ol} \; (N=' + str(self.N) +')$')

                subdivisions = 20
                for j in range(0, subdivisions):
                    g = (self.n_disc - 1) / subdivisions * j
                    g1 = (self.n_disc - 1) / subdivisions * (j + 1)
                    lx = (self.n_disc - 1) / subdivisions * h
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
    #     domain = np.arange(0.0, 1.0, 1.0/self.n_disc)
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

    L = results[0].L
    domain = np.arange(0.0, L)
    ax.plot(domain, reference.closed_loop_cost[:L], label='N='+str(reference.N + reference.L), linewidth=linewidth_global)

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

def run_simulations(Ns, L, dim, exec_folder, result_folder, prefix="", ref=False, inexact_data = False, std_dev = 0):
    # run simulations

    for N in Ns:
        folder_prefix = prefix + "N=" + str(N) + "_L=" + str(L) + "_"
        if dim == 1:
            if inexact_data:
                if ref:
                    call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                          "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                          "--result_folder_prefix=" + folder_prefix, "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=0", "--fi", "--id", "--std_dev=" + str(std_dev)])
                else:
                    call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                          "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder, "--output=0",
                          "--result_folder_prefix=" + folder_prefix,
                          "--y_lower=-0.15", "--y_upper=0.15",
                          "--pythonparam=python_parameters.txt", "--id", "--std_dev=" + str(std_dev)])
            else:
                if ref:
                    call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                          "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                          "--result_folder_prefix=" + folder_prefix, "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=0", "--fi"])
                else:
                    call([exec_folder + "heat", "-c", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                          "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder, "--output=0",
                          "--result_folder_prefix=" + folder_prefix,
                          "--y_lower=-0.15", "--y_upper=0.15",
                          "--pythonparam=python_parameters.txt"])


        elif dim == 2:
            if ref:
                call([exec_folder + "heat", "-c", "-d", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                      "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                      "--result_folder_prefix=" + folder_prefix, "--fi", "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=5"])
            else:
                call([exec_folder + "heat", "-c", "-d", "-L " + str(L), "-N" + str(N), "--ov", "--cv", "--matA=A.mtx", "--matB_w=B_w.mtx",
                      "--matB_y=B_y.mtx", "--b_u=b_u.txt", "--b_y_out=b_y_out.txt", "--result_folder=" + result_folder,
                      "--result_folder_prefix=" + folder_prefix, "--pythonparam=python_parameters.txt", "--dof_x=dof_x.txt", "--dof_y=dof_y.txt", "--output=0"])


if __name__ == "__main__":
    exec_folder = 'cpp/'#cmake-build-debug/'  # folder with executable
    result_folder = 'results99/'              # folder where results are stored

    sim = True
    #sim = False
    dim = 1
    inexact_data = True
    std_dev = 0.02
    if sim == True:
        # generate results
        min_N = 40
        max_N = 20
        #Ns = range(min_N,max_N+1)
        Ns = [2,5,10,20]
        L = 100
        run_simulations(Ns, L, dim, exec_folder, result_folder, prefix="mpc_", inexact_data = inexact_data, std_dev = std_dev)

        # reference solution (= open loop simulation with long horizon)
        run_simulations([max_N + L], 1, dim, exec_folder, result_folder, prefix="ref_", ref=True, inexact_data = inexact_data, std_dev = std_dev)

        # overtaking solution (= open loop simulation with long horizon)
        #run_simulations([max_N + L], 1, dim, exec_folder, result_folder, prefix="opt_")

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
    #opt_result_folder = map(str, list(p.glob('opt_*')))[0]
    #opt_result = SimulationResult(opt_result_folder)

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
    plot_cumulative_closed_loop_cost(mpc_list, ref_result, output_file='figures2/cumulative_cost.pdf')

    #unc_result.plot_closed_loop('figures/uncontrolled.mp4', L=100)

    # create videos of mpc closed loop simulations
    for r in mpc_list:
        r.plot_closed_loop('figures2/heat_N={}_2_stddev={}.mp4'.format(r.N, r.std_dev), L=100, reference=ref_result)
        #r.plot_open_loop('figures/ol_heat_N={}.mp4'.format(r.N), k=0, reference=ref_result)

    pass



    #res1.plot_open_loop()
    #res5.plot_closed_loop(ref_ol_2)
