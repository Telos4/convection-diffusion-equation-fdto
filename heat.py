from __future__ import print_function
from fenics import *
import numpy as np
import os, sys
import scipy.io

from collections import OrderedDict

# Set log level
set_log_level(WARNING)

# Prepare a mesh
#mesh = UnitSquareMesh(10,10)
mesh = UnitIntervalMesh(50)

# Choose a time step size
k = Constant(1e-1)

# MPC horizon length
N = 10

# boundary heat conductivity parameters
alpha = Constant(1.0e3)
beta = Constant(1.0e3)

# Compile sub domains for boundaries
left = CompiledSubDomain("near(x[0], 0.)")
right = CompiledSubDomain("near(x[0], 1.)")

# Label boundaries, required for the objective
boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
left.mark(boundary_parts, 0)    # boundary part where control is applied
right.mark(boundary_parts, 1)   # boundary part for outside temperature
ds = Measure("ds", subdomain_data=boundary_parts)

def output_matrices(us, y_outs):
    # Define function space
    U = FunctionSpace(mesh, "Lagrange", 1)

    # Define test and trial functions
    v = TestFunction(U)
    y = TrialFunction(U)
    y0 = TrialFunction(U)

    u = Constant(1.0)
    y_out = Constant(1.0)

    # Define variational formulation
    a = (y / k * v + inner(grad(y), grad(v))) * dx + alpha * y * v * ds
    f_y = y0 / k * v * dx

    f_u = beta * u * v * ds(0)

    f_y_out = beta * y_out * v * ds(1)

    A = assemble(a)
    B_y = assemble(f_y)

    b_u = assemble(f_u)
    b_y_out = assemble(f_y_out)

    scipy.io.savemat('sys.mat', {'A': A.array(), 'B_y': B_y.array(), 'b_u': b_u.array(), 'b_y_out': b_y_out.array(),
                                 'N': N, 'u': us, 'y_out': y_outs})


def solve_forward(us, y_outs, record=False):
    """ The forward problem """
    ofile = File("results/y.pvd")

    # Define function space
    U = FunctionSpace(mesh, "Lagrange", 1)

    # Set up initial values
    y0 = Function(U, name="y0")

    # Define test and trial functions
    v = TestFunction(U)
    y = TrialFunction(U)
    u = Constant(1.0)
    y_out = Constant(1.0)

    # Define variational formulation
    #F = ((y - y0)/k * v + inner(grad(y), grad(v))) * dx - (1.0e3 * (u - y) * v) * ds
    #a = lhs(F)
    #L = rhs(F)
    a = (y/k * v + inner(grad(y), grad(v))) * dx + alpha * y * v * ds
    f_y = y0 / k * v * dx

    f_u = beta * u * v * ds(0)

    f_y_out = beta * y_out * v * ds(1)

    # Prepare solution
    y = Function(U, name="y")

    i = 0

    ys = OrderedDict()
    y_omegas = OrderedDict()
    #ys[i] = y0.copy(deepcopy=True)
    y_omegas[i] = Function(U, name="y_omega[0]")

    As = []
    bs = []

    while i < N:
        #plot(y0)
        t = i*float(k)
        #ofile << (y, t)
        u.assign(us[i])
        y_out.assign(y_outs[i])

        solve(a == f_u + f_y + f_y_out, y)
        y0.assign(y)

        i += 1

    scipy.io.savemat('y.mat', {'y': y.vector().array()})


    return y, ys, y_omegas


# Callback function for the optimizer
# Writes intermediate results to a logfile
def eval_cb(j, ms):
    """ The callback function keeping a log """

    for m in ms:
        print("m = %15.10e " % float(m))
    print("objective = %15.10e " % j)


# Prepare the objective function
def objective(times, y):
    """ The objective """

    y_ref = Constant(0.5)
    #combined = zip(times, observations)
    #area = times[-1] - times[0]
    #M = len(times)
    #I = area / M * sum(inner(y - y_obs, y - y_obs) * ds(1) * dt[t]
    #                   for (t, y_obs) in combined)
    I = sum(inner(y-y_ref,y-y_ref) * dx * dt[t]
            for t in times)

    return I


def optimize(dbg=False):
    """ The optimization routine """

    # Define the control
    #us = [Constant(float(i)) for i in range(0,10)]
    #us = [0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    us = OrderedDict()
    i = 0
    while i < N:
        us[i] = Constant(1.0)
        i += 1

    #source = Source(us, degree=3, name="source")

    # provide the coefficient on which this expression depends and its derivative
    #source.dependencies = us
    #source.user_defined_derivatives = {}
    #for u in us:
    #    source.user_defined_derivatives[u] = Source(us, Source=source, derivative=us, degree=3)
    #source.user_defined_derivatives = {us: Source(us, Source=source, derivative=us, degree=3)}

    # Execute first time to annotate and record the tape
    y,  ys, y_omegas = solve_forward(us, record=False, annotate=False)

    p, ps = solve_adjoint(ys, y_omegas)

    if dbg:
        # Check the recorded tape
        success = replay_dolfin(tol=0.0, stop=True)
        print("replay: ", success)

        # for the equations recorded on the forward run
        adj_html("forward.html", "forward")
        # for the equations to be assembled on the adjoint run
        adj_html("adjoint.html", "adjoint")

    # Load references
    #refs = np.loadtxt("recorded.txt")

    # create noise to references
    # gamma = 1.e-5
    # if gamma > 0:
    #     noise = np.random.normal(0, gamma, refs.shape[0])
    #
    #     # add noise to the refs
    #     refs += noise
    #
    # # map refs to be constant
    # refs = list(map(Constant, refs))

    # Define the control
    controls = [Control(u) for (t,u) in us.items()]

    Jform = objective(times, y)
    J = Functional(Jform)

    # compute the gradient
    dJd0 = compute_gradient(J, controls)
    #for elem in dJd0:
    #    print("gradient = ", float(elem))

    # Prepare the reduced functional
    reduced_functional = ReducedFunctional(J, controls, eval_cb_post=eval_cb)

    # Run the optimisation
    omega_opt = minimize(reduced_functional, method="L-BFGS-B", \
                         tol=1.0e-12, options={"disp": True, "gtol": 1.0e-12})

    # Print the obtained optimal value for the controls
    print("omega = %f" % float(omega_opt))


if __name__ == "__main__":
    us = np.array([0.1 * (i+1) for i in range(0,N)])
    y_outs = np.array([1.0 - 0.1 * (i+1) for i in range(0,N)])

    output_matrices(us, y_outs)

    solve_forward(us, y_outs)
