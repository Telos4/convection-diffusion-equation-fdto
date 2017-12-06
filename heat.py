from __future__ import print_function
from fenics import *
import numpy as np
import os, sys
import scipy.io

from collections import OrderedDict

# Set log level
set_log_level(WARNING)

# Prepare a mesh
mesh = UnitIntervalMesh(100)

# Choose a time step size
k = Constant(1e-2)

# MPC horizon length
N = 10

# boundary heat conductivity parameters
alpha = Constant(1.0)
beta = Constant(1.0)
gamma = Constant(1.0e3)

# Compile sub domains for boundaries
left = CompiledSubDomain("near(x[0], 0.)")
right = CompiledSubDomain("near(x[0], 1.)")

# Label boundaries, required for the objective
boundary_parts = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
left.mark(boundary_parts, 0)    # boundary part where control is applied
right.mark(boundary_parts, 1)   # boundary part for outside temperature
ds = Measure("ds", subdomain_data=boundary_parts)

def output_matrices():
    # Define function space
    U = FunctionSpace(mesh, "Lagrange", 1)

    # Define test and trial functions
    v = TestFunction(U)
    y = TrialFunction(U)
    y0 = TrialFunction(U)

    u = Constant(1.0)
    y_out = Constant(1.0)

    # Define variational formulation
    a = (y / k * v + alpha * inner(grad(y), grad(v))) * dx + alpha * gamma/beta * y * v * ds
    f_y = y0 / k * v * dx

    f_u = alpha * gamma/beta * u * v * ds(0)

    f_y_out = alpha * gamma/beta * y_out * v * ds(1)

    A = assemble(a)
    B_y = assemble(f_y)

    b_u = assemble(f_u)
    b_y_out = assemble(f_y_out)

    # output matrices for use in matlab optimization
    scipy.io.savemat('sys.mat', {'A': A.array(), 'B_y': B_y.array(), 'b_u': b_u.array(), 'b_y_out': b_y_out.array(),
                                 'N': N, 'u': us, 'y_out': y_outs})

    return b_u, b_y_out


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


if __name__ == "__main__":
    us = np.array([0.1 * (i+1) for i in range(0,N)])
    y_outs = np.array([1.0 - 0.1 * (i+1) for i in range(0,N)])

    output_matrices()

    solve_forward(us, y_outs)
