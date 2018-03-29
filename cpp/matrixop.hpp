/*
 * hs071_nlp.hpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#ifndef MATRIXOP_HPP_
#define MATRIXOP_HPP_


#include "boost/filesystem.hpp"

#include <valarray>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <ctime>
#include <random>

using namespace std;

class MATRIXOP {
public:

    MATRIXOP();
    MATRIXOP(int N_, int steps, string file_A, string file_B_y, string file_B_w, string file_b_u, string file_b_y,
            string file_dof_x, string file_dof_y, double eps_, double y_ref_, double u_ref_, double u_upper_, double u_lower_, double y_upper_, double y_lower_, double w_upper_, double w_lower_,
            double boundary_right, double boundary_left, double boundary_top, double boundary_bot, double std_dev,
            bool dim2, bool convection, bool closed_values, bool open_values, bool free_init_value, bool inexact_data, string result_folder, string result_folder_name);



    /** Default destructor */
    ~MATRIXOP();

    void read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals);
    void read_vector(string filename, valarray<double> &vals);
    void print_matrix(valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) const;
    void print_vector(valarray<double> &vals) const;
    void create_folder(string folder_prefix);

    double eval_f(valarray<double> &x);
    double vec_Q_vec(valarray<double> y, double y_ref);
    double vec_R_vec(valarray<double> y, double y_ref);
    double vec_W_vec(valarray<double> y);

    valarray<double> eval_grad_f(valarray<double> z);
    valarray<double> Q_vec(valarray<double> y);
    valarray<double> R_vec(valarray<double> u);
    valarray<double> W_vec(valarray<double> w);

    void A_eq();
    void initialize_order();
    valarray<double> matrix_vector_mult(valarray<int> &rows,
            valarray<int> &cols, valarray<double> &vals, valarray<double> &vec, unsigned length);




    /*
     * Save indices of matrix X in X_rows and X_cols and the corresponding value in X_vals.
     * Matrix A:    state variables in the next time step, compare with the weak form of the variational formulation
     * Matrix B_y:  state variables in the current time step, compare with the weak form of the variational formulation
     * Matrix B_W:  convection and state variable in the next time step, compare with the weak form of the variational formulation
     * Matrix A_eq: equation (12) in docu_older.pdf
     * 
     * b_u/b_y:     vectors from the discretization of the weak formulation
     * x_old:       state, control, convection of the previous time step. Used as starting point for the next step
     * 
     * order:       In the 2d case the state variables are not ordered nicely, because of the way Fenics discretizes the variational formulation in heat2d.py. 
     *              order[i] gives the position the i-th state variable should be written to. See also MATRIXOP::initialize_order and HEAT_NLP::finalize_solution
     * dof_x/dof_y: 2d only, gives the x and y indices of every state variable which are then used to initialize the order array
     * error:       length steps+MPC_horizon, save deviation for every time step
     * 
     * n_y:         number of state variables
     * n_u:         number of control variables, always 1 as of now
     * n_w:         number of convection variables, always 1 as of now (if convection is enabled)
     * n_z:         total number of variables in the optimization, calculated by n_z = (N + 1) * n_y + N * n_u + N * n_w + N * n_y
     * N:           MPC_horizon 
     * discretization_n:    discretization parameter in heat.py or heat2d.py. This means we have a (n x 1) or (n x n) grid and accordingly n+1 or (n+1)^2 degrees of freedom on our grid
     * steps:       number of closed loop steps
     * 
     * eps:         epsilon of the cost functional
     * y_ref:       reference solution for y from cost functional
     * u_ref:       reference solution for u from cost functional
     * closed_loop_cost:    gives the total closed loop cost after the closed loop has finished, update in HEAT_NLP::finalize_solution
     * x_upper/lower:       upper/lower boundaries for state y, control u, convection w
     * boundary_*:  specifies the edges of the rectangle where the state constraints are applied in the two dimensional case
     * std_dev:    set the standard deviation in the case of inexact_data
     * 
     * dim2:            true: two dimensional false: one dimensional
     * convection:      convection on or off
     * closed_values:   write the values of the closed loop optimization to a file (just the actual trajectory and not all the open loop values) 
     * open_values:     create a file for every optimization step and write all parameters of this step
     * free_init_value: the first optimization step has a free initial value and doesn't start at the reference solution
     * inexact_data:    indicates whether the boundary data y_out should be exact or have a gaussian distributed error with mean 0 and standard deviation set by std_dev
     * 
     * foldername:      folder where results shall be stored
     * result_folder:   folder where the above folder shall be placed
     */

    valarray<int> A_rows, A_cols, B_y_rows, B_y_cols, A_eq_rows, A_eq_cols, B_w_rows, B_w_cols, order;
    valarray<double> A_vals, B_y_vals, A_eq_vals, B_w_vals, b_u, b_y, y_old, u_old, w_old, dof_x, dof_y, error;
    int n_y, n_u, n_w, n_z, N, iter, discretization_n, steps;
    double eps, y_ref, u_ref, closed_loop_cost, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, boundary_right, boundary_left, boundary_top, boundary_bot, std_dev;

    bool dim2, convection, closed_values, open_values, free_init_value, inexact_data;
    string foldername;
    string result_folder;
private:


};





#endif /* MATRIXOP_HPP_ */
