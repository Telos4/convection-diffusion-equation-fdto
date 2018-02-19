/*
 * HEAT_NLP.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include "heat_nlp.hpp"

#include <cassert>

using namespace Ipopt;

double eval_y_out(int i);

// constructor

//HEAT_NLP::HEAT_NLP() {
//}

//destructor

HEAT_NLP::~HEAT_NLP() {
}

HEAT_NLP::HEAT_NLP(MATRIXOP &data_)
: data(data_) {
}

bool HEAT_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	Index& nnz_h_lag, IndexStyleEnum& index_style) {

    n = data.n_z;
    m = data.N * data.n_y;

    if (data.convection) {
	//need to compare A and B_w, always the same size?
	nnz_jac_g = (data.B_y_rows.size() + data.A_rows.size() + data.b_u.size() + data.n_y) * data.N;
    }
    else {
	//jac(g) = A_eq
	nnz_jac_g = (data.B_y_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;
    }

    //mit Q = I, R = I for no convection
    //don't need this when using approximation of hessian of the lagrangian with convection on
    nnz_h_lag = data.n_z;


    index_style = TNLP::C_STYLE;

    return true;
}

bool HEAT_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
	Index m, Number* g_l, Number* g_u) {


    double inf = 1e19;

    //initial value
    if (data.free_init_value && data.iter == 0) {
	for (int i = 0; i < data.n_y; ++i) {
	    x_u[i] = inf; //data.y_old[data.n_y + i];
	    x_l[i] = -inf; //data.y_old[data.n_y + i];
	}
    }
    else {
	for (int i = 0; i < data.n_y; ++i) {
	    x_u[i] = data.y_old[data.n_y + i];
	    x_l[i] = data.y_old[data.n_y + i];
	}
    }

    //state constraints in 2d
    if (data.dim2) {

	double left = 0.3;
	double right = 0.7;
	double top = 0.7;
	double bot = 0.3;


	for (int k = 1; k < data.N + 1; ++k) {
	    for (int i = 0; i < data.n_y; ++i) {
		double x = data.dof_x[i];
		double y = data.dof_y[i];
		if (x >= left && x <= right && y >= bot && y <= top) {
		    x_u[k * data.n_y + i] = data.y_upper;
		    x_l[k * data.n_y + i] = data.y_lower;
		}
		else {
		    x_u[k * data.n_y + i] = inf;
		    x_l[k * data.n_y + i] = -inf;
		}
	    }
	}
    }

	//state constraints in 1d
    else {

	int left = floor(data.n_y / 4);
	//int right = (int) (data.n_y * 3 / 4);

	for (int k = 1; k < data.N + 1; ++k) {
	    //bound for y left, right
	    for (int i = 0; i < left; ++i) {
		x_u[k * data.n_y + i] = x_u[(k + 1) * data.n_y - i - 1] = inf;
		x_l[k * data.n_y + i] = x_l[(k + 1) * data.n_y - i - 1] = -inf;
	    }

	    //bound for y center
	    for (int i = left; i < data.n_y - left + 1; ++i) {
		x_u[k * data.n_y + i] = data.y_upper;
		x_l[k * data.n_y + i] = data.y_lower;
	    }
	    /*
		    //bound for y right
		    for (Index i = right; i < data.n_y; ++i) {
			x_u[k * data.n_y + i] = inf;
			x_l[k * data.n_y + i] = -inf;
		    }
	     */
	}
    }

    //bound for u
    for (int i = (data.N + 1) * data.n_y; i < (data.N + 1) * data.n_y + data.N * data.n_u; ++i) {
	x_u[i] = data.u_upper;
	x_l[i] = data.u_lower;
    }

    if (data.convection) {
	//bound for w
	for (int i = (data.N + 1) * data.n_y + data.N * data.n_u; i < data.n_z; ++i) {
	    x_u[i] = data.w_upper;
	    x_l[i] = data.w_lower;
	}
    }

    //equality constraints for g
    for (int k = 0; k < data.N; ++k) {
	for (int i = 0; i < data.n_y; ++i) {
	    g_u[data.n_y * k + i] = g_l[data.n_y * k + i] = data.b_y[i] * eval_y_out(data.iter + k);
	}
    }

    return true;
}

//where to put this function?

double eval_y_out(int i) {
    return 0.3 * sin(0.1 * i);
}
//Lösung vom letzten

bool HEAT_NLP::get_starting_point(Index n, bool init_x, Number* x,
	bool init_z, Number* z_L, Number* z_U,
	Index m, bool init_lambda,
	Number * lambda) {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);


    //use the shifted solution from the previous solution,
    //for the uninitialized values use the last step again
    for (int i = 0; i < data.N * data.n_y; ++i) {
	x[i] = data.y_old[data.n_y + i];
    }
    for (int i = data.N * data.n_y; i < (data.N + 1) * data.n_y; ++i) {
	x[i] = data.y_old[i];
    }

    for (int i = 0; i < (data.N - 1) * data.n_u; ++i) {
	x[(data.N + 1) * data.n_y + i] = data.u_old[i + data.n_u];
    }
    for (int i = (data.N - 1) * data.n_u; i < data.N * data.n_u; ++i) {
	x[(data.N + 1) * data.n_y + i] = data.u_old[i];
    }

    if (data.convection) {
	for (int i = 0; i < (data.N - data.n_w) * data.n_w; ++i) {
	    x[(data.N + 1) * data.n_y + data.N * data.n_u + i] = data.w_old[i + data.n_w];
	}
	for (int i = (data.N - data.n_w) * data.n_w; i < data.N * data.n_w; ++i) {
	    x[(data.N + 1) * data.n_y + data.N * data.n_u + i] = data.w_old[i];
	}
    }

    return true;
}

bool HEAT_NLP::eval_f(Index n, const Number* x, bool new_x, Number & obj_value) {

    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }

    obj_value = data.eval_f(z);
    return true;
}

//valarray Konstruktor mit Pointer

bool HEAT_NLP::eval_grad_f(Index n, const Number* x, bool new_x,
	Number * grad_f) {

    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }
    valarray<double> grad = data.eval_grad_f(z);
    for (int i = 0; i < n; ++i) {
	grad_f[i] = grad[i];
    }
    return true;
}

bool HEAT_NLP::eval_g(Index n, const Number* x, bool new_x,
	Index m, Number * g) {


    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }
    //this works with convection on as well. You're then looking at the matrix [A_eq, 0].
    //Because of the 0 matrix the matrix vector multiplication works, despite the matrix
    //having the wrong dimension. Add the non linear part afterwards
    valarray<double> eval1 = data.matrix_vector_mult(data.A_eq_rows, data.A_eq_cols,
	    data.A_eq_vals, z);

    valarray<double> eval2(0.0, m);
    //nonlinear part
    if (data.convection) {
	//valarray<double> eval2(m);
	for (int i = 0; i < data.N; ++i) {
	    valarray<double> temp = z[slice((i + 1) * data.n_y, data.n_y, 1)];
	    eval2[slice(i * data.n_y, data.n_y, 1)] = z[(data.N + 1) * data.n_y + data.N * data.n_u + i]
		    * data.matrix_vector_mult(data.B_w_rows, data.B_w_cols, data.B_w_vals, temp);
	}
    }

    //eval2 = 0 if convection is off
    for (int i = 0; i < data.N * data.n_y; ++i) {
	g[i] = eval1[i] + eval2[i];
    }
    return true;
}


//copy routine from A_eq and change
//works only if A and B_w have the same pattern.
//enough for testing for now
//implement general case later
//as long as A and B_w should always have the same pattern, because
//need to change nnz_jacobian in get_nlp_info as well

bool HEAT_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
	Index m, Index nele_jac, Index* iRow,
	Index *jCol, Number * values) {

    if (data.convection) {
	long size_B_y = data.B_y_rows.size();
	//size A = size B_w for now
	long size_A = data.A_rows.size();

	assert(size_A == data.B_w_rows.size());

	//set pattern
	if (values == NULL) {
	    int count = 0;
	    for (int i = 0; i < data.N; ++i) {
		int countA = 0;
		int countB_y = 0;

		for (int k = 0; k < data.n_y; ++k) {
		    while (countB_y != size_B_y && data.B_y_rows[countB_y] == k) {
			iRow[count] = i * data.n_y + k;
			jCol[count] = i * data.n_y + data.B_y_cols[countB_y];
			++count;
			++countB_y;

		    }
		    //change this part
		    while (countA != size_A && data.A_rows[countA] == k) {
			iRow[count] = i * data.n_y + k;
			jCol[count] = (i + 1) * data.n_y + data.A_cols[countA];
			++count;
			++countA;
		    }

		    iRow[count] = i * data.n_y + k;
		    jCol[count] = (data.N + 1) * data.n_y + i * data.n_u;
		    ++count;

		    iRow[count] = i * data.n_y + k;
		    jCol[count] = (data.N + 1) * data.n_y + data.N * data.n_u + i * data.n_w;
		    ++count;
		}
	    }
	}
	    //set values
	else {
	    int count = 0;
	    for (int i = 0; i < data.N; ++i) {
		int countA = 0;
		int countB_y = 0;

		//w_i
		double w = x[(data.N + 1) * data.n_y + data.N * data.n_u + i];

		//B_w . y_i
		valarray<double> temp(data.n_y);
		for (int k = 0; k < data.n_y; ++k) {
		    temp[k] = x[(i + 1) * data.n_y + k];
		}
		valarray<double> y = data.matrix_vector_mult(data.B_w_rows, data.B_w_cols, data.B_w_vals, temp);


		for (int k = 0; k < data.n_y; ++k) {
		    while (countB_y != size_B_y && data.B_y_rows[countB_y] == k) {
			values[count] = -1 * data.B_y_vals[countB_y];
			++count;
			++countB_y;
		    }
		    //change this part
		    while (countA != size_A && data.A_rows[countA] == k) {
			values[count] = data.A_vals[countA] + w * data.B_w_vals[countA];
			++count;
			++countA;
		    }

		    values[count] = -data.b_u[k];
		    ++count;

		    values[count] = y[k];
		    ++count;
		}
	    }
	}
    }


	//no convection
    else {
	//copy A_eq
	long size_A_eq = (data.B_y_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;

	assert(nele_jac == size_A_eq);

	if (values == NULL) {
	    for (int i = 0; i < size_A_eq; ++i) {
		iRow[i] = data.A_eq_rows[i];
		jCol[i] = data.A_eq_cols[i];
	    }
	}
	else {
	    for (int i = 0; i < size_A_eq; ++i) {
		values[i] = data.A_eq_vals[i];
	    }
	}
    }
    return true;
}

//approximate hessian of the lagrange function by hessian of the objective function

bool HEAT_NLP::eval_h(Index n, const Number* x, bool new_x,
	Number obj_factor, Index m, const Number* lambda,
	bool new_lambda, Index nele_hess,
	Index* iRow, Index* jCol, Number * values) {
    assert(data.n_z == n);


    //test if approximation of hessian using only the hessian of the objective function is good enough
    /*
    if (data.convection) {
	return false;
    }
     */
    if (values == NULL) {
	for (int i = 0; i < data.n_z; ++i) {
	    iRow[i] = i;
	    jCol[i] = i;
	}
    }
    else {
	for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
	    values[i] = obj_factor * data.eps;
	}
	//true for convection aswell, as long as W = I
	for (int i = (data.N + 1) * data.n_y; i < data.n_z; ++i) {
	    values[i] = obj_factor * 1;
	}
    }

    return true;
}

void HEAT_NLP::finalize_solution(SolverReturn status,
	Index n, const Number* x, const Number* z_L, const Number* z_U,
	Index m, const Number* g, const Number* lambda,
	Number obj_value,
	const IpoptData* ip_data,
	IpoptCalculatedQuantities * ip_cq) {



    //save result for next optimization step
    for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
	data.y_old[i] = x[i];
    }
    for (int i = 0; i < data.N * data.n_u; ++i) {
	data.u_old[i] = x[(data.N + 1) * data.n_y + i];
    }

    if (data.convection) {
	for (int i = 0; i < data.N * data.n_w; ++i) {
	    data.w_old[i] = x[(data.N + 1) * data.n_y + data.N * data.n_u + i];
	}
    }

    //closed loop cost
    valarray<double> x_new(data.n_y);
    for (int i = 0; i < data.n_y; ++i) {
	x_new[i] = x[data.n_y + i];
    }
    valarray<double> u_new(data.n_u);
    for (int i = 0; i < data.n_u; ++i) {
	u_new[i] = x[(data.N + 1) * data.n_y + i];
    }

    data.closed_loop_cost += data.eps * 0.5 * data.vec_Q_vec(x_new, data.y_ref) + 0.5 * data.vec_R_vec(u_new, data.u_ref);

    if (data.convection) {
	valarray<double> w_new(data.n_w);
	for (int i = 0; i < data.n_w; ++i) {
	    w_new[i] = x[(data.N + 1) * data.n_y + data.N * data.n_u + i];
	}

	data.closed_loop_cost += 0.5 * data.vec_W_vec(w_new);
    }



    //closed loop values
    if (data.closed_values) {
	ofstream ofs_y(data.foldername + "closedloop_y.txt", ofstream::app);
	ofstream ofs_u(data.foldername + "closedloop_u.txt", ofstream::app);
	ofstream ofs_cost(data.foldername + "closedloop_cost.txt", ofstream::app);

	//write initial values
	if (data.iter == 0) {
	    if (data.dim2) {
		valarray<double> temp(data.n_y);

		for (int i = 0; i < data.n_y; ++i) {
		    temp[data.order[i]] = x[i];
		}
		for (int i = 0; i < data.n_y; ++i) {
		    ofs_y << temp[i] << " ";
		}
		ofs_y << endl;
	    }

	    else {
		for (int i = 0; i < data.n_y; ++i) {
		    ofs_y << x[i] << " ";
		}
		ofs_y << endl;
	    }
	}

	//order data, then write
	if (data.dim2) {
	    valarray<double> temp(data.n_y);

	    for (int i = 0; i < data.n_y; ++i) {
		temp[data.order[i]] = x[data.n_y + i];
	    }
	    for (int i = 0; i < data.n_y; ++i) {
		ofs_y << temp[i] << " ";
	    }
	    ofs_y << endl;
	}

	else {
	    for (int i = 0; i < data.n_y; ++i) {
		ofs_y << x[data.n_y + i] << " ";
	    }
	    ofs_y << endl;
	}


	for (int i = 0; i < data.n_u; ++i) {
	    ofs_u << x[(data.N + 1) * data.n_y + i] << " ";
	}
	ofs_u << endl;

	ofs_cost << data.closed_loop_cost << endl;

	ofs_y.close();
	ofs_u.close();
	ofs_cost.close();



	if (data.convection) {
	    ofstream ofs_w(data.foldername + "closedloop_w.txt", ofstream::app);
	    for (int i = 0; i < data.n_w; ++i) {
		ofs_w << x[(data.N + 1) * data.n_y + data.N * data.n_u + i] << " ";
	    }
	    ofs_w << endl;
	    ofs_w.close();
	}
    }


    //open loop values
    if (data.open_values) {
	string filename_y = data.foldername + "openloop_y" + std::to_string(data.iter) + ".txt";
	string filename_u = data.foldername + "openloop_u" + std::to_string(data.iter) + ".txt";
	ofstream ofs_y_open(filename_y, ofstream::trunc);
	ofstream ofs_u_open(filename_u, ofstream::trunc);


	if (data.dim2) {
	    valarray<double> temp((data.N + 1) * data.n_y);

	    for (int k = 0; k < data.N + 1; ++k) {
		for (int i = 0; i < data.n_y; ++i) {
		    temp[data.order[i] + k * data.n_y] = x[k * data.n_y + i];
		}
	    }

	    for (int k = 0; k < data.N + 1; ++k) {
		for (int i = 0; i < data.n_y; ++i) {
		    ofs_y_open << temp[k * data.n_y + i] << " ";
		}
		ofs_y_open << endl;
	    }
	}

	else {
	    for (int k = 0; k < data.N + 1; ++k) {
		for (int i = 0; i < data.n_y; ++i) {
		    ofs_y_open << x[k * data.n_y + i] << " ";
		}
		ofs_y_open << endl;
	    }
	}




	for (int k = 0; k < data.N; ++k) {
	    for (int i = 0; i < data.n_u; ++i) {
		ofs_u_open << x[(data.N + 1) * data.n_y + k * data.n_u + i] << endl;
	    }
	}

	ofs_y_open.close();
	ofs_u_open.close();

	if (data.convection) {
	    string filename_w = data.foldername + "openloop_w" + std::to_string(data.iter) + ".txt";
	    ofstream ofs_w_open(filename_w, ofstream::trunc);
	    for (int k = 0; k < data.N; ++k) {
		for (int i = 0; i < data.n_w; ++i) {
		    ofs_w_open << x[(data.N + 1) * data.n_y + data.N * data.n_u + k * data.n_u + i] << endl;
		}
	    }
	    ofs_w_open.close();
	}
    }




    //count for closed loop
    ++data.iter;
}
