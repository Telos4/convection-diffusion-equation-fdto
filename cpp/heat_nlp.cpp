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

/*
HEAT_NLP::HEAT_NLP(int N_, string file_A, string file_B, string file_b_u, string file_b_y,
	double eps, double y_ref, double u_ref)
: data(N_, file_A, file_B, file_b_u, file_b_y, eps, y_ref, u_ref) {
}
 */

HEAT_NLP::HEAT_NLP(MATRIXOP &data_)
: data(data_) {

}

bool HEAT_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	Index& nnz_h_lag, IndexStyleEnum& index_style) {
    n = data.n_z;
    m = data.N * data.n_y;

    //jac(g) = A_eq
    nnz_jac_g = (data.B_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;

    //mit Q = I, R = I
    nnz_h_lag = data.n_z;
    index_style = TNLP::C_STYLE;

    return true;
}

bool HEAT_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
	Index m, Number* g_l, Number* g_u) {
    double u_upper = 0.75;
    double u_lower = 0.25;
    double y_upper = 0.65;
    double y_lower = 0.35;
    double inf = 1e19;
    Index left = floor(data.n_y / 4);
    //Index right = (int) (data.n_y * 3 / 4);


    //initial value, testing
    for (Index i = 0; i < data.n_y; ++i) {
	x_u[i] = data.y_old[data.n_y + i];
	x_l[i] = data.y_old[data.n_y + i];
    }


    for (Index k = 1; k < data.N + 1; ++k) {
	//bound for y left, right
	for (Index i = 0; i < left; ++i) {
	    x_u[k * data.n_y + i] = x_u[(k + 1) * data.n_y - i - 1] = inf;
	    x_l[k * data.n_y + i] = x_l[(k + 1) * data.n_y - i - 1] = -inf;
	}

	//bound for y center
	for (Index i = left; i < data.n_y - left + 1; ++i) {
	    x_u[k * data.n_y + i] = y_upper;
	    x_l[k * data.n_y + i] = y_lower;
	}
	/*
		//bound for y right
		for (Index i = right; i < data.n_y; ++i) {
		    x_u[k * data.n_y + i] = inf;
		    x_l[k * data.n_y + i] = -inf;
		}
	 */
    }
    //bound for u
    for (Index i = (data.N + 1) * data.n_y; i < data.n_z; ++i) {
	x_u[i] = u_upper;
	x_l[i] = u_lower;
    }

    //equality constraints for g
    for (Index i = 0; i < data.N; ++i) {
	for (Index k = 0; k < data.n_y; ++k) {
	    g_u[data.n_y * i + k] = g_l[data.n_y * i + k] = data.b_y[k] * eval_y_out(data.iter + i);
	}
    }
    return true;
}

//where to put this function?

double eval_y_out(int i) {
    return 0.5 + 0.292 * sin(0.1 * i);
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
    //for the uninitialized values use the last known value y_old[data.N*data.n_y - 1]
    for (int i = 0; i < data.N * data.n_y; ++i) {
	x[i] = data.y_old[data.n_y + i];
    }
    for (int i = data.N * data.n_y; i < (data.N + 1) * data.n_y; ++i) {
	x[i] = data.y_old[data.N * data.n_y - 1];
    }

    //wrong for 2 dimensions
    for (int i = 0; i < data.N * data.n_u - 1; ++i) {
	x[(data.N + 1) * data.n_y + i] = data.u_old[i + 1];
    }
    for (int i = data.N * data.n_u - 1; i < data.N * data.n_u; ++i) {
	x[(data.N + 1) * data.n_y + i] = data.u_old[data.N * data.n_u - 1];
    }

    return true;
}

//hier oder in data?

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
    //data.print_matrix(data.A_eq_rows, data.A_eq_cols, data.A_eq_vals);
    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }
    valarray<double> eval = data.matrix_vektor_mult(data.A_eq_rows, data.A_eq_cols,
	    data.A_eq_vals, z);
    //data.print_vector(z);
    //cout << endl;
    //data.print_vector(eval);
    assert(data.N * data.n_y == m);
    for (int i = 0; i < data.N * data.n_y; ++i) {
	g[i] = eval[i];
    }
    return true;
}

bool HEAT_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
	Index m, Index nele_jac, Index* iRow,
	Index *jCol, Number * values) {

    int size_A_eq = (data.B_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;

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
    return true;
}

bool HEAT_NLP::eval_h(Index n, const Number* x, bool new_x,
	Number obj_factor, Index m, const Number* lambda,
	bool new_lambda, Index nele_hess,
	Index* iRow, Index* jCol, Number * values) {
    assert(data.n_z == n);
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

    //count for outer loop
    ++data.iter;

    for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
	data.y_old[i] = x[i];
    }
    for (int i = 0; i < data.N * data.n_u; ++i) {
	data.u_old[i] = x[(data.N + 1) * data.n_y + i];
    }



    //open loop values
    
    string filename_y = "results50/openloop_y" + std::to_string(data.iter) + ".txt";
    //string filename_u = "results/openloop_u" + std::to_string(data.iter) + ".txt";
    ofstream ofs_y_open(filename_y);
    //ofstream ofs_u_open(filename_u);
    
    for (int k = 0; k < data.N + 1; ++k) {
	for (int i = 0; i < data.n_y; ++i) {
	    ofs_y_open << x[k * data.n_y + i] << " ";
	}
	ofs_y_open << endl;
    }
    
    /*
    for (int k = 0; k < data.N; ++k) {
	for (int i = 0; i < data.n_u; ++i) {
	    ofs_u_open << x[(data.N + 1) * data.n_y + k * data.n_u + i] << " ";
	}
	ofs_u_open << endl;
    }

    ofs_u_open.close();
    */
    ofs_y_open.close();
    
    
    /*   
    ofstream ofs_y("solution_y.txt", ofstream::app);
    ofstream ofs_u("solution_u.txt", ofstream::app);
    //closed loop values

    for (int i = 0; i < data.n_y; ++i) {
	ofs_y << x[data.n_y + i] << " ";
    }
    ofs_y << endl;

    for (int i = 0; i < data.n_u; ++i) {
	ofs_u << x[(data.N + 1) * data.n_y + i] << " ";
    }
    ofs_u << endl;

    ofs_y.close();
    ofs_u.close();
     */

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
}

