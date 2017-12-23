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

HEAT_NLP::HEAT_NLP() {
}

//destructor

HEAT_NLP::~HEAT_NLP() {
}
double eps = 0.001;
double y_ref = 0.5;
double u_ref = 0.5;

HEAT_NLP::HEAT_NLP(int N_, string file_A, string file_B, string file_b_u, string file_b_y,
        double eps, double y_ref, double u_ref)
: data(N_, file_A, file_B, file_b_u, file_b_y, eps, y_ref, u_ref) {
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
    double y_upper = 1;
    double y_lower = -1;
    double inf = 1e19;
    Index left = (int) (data.n_y / 4);
    Index right = (int) (data.n_y * 3 / 4);

    //bound for y left
    for (Index i = 0; i < left; ++i) {
        x_u[i] = inf;
        x_l[i] = -inf;
    }

    //bound for y center
    for (Index i = left; i < right; ++i) {
        x_u[i] = y_upper;
        x_l[i] = y_lower;
    }

    //bound for y right
    for (Index i = right; i < (data.N + 1) * data.n_y; ++i) {
        x_u[i] = inf;
        x_l[i] = -inf;
    }

    //bound for u
    for (Index i = (data.N + 1) * data.n_y; i < data.n_z; ++i) {
        x_u[i] = u_upper;
        x_l[i] = u_lower;
    }

    //equality constraints for g
    for (Index i = 0; i < data.N; ++i) {
        for (Index k = 0; i < data.n_y; ++i) {
            g_u[data.N * i + k] = g_l[data.N * i + k] = data.b_y[k] * eval_y_out(i);
        }
    }
    return true;
}

//where to put this function?

double eval_y_out(int i) {
    return 0.5 + 0.3 * sin(0.1 * i);
}

bool HEAT_NLP::get_starting_point(Index n, bool init_x, Number* x,
        bool init_z, Number* z_L, Number* z_U,
        Index m, bool init_lambda,
        Number * lambda) {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    for (Index i = 0; i < data.n_y; ++i) {
        x[i] = 0;
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
    valarray<double> eval = data.matrix_vektor_mult(data.A_eq_rows, data.A_eq_cols,
            data.A_eq_vals, z);

    for (int i = 0; i < n; ++i) {
        g[i] = eval[i];
    }
    return true;
}

bool HEAT_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number * values) {

    int size_A_eq = (data.B_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;
    if (values == NULL) {
        for (int i = 0; i < size_A_eq; ++i) {
            iRow[i] = data.A_eq_rows[i];
            jCol[i] = data.A_eq_cols[i];
        }
    }
    else {
        for (int i = 0; i < size_A_eq; ++i) {
            values[i] = data.A_eq_rows[i];
        }
    }
    return true;
}

bool HEAT_NLP::eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number * values) {
    if (values == NULL) {
        for (int i = 0; i < data.n_z; ++i) {
            iRow[i] = i;
            jCol[i] = i;
        }
    }
    else {
        for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
            values[i] = data.eps;
        }
        for (int i = 0; i < data.N * data.n_u; ++i) {
            values[i] = 1;
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
    // here is where we would store the solution to variables, or write to a file, etc
    // so we could use the solution.
    ofstream ofs_y("solution_y.txt");
    ofstream ofs_u("solution_u.txt");

    for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
        ofs_y << x[i] << endl;
    }
    for (int i = 0; i < data.N * data.n_u; ++i) {
        ofs_y << x[(data.N + 1) * data.n_y + i] << endl;
    }
    ofs_y.close();
    ofs_u.close();
}

