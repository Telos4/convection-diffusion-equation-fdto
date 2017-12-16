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

HEAT_NLP::HEAT_NLP(int N_, string file_A, string file_B, string file_b_u, string file_b_y)
: data(N_, file_A, file_B, file_b_u, file_b_y) {
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
    double eps = 0.001;
    obj_value = data.eval_f(x, n);
    return true;
}

bool HEAT_NLP::eval_grad_f(Index n, const Number* x, bool new_x,
        Number * grad_f) {
    assert(n == 4);

    grad_f[0] = 2 * x[0] * x[3] + x[3]* (x[1] + x[2]);
    grad_f[1] = x[0] * x[3];
    grad_f[2] = x[0] * x[3] + 1;
    grad_f[3] = x[0] * (x[0] + x[1] + x[2]);

    return true;
}

bool HEAT_NLP::eval_g(Index n, const Number* x, bool new_x,
        Index m, Number * g) {
    g[0] = x[0] * x[1] * x[2] * x[3];
    g[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

    return true;
}

bool HEAT_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
        Index m, Index nele_jac, Index* iRow,
        Index *jCol, Number * values) {
    if (values == NULL) {
        iRow[0] = 0;
        jCol[0] = 0;
        iRow[1] = 0;
        jCol[1] = 1;
        iRow[2] = 0;
        jCol[2] = 2;
        iRow[3] = 0;
        jCol[3] = 3;
        iRow[4] = 1;
        jCol[4] = 0;
        iRow[5] = 1;
        jCol[5] = 1;
        iRow[6] = 1;
        jCol[6] = 2;
        iRow[7] = 1;
        jCol[7] = 3;
    }
    else {
        values[0] = x[1] * x[2] * x[3]; // 0,0
        values[1] = x[0] * x[2] * x[3]; // 0,1
        values[2] = x[0] * x[1] * x[3]; // 0,2
        values[3] = x[0] * x[1] * x[2]; // 0,3

        values[4] = 2 * x[0]; // 1,0
        values[5] = 2 * x[1]; // 1,1
        values[6] = 2 * x[2]; // 1,2
        values[7] = 2 * x[3]; // 1,3
    }
    return true;
}

bool HEAT_NLP::eval_h(Index n, const Number* x, bool new_x,
        Number obj_factor, Index m, const Number* lambda,
        bool new_lambda, Index nele_hess,
        Index* iRow, Index* jCol, Number * values) {
    if (values == NULL) {
        Index iter = 0;
        for (Index iInd = 0; iInd < n; ++iInd) {
            for (Index jInd = 0; jInd <= iInd; ++jInd) {
                iRow[iter] = iInd;
                jCol[iter] = jInd;
                ++iter;
            }
        }

        assert(iter == nele_hess);
    }
    else {
        // return the values. This is a symmetric matrix, fill the lower left
        // triangle only

        // fill the objective portion
        values[0] = obj_factor * (2 * x[3]); // 0,0

        values[1] = obj_factor * (x[3]); // 1,0
        values[2] = 0; // 1,1

        values[3] = obj_factor * (x[3]); // 2,0
        values[4] = 0; // 2,1
        values[5] = 0; // 2,2

        values[6] = obj_factor * (2 * x[0] + x[1] + x[2]); // 3,0
        values[7] = obj_factor * (x[0]); // 3,1
        values[8] = obj_factor * (x[0]); // 3,2
        values[9] = 0; // 3,3


        // add the portion for the first constraint
        values[1] += lambda[0] * (x[2] * x[3]); // 1,0

        values[3] += lambda[0] * (x[1] * x[3]); // 2,0
        values[4] += lambda[0] * (x[0] * x[3]); // 2,1

        values[6] += lambda[0] * (x[1] * x[2]); // 3,0
        values[7] += lambda[0] * (x[0] * x[2]); // 3,1
        values[8] += lambda[0] * (x[0] * x[1]); // 3,2

        // add the portion for the second constraint
        values[0] += lambda[1] * 2; // 0,0

        values[2] += lambda[1] * 2; // 1,1

        values[5] += lambda[1] * 2; // 2,2

        values[9] += lambda[1] * 2; // 3,3
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

    // For this example, we write the solution to the console
    printf("\n\nSolution of the primal variables, x\n");
    for (Index i = 0; i < n; i++) {
        printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (Index i = 0; i < n; i++) {
        printf("z_L[%d] = %e\n", i, z_L[i]);
    }
    for (Index i = 0; i < n; i++) {
        printf("z_U[%d] = %e\n", i, z_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj_value);
}

