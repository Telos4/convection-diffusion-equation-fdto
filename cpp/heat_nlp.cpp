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

/*
 * n: number of variables
 * m: number of constraints 
 */
bool HEAT_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	Index& nnz_h_lag, IndexStyleEnum& index_style) {

    n = data.n_z;
    m = data.N * data.n_y;

    if (data.convection) {
	//we assume that A and B_w have the same size and pattern
	nnz_jac_g = (data.B_y_rows.size() + data.A_rows.size() + data.b_u.size() + data.n_y) * data.N;
    }
    else {
	//jac(g) = A_eq, compare docu_older.pdf
	nnz_jac_g = (data.B_y_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N;
    }

    //for the non zero entries in the hessian of the lagrange fucntion we assume that Q, R, W are unit matrices
    //don't need this when using approximation of hessian of the lagrangian with convection on
    nnz_h_lag = data.n_z;


    index_style = TNLP::C_STYLE;

    return true;
}

bool HEAT_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
	Index m, Number* g_l, Number* g_u) {


    double inf = 1e19;

    //The initial state is calculated in the previous step.
    //free_init_value allows the first step of the optimization to start with a 
    //free initial value and not the reference solution (see the constructor in matrixopp.cpp.

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
	//all between 0 and 1, assuming an [0,1]x[0,1] grid
	double left = data.boundary_left;
	double right = data.boundary_right;
	double top = data.boundary_top;
	double bot = data.boundary_bot;

	for (int k = 1; k < data.N + 1; ++k) {
	    for (int i = 0; i < data.n_y; ++i) {
		//dof_x and dof_y give the indices of the state variables.
		//From these we need to reconstruct the actual (x,y) position in our grid.
		//We achieve this by dividing by our discretization parameter.
		double x = data.dof_x[i] / data.discretization_n;
		double y = data.dof_y[i] / data.discretization_n;

		//test if the point (x,y) is in our constrain area
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
	//there are no state boundaries in the left and right quarter of the [0,1] domain
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

	    //free_init_value means we want to calculate the reference solution, where we know the exact values for y_out
	    if (data.free_init_value && data.inexact_data) {
		g_u[data.n_y * k + i] = g_l[data.n_y * k + i] = data.b_y[i] * (eval_y_out(data.iter + k) + data.error[data.iter + k]);
	    }
	    else {
		g_u[data.n_y * k + i] = g_l[data.n_y * k + i] = data.b_y[i] * eval_y_out(data.iter + k);
	    }
	}
    }

    return true;
}

//where to put this function?

double eval_y_out(int i) {
    return 0.3 * sin(0.1 * i);
}

/*
 * The point we want to start the optimization at is the (one timestep) shifted solution from the previous step.
 * For the uninitialized values at the end we use the last values of the previous step for a second time.
 */
bool HEAT_NLP::get_starting_point(Index n, bool init_x, Number* x,
	bool init_z, Number* z_L, Number* z_U,
	Index m, bool init_lambda,
	Number * lambda) {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);



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

/*
 * see matrixop.cpp
 * need to set obj_value
 * Evaluate the cost functional 
 * J(y,u,w) = 0.5 * sum(k=0..N)[eps*(y_k-y_ref).Q.(y_k-y_ref)] + 0.5 * sum(k=0..N-1)[(u_k-u_ref.R.(u_k-u_ref) + (w_k.W.w_k)]
 */
bool HEAT_NLP::eval_f(Index n, const Number* x, bool new_x, Number & obj_value) {

    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }

    obj_value = data.eval_f(z);
    return true;
}

/*
 * see matrixop.cpp
 * need to set grad_f
 * Evaluate the gradient of the cost functional
 * gradJ(y,u,w) =   [eps*Q.(y_0-y_ref)	]
 *		    [	    ...		]
 *		    [eps*Q.(y_N-y_ref	]
 *		    [R.(u_0-u_ref)	]
 *		    [	    ...		]
 *		    [R.(u_(N-1)-u_ref	]
 *		    [W.w_0		]
 *		    [	    ...		]
 *		    [W.w_(N-1)		]
 */
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

/*
 * evaluate the constraints, equation (18) in docu_old.pdf
 * The constant terms b_(y,out)*y_out are implemented as equality constraints for g in get_bounds_info().
 * Split g(z) = A_eq*(y,u) + h(y,w) in the linear part A_eq*(y,u) and in the case of convection 
 * the nonlinear part h(y,w) =	[w_0 B_w*y_1]
 *				[     .	    ]
 *				[w_N B_w*y_N]
 * where y_i is the i-th state of the open loop.
 * 
 */
bool HEAT_NLP::eval_g(Index n, const Number* x, bool new_x,
	Index m, Number * g) {


    valarray<double> z(n);
    for (int i = 0; i < n; ++i) {
	z[i] = x[i];
    }

    //The dimension do not fit in case of convection, but the way the matrix vector multiplication is implemented allows us to ignore this.
    //The convection will therefore have no effect here.
    valarray<double> eval1 = data.matrix_vector_mult(data.A_eq_rows, data.A_eq_cols,
	    data.A_eq_vals, z, data.N * data.n_y);

    //nonlinear part
    valarray<double> eval2(0.0, m);
    if (data.convection) {
	for (int i = 0; i < data.N; ++i) {
	    valarray<double> temp = z[slice((i + 1) * data.n_y, data.n_y, 1)];
	    double w_i = z[(data.N + 1) * data.n_y + data.N * data.n_u + i];
	    eval2[slice(i * data.n_y, data.n_y, 1)] = w_i * data.matrix_vector_mult(data.B_w_rows, data.B_w_cols, data.B_w_vals, temp, data.n_y);
	}
    }

    //eval2 = 0 if convection is off
    for (int i = 0; i < data.N * data.n_y; ++i) {
	g[i] = eval1[i] + eval2[i];
    }
    return true;
}

/*
 * jac_g = 
 * [-B_y    A+w_0*B_w					    -b_u			B_w*y_1				    ]
 * [	    -B_y	A+w_1*B_w				    -b_u			    B_w*y_2		    ]
 * [			       ..					..				    ..		    ]
 * [				    ..					    ..					..	    ]
 * [				    -B_y    A+w_(N-1)*B_w			-b_u				    B_w*y_N ]
 * 
 * 
 * Assumption: A and B_w have the same pattern
 * 
 * The idea is the same as the creation of A_eq in matrixop.cpp.
 * The only new thing here are the w_i*B_w and B_w*y_i terms
 */

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
		valarray<double> y = data.matrix_vector_mult(data.B_w_rows, data.B_w_cols, data.B_w_vals, temp, data.n_y);


		for (int k = 0; k < data.n_y; ++k) {
		    while (countB_y != size_B_y && data.B_y_rows[countB_y] == k) {
			values[count] = -1 * data.B_y_vals[countB_y];
			++count;
			++countB_y;
		    }

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

    /*
     * The control and convection term for the next time step are calculated with exact boundary data. 
     * In the case of inexcact boundary data we need to recalculate the resulting state.
     * The formula can be found in docu_old.pdf (10). 
     * A*y_(k+1) + w_k*B_w.y_(k+1) = B_y.y_k + b_u*u_k + b_y*y_(out,k) , 
     * with (.) as matrix vector multiplication and (*) as scalar multiplication
     * 
     * save this new state in y_new. If inexact_data is off, y_new contains the new state from the open loop 
     * 
     * The linear system is solved with Eigen: https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
     * Currently a sparse LU decompositoin is used, because we have a one dimensional small problem
     */
    valarray<double> y_new(data.n_y);
    if (data.inexact_data && !data.free_init_value) {
	double u = x[(data.N + 1) * data.n_y];

	//I don't know why I need this minus, but without it the convection works in the wrong direction.
	//This is especially strange because this minus is not needed in eval_g() and eval_jac_g
	//And it shouldn't be needed because the direction of the convection is set when creating B_w in heat.py
	double w = -x[(data.N + 1) * data.n_y + data.N * data.n_u];
	double y_out = eval_y_out(data.iter) + data.error[data.iter];

	typedef Eigen::Triplet<double> T;

	std::vector<T> tripletList;
	tripletList.reserve(2 * data.A_rows.size());

	for (int i = 0; i < data.A_rows.size(); ++i) {
	    tripletList.push_back(T(data.A_cols[i], data.A_rows[i], data.A_vals[i]));
	}
	for (int i = 0; i < data.B_w_rows.size(); ++i) {
	    tripletList.push_back(T(data.B_w_cols[i], data.B_w_rows[i], w * data.B_w_vals[i]));
	}

	Eigen::SparseMatrix<double> mat(data.n_y, data.n_y);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());

	Eigen::VectorXd b(data.n_y), y;

	valarray<double> y_k(data.n_y);

	for (int i = 0; i < data.n_y; ++i) {
	    y_k[i] = x[i];
	}

	valarray<double> rhs = data.matrix_vector_mult(data.B_y_rows, data.B_y_cols, data.B_y_vals, y_k, data.n_y);

	for (int i = 0; i < data.n_y; ++i) {
	    rhs[i] += data.b_u[i] * u + data.b_y[i] * y_out;
	    b[i] = rhs[i];
	}

	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	//compute LU decomposition
	solver.compute(mat);
	if (solver.info() != Eigen::Success) {
	    cout << "decomposition failed, in finalize_solution -> Eigen" << endl;
	    exit(1);
	}
	//solve system
	y = solver.solve(b);
	if (solver.info() != Eigen::Success) {
	    cout << "solving failed, in finalize_solution -> Eigen" << endl;
	    exit(1);
	}

	for (int i = 0; i < data.n_y; ++i) {
	    y_new[i] = y[i];
	}
    }

    else {
	for (int i = 0; i < data.n_y; ++i) {
	    y_new[i] = x[data.n_y + i];
	}
    }

    //save results for next optimization step
    for (int i = 0; i < (data.N + 1) * data.n_y; ++i) {
	data.y_old[i] = x[i];
    }
    //set first state in the next step
    if (data.inexact_data) {
	for (int i = 0; i < data.n_y; ++i) {
	    data.y_old[data.n_y + i] = y_new[i];
	}
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
	x_new[i] = y_new[i];
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

    //hypothetical closed loop cost for the reference solution
    if (data.free_init_value && data.closed_values) {
	
	ofstream ofs_cost(data.foldername + "closedloop_cost.txt", ofstream::app);
	
	double closedloop_reference = 0;
	
	for (int k = 0; k < data.N; ++k) {
	    valarray<double> y(data.n_y);
	    for (int i = 0; i < data.n_y; ++i) {
		y[i] = x[k * data.n_y + i];
	    }
	    valarray<double> u(data.n_u);
	    for (int i = 0; i < data.n_u; ++i) {
		u[i] = x[(data.N + 1) * data.n_y + k * data.n_u + i];
	    }

	    closedloop_reference += data.eps * 0.5 * data.vec_Q_vec(y, data.y_ref) + 0.5 * data.vec_R_vec(u, data.u_ref);

	    if (data.convection) {
		valarray<double> w(data.n_w);
		for (int i = 0; i < data.n_w; ++i) {
		    w[i] = x[(data.N + 1) * data.n_y + data.N * data.n_u + k * data.n_u + i];
		}

		closedloop_reference += 0.5 * data.vec_W_vec(w);
	    }
	    ofs_cost << closedloop_reference << endl;

	}

	ofs_cost.close();
    }

    //write closed loop values
    if (data.closed_values) {
	ofstream ofs_y(data.foldername + "closedloop_y.txt", ofstream::app);
	ofstream ofs_u(data.foldername + "closedloop_u.txt", ofstream::app);

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
		temp[data.order[i]] = y_new[i];
	    }
	    for (int i = 0; i < data.n_y; ++i) {
		ofs_y << temp[i] << " ";
	    }
	    ofs_y << endl;
	}

	//1d
	else {
	    for (int i = 0; i < data.n_y; ++i) {
		ofs_y << y_new[i] << " ";
	    }
	    ofs_y << endl;
	}


	for (int i = 0; i < data.n_u; ++i) {
	    ofs_u << x[(data.N + 1) * data.n_y + i] << " ";
	}
	ofs_u << endl;

	//the hypothetical closed loop cost for the reference solution is already written above
	if (!data.free_init_value) {
	    ofstream ofs_cost(data.foldername + "closedloop_cost.txt", ofstream::app);
	    ofs_cost << data.closed_loop_cost << endl;
	    ofs_cost.close();
	}

	ofs_y.close();
	ofs_u.close();



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

	//write state in 2d
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

	    //write state in 1d
	else {
	    for (int k = 0; k < data.N + 1; ++k) {
		for (int i = 0; i < data.n_y; ++i) {
		    ofs_y_open << x[k * data.n_y + i] << " ";
		}
		ofs_y_open << endl;
	    }
	}

	//write u
	for (int k = 0; k < data.N; ++k) {
	    for (int i = 0; i < data.n_u; ++i) {
		ofs_u_open << x[(data.N + 1) * data.n_y + k * data.n_u + i] << endl;
	    }
	}

	ofs_y_open.close();
	ofs_u_open.close();

	//write w
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
