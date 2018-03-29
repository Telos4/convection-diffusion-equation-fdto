/*
 * matrixop.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include "matrixop.hpp"

using namespace std;
// check matrixop.hpp for an explanation of the variables
// constructor

MATRIXOP::MATRIXOP() {
}

MATRIXOP::MATRIXOP(int N_, int steps_, string file_A, string file_B_y, string file_B_w, string file_b_u, string file_b_y,
	string file_dof_x, string file_dof_y, double eps_, double y_ref_, double u_ref_, double u_upper_, double u_lower_, double y_upper_, double y_lower_, double w_upper_, double w_lower_,
	double boundary_right_, double boundary_left_, double boundary_top_, double boundary_bot_, double std_dev_,
	bool dim2_, bool convection_, bool closed_values_, bool open_values_, bool free_init_value_, bool inexact_data_, string result_folder_, string result_folder_prefix_) {

    dim2 = dim2_;
    convection = convection_;
    inexact_data = inexact_data_;

    if (N_ == 0) {
	cout << "N = 0 does not make sense, exiting" << endl;
	exit(1);
    }

    N = N_;
    steps = steps_;

    read_matrix(file_A, A_rows, A_cols, A_vals);
    read_matrix(file_B_y, B_y_rows, B_y_cols, B_y_vals);
    read_vector(file_b_u, b_u);
    read_vector(file_b_y, b_y);

    if (dim2) {
	read_vector(file_dof_x, dof_x);
	read_vector(file_dof_y, dof_y);
    }
    if (convection) {
	read_matrix(file_B_w, B_w_rows, B_w_cols, B_w_vals);
    }


    n_y = b_u.size();
    n_u = 1;

    if (convection) {
	n_w = 1;
	w_upper = w_upper_;
	w_lower = w_lower_;
    }
    else {
	n_w = 0;
    }

    n_z = (N + 1) * n_y + N * n_u + N * n_w;

    eps = eps_;
    y_ref = y_ref_;
    u_ref = u_ref_;
    u_upper = u_upper_;
    u_lower = u_lower_;
    y_upper = y_upper_;
    y_lower = y_lower_;
    boundary_right = boundary_right_;
    boundary_left = boundary_left_;
    boundary_top = boundary_top_;
    boundary_bot = boundary_bot_;
    std_dev = std_dev_;

    iter = 0;
    closed_loop_cost = 0;

    //assuming an n x n grid in heat2d.py, not used in 1d case 
    discretization_n = (int) sqrt(n_y);

    closed_values = closed_values_;
    open_values = open_values_;

    free_init_value = free_init_value_;

    //initialize A_eq
    A_eq();

    if (dim2) {
	initialize_order();
    }


    y_old.resize((N + 1) * n_y);
    u_old.resize(N * n_u);

    //start the optimization at the reference solution
    for (int i = 0; i < (N + 1) * n_y; ++i) {
	y_old[i] = y_ref;
    }
    for (int i = 0; i < N * n_u; ++i) {
	u_old[i] = u_ref;
    }
    if (convection) {
	w_old.resize(N * n_w);
	for (int i = 0; i < N * n_w; ++i) {
	    w_old[i] = 0.0;
	}
    }

    //fill the error vector
    if (inexact_data) {
	error.resize(N + steps);
	
	//seed 2018
	std::default_random_engine generator(2018);
	std::normal_distribution<double> distribution(0.0, std_dev);
	
	for (int i = 0; i < N + steps; ++i) {
	    error[i] = distribution(generator);
	}
    }

    //create folder in which the result will be stored
    if (closed_values || open_values) {
	result_folder = result_folder_;
	create_folder(result_folder_prefix_);
    }
}

//destructor

MATRIXOP::~MATRIXOP() {
}

/*
 * Read the matrix saved in 'filename' to the three equally sized valarray rows, cols, vals.
 * The expected format is:
 * 2 lines of comments (info about matrix etc.)
 * x, y, z (x, y are the dimension of the matrix, z is the number of saved entries)
 * actual matrix
 */
void MATRIXOP::read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) {
    ifstream ifs(filename.c_str(), ifstream::in);

    if (ifs.is_open()) {
	//skip comments in .mtx files (first two lines)
	string trash;

	getline(ifs, trash);
	getline(ifs, trash);

	int n, i = 0;
	//third value in file is size
	ifs >> n >> n >> n;
	rows.resize(n);
	cols.resize(n);
	vals.resize(n);

	//then the actual values matrix matrix follows
	while (ifs >> rows[i] >> cols[i] >> vals[i]) {
	    ++i;
	}
	//counting in file starts at 1
	for (int i = 0; i < n; ++i) {
	    rows[i] -= 1;
	    cols[i] -= 1;
	}
    }
    else {
	cout << endl << "could'nt open  " << filename << " ... terminating" << endl;
	exit(1);
    }
}

/*
 * Read the vector saved in 'filename' to the valarray vals.
 * The first entry is expected to be the length of the vector
 */
void MATRIXOP::read_vector(string filename, valarray<double> &vals) {
    //ifstream input(filename);
    ifstream ifs(filename.c_str(), ifstream::in);

    if (ifs.is_open()) {
	int n, i = 0;
	ifs >> n;
	vals.resize(n);
	while (ifs >> vals[i]) {
	    ++i;
	}
    }
    else {
	cout << endl << "could'nt open  " << filename << " ... terminating" << endl;
	exit(1);
    }
}

void MATRIXOP::print_matrix(valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) const {
    for (unsigned i = 0; i < rows.size(); ++i) {
	cout << rows[i] << " " << cols[i] << " " << vals[i] << endl;
    }
}

void MATRIXOP::print_vector(valarray<double> &vals) const {
    for (unsigned i = 0; i < vals.size(); ++i) {
	cout << vals[i] << endl;
    }
}

/*
 * Create the folder where the folders with the results shall be stored in.
 * Then create the folders with the actual results inside this one.
 * result_folder_prefix will be applied to the folder with the results
 */
void MATRIXOP::create_folder(string result_folder_prefix) {
    std::cout << "creating folder!" << std::endl;
    time_t t = time(0); // get time now
    struct tm * now = localtime(& t);


    foldername = result_folder + result_folder_prefix + to_string(now->tm_year + 1900) + "-" + to_string(now->tm_mon + 1) + "-" + to_string(now->tm_mday)
	    + "_" + to_string(now->tm_hour) + "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec) + "/";

    //cout << boost::filesystem::current_path().string() << endl;
    boost::filesystem::path p(foldername);

    //try to create the folder, where the folders with the results shall be stored in, no error if folder exists already
    try {
	boost::filesystem::create_directory(result_folder);
    }
    catch (boost::filesystem::filesystem_error &e) {
	std::cerr << e.what() << '\n';
	exit(1);
    }
    //try to create the folder where the actual results will be stored
    try {
	boost::filesystem::create_directory(p);
    }
    catch (boost::filesystem::filesystem_error &e) {
	std::cerr << e.what() << '\n';
	exit(1);
    }
}

/*
 * Evaluate the cost functional 
 * J(y,u,w) = 0.5 * sum(k=0..N)[eps*(y_k-y_ref).Q.(y_k-y_ref)] + 0.5 * sum(k=0..N-1)[(u_k-u_ref.R.(u_k-u_ref) + (w_k.W.w_k)]
 * The matrices Q, R, W in the cost functional are assumed to be unit matrices.
 * Should this change adapt the functions vec_X_vec which realize the vector-matrix-vector multiplication
 */
double MATRIXOP::eval_f(valarray<double> &x) {
    //Q, R, W unit matrices
    double re = 0;

    //sum over yQy
    for (int k = 0; k < N + 1; ++k) {
	re += eps * 0.5 * vec_Q_vec(x[slice(k*n_y, n_y, 1)], y_ref);
    }
    //sum over uRu
    for (int k = 0; k < N; ++k) {
	re += 0.5 * vec_R_vec(x[slice((N + 1) * n_y + k*n_u, n_u, 1)], u_ref);
    }

    if (convection) {
	//sum over wWw
	for (int k = 0; k < N; ++k) {
	    re += 0.5 * vec_W_vec(x[slice((N + 1) * n_y + N * n_u + k*n_w, n_w, 1)]);
	}
    }
    return re;
}

double MATRIXOP::vec_Q_vec(valarray<double> y, double y_ref) {
    double re = 0;
    for (unsigned i = 0; i < y.size(); ++i) {
	re += (y[i] - y_ref)*(y[i] - y_ref);
    }

    return re;
}

double MATRIXOP::vec_R_vec(valarray<double> y, double y_ref) {
    double re = 0;
    for (unsigned i = 0; i < y.size(); ++i) {
	re += (y[i] - y_ref)*(y[i] - y_ref);
    }

    return re;
}

double MATRIXOP::vec_W_vec(valarray<double> y) {
    double re = 0;
    for (unsigned i = 0; i < y.size(); ++i) {
	re += y[i] * y[i];
    }

    return re;
}

/*
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
 * The matrices Q, R, W in the cost functional are assumed to be unit matrices.
 * Should this change adapt the functions X_vec which realize the matrix-vector multiplication 
 * and in the case of X = Q, R also the difference to the reference solution.
 */
valarray<double> MATRIXOP::eval_grad_f(valarray<double> z) {
    //Q, R, W unit matrices for now
    valarray<double> re(n_z);

    for (int i = 0; i < N + 1; ++i) {
	re[slice(i*n_y, n_y, 1)] = eps * Q_vec(z[slice(i*n_y, n_y, 1)]);
    }

    for (int i = 0; i < N; ++i) {
	re[slice((N + 1) * n_y + i*n_u, n_u, 1)] = R_vec(z[slice((N + 1) * n_y + i*n_u, n_u, 1)]);
    }

    if (convection) {
	for (int i = 0; i < N; ++i) {
	    re[slice((N + 1) * n_y + N * n_u + i * n_w, n_w, 1)] = W_vec(z[slice((N + 1) * n_y + N * n_u + i*n_w, n_w, 1)]);
	}
    }
    return re;
}

valarray<double> MATRIXOP::Q_vec(valarray<double> y) {
    valarray<double> re(y.size());

    for (unsigned i = 0; i < y.size(); ++i) {
	re[i] = y[i] - y_ref;
    }

    return re;
}

valarray<double> MATRIXOP::R_vec(valarray<double> u) {
    valarray<double> re(u.size());

    for (unsigned i = 0; i < u.size(); ++i) {
	re[i] = u[i] - u_ref;
    }
    return re;
}

valarray<double> MATRIXOP::W_vec(valarray<double> w) {
    valarray<double> re(w.size());

    for (unsigned i = 0; i < w.size(); ++i) {
	re[i] = w[i];
    }
    return re;
}


//parameters have size (data.B_rows.size() + data.A_rows.size() + data.b_u.size()) * data.N

/*
 * Create the matrix A_eq found in equation (12) in docu_older.pdf 
 * A_eq =   [ -B_y    A                     -b_u             ]
 *          [      -B_y    A                    -b_u         ]
 *          [                 ...                     ...    ]
 *          [                   -B_y    A                -b_u]
 * 
 * The matrix is created in N blocks [.. , -B_y, A, .. , -b_u, ..] and in each block row by row
 */
void MATRIXOP::A_eq() {
    long size_A_eq = (B_y_rows.size() + A_rows.size() + b_u.size()) * N;
    long size_B = B_y_rows.size();
    long size_A = A_rows.size();
    A_eq_rows.resize(size_A_eq);
    A_eq_cols.resize(size_A_eq);
    A_eq_vals.resize(size_A_eq);

    //whenever we write a value to A_eq we increment count, to know where the next value shall be written to
    int count = 0;
    for (int i = 0; i < N; ++i) {
	//countA and countB are used to keep track of the current position in A and B_y
	int countA = 0;
	int countB = 0;

	//iterate over the n_y lines of one block of A_eq
	for (int k = 0; k < n_y; ++k) {

	    //check if every value of B_y is written already and if we are still in the right row
	    while (countB != size_B && B_y_rows[countB] == k) {
		A_eq_rows[count] = i * n_y + k;
		A_eq_cols[count] = i * n_y + B_y_cols[countB];
		A_eq_vals[count] = -1 * B_y_vals[countB];
		++count;
		++countB;
	    }

	    //check if every value of A is written already and if we are still in the right row
	    while (countA != size_A && A_rows[countA] == k) {
		A_eq_rows[count] = i * n_y + k;
		A_eq_cols[count] = (i + 1) * n_y + A_cols[countA];
		A_eq_vals[count] = A_vals[countA];
		++count;
		++countA;
	    }

	    A_eq_rows[count] = i * n_y + k;
	    A_eq_cols[count] = (N + 1) * n_y + i;
	    A_eq_vals[count] = -b_u[k];
	    ++count;
	}
    }
    cout << endl;
}

void MATRIXOP::initialize_order() {
    //assumption: we have a n x n grid
    //in dof_x and dof_y we have the x and y indices of every degree of freedom (the state variables) of the unit square mesh from the pde. (not ordered)
    //to be able to have a ordered output of our state we need to know the position every state variable should be written to
    //order is given by (n+1)*y+x, where n is the number of degrees of freedom per column minus one

    order.resize(n_y);
    for (int i = 0; i < n_y; ++i) {
	order[i] = dof_x[i] + discretization_n * dof_y[i];
    }
}

/*
 * Multiply the matrix given by rows, cols, vals with a vector vec
 * Used in heat_nlp.cpp
 */
valarray<double> MATRIXOP::matrix_vector_mult(valarray<int> &rows,
	valarray<int> &cols, valarray<double> &vals, valarray<double> &vec, unsigned length) {
    long size_matrix = rows.size();
    valarray<double> re(length);
    int iter = 0;
    for (int i = 0; i < length; ++i) {
	double sum = 0;
	while (iter < size_matrix && rows[iter] == i) {
	    sum += vals[iter] * vec[cols[iter]];
	    ++iter;
	}
	re[i] = sum;
    }
    return re;
}
