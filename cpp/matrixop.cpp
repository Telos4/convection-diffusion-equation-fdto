/*
 * matrixop.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include "matrixop.hpp"

using namespace std;

// constructor

MATRIXOP::MATRIXOP() {
}

MATRIXOP::MATRIXOP(int N_, string file_A, string file_B_y, string file_B_w, string file_b_u, string file_b_y, string file_dof_x, string file_dof_y, 
	double eps_, double y_ref_, double u_ref_, bool dim2_, bool convection_, bool closed_values_, bool open_values_) {
    
    dim2 = dim2_;
    convection = convection_;

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

    if (N_ == 0) {
	cout << "N = 0 does not make sense, exiting" << endl;
	exit(1);
    }
    N = N_;

    n_y = b_u.size();
    n_u = 1;

    if (convection) {
	n_w = 1;
    }
    else {
	n_w = 0;
    }

    n_z = (N + 1) * n_y + N * n_u + N * n_w;

    eps = eps_;
    y_ref = y_ref_;
    u_ref = u_ref_;

    iter = 0;
    closed_loop_cost = 0;

    closed_values = closed_values_;
    open_values = open_values_;

    //initialize A_eq
    A_eq();

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

    if (closed_values || open_values) {
	create_folder();
    }
}

//destructor

MATRIXOP::~MATRIXOP() {
}

void MATRIXOP::read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) {
    ifstream ifs(filename.c_str(), ifstream::in);

    if (ifs.is_open()) {
	//skip comments in .mtx files 
	string trash;

	getline(ifs, trash);
	getline(ifs, trash);

	int n, i = 0;
	//third value in file is size
	ifs >> n >> n >> n;
	rows.resize(n);
	cols.resize(n);
	vals.resize(n);

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
	cout << endl << "File " << filename << " is not open" << endl;
    }
}

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
	cout << endl << "File " << filename << " is not open" << endl;
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

void MATRIXOP::create_folder() {

    time_t t = time(0); // get time now
    struct tm * now = localtime(& t);

    foldername = "../results/" + to_string(now->tm_year + 1900) + "-" + to_string(now->tm_mon + 1) + "-" + to_string(now->tm_mday)
	    + "_" + to_string(now->tm_hour) + "-" + to_string(now->tm_min) + "-" + to_string(now->tm_sec) + "/";

    //cout << boost::filesystem::current_path().string() << endl;
    boost::filesystem::path p(foldername);

    //create results folder first, cant create subdirectories directly?
    try {
	boost::filesystem::create_directory("../results");
    }
    catch (boost::filesystem::filesystem_error &e) {
	std::cerr << e.what() << '\n';
    }
    try {
	boost::filesystem::create_directory(p);
    }
    catch (boost::filesystem::filesystem_error &e) {
	std::cerr << e.what() << '\n';
    }
}


//todo insert matrices Q, R, ...  y_ref, u_ref as vector?

double MATRIXOP::eval_f(valarray<double> &x) {
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

//Q, R schwach besetzt?, dann später zu schreibende matrix vektor multiplikation nutzen

valarray<double> MATRIXOP::eval_grad_f(valarray<double> z) {
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

void MATRIXOP::A_eq() {
    long size_A_eq = (B_y_rows.size() + A_rows.size() + b_u.size()) * N;
    long size_B = B_y_rows.size();
    long size_A = A_rows.size();
    A_eq_rows.resize(size_A_eq);
    A_eq_cols.resize(size_A_eq);
    A_eq_vals.resize(size_A_eq);

    int count = 0;
    for (int i = 0; i < N; ++i) {
	int countA = 0;
	int countB = 0;
	for (int k = 0; k < n_y; ++k) {
	    while (countB != size_B && B_y_rows[countB] == k) {
		A_eq_rows[count] = i * n_y + k;
		A_eq_cols[count] = i * n_y + B_y_cols[countB];
		A_eq_vals[count] = -1 * B_y_vals[countB];
		++count;
		++countB;
	    }

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

valarray<double> MATRIXOP::matrix_vector_mult(valarray<int> &rows,
	valarray<int> &cols, valarray<double> &vals, valarray<double> &vec) {
    long size_matrix = rows.size();
    valarray<double> re(N * n_y);
    int iter = 0;
    for (int i = 0; i < N * n_y; ++i) {
	double sum = 0;
	while (iter < size_matrix && rows[iter] == i) {
	    sum += vals[iter] * vec[cols[iter]];
	    ++iter;
	}
	re[i] = sum;
    }
    return re;
}