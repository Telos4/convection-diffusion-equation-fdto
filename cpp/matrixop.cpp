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

MATRIXOP::MATRIXOP(int N_, string file_A, string file_B, string file_b_u, string file_b_y) {
    read_matrix(file_A, A_rows, A_cols, A_vals);
    read_matrix(file_B, B_rows, B_cols, B_vals);
    read_vector(file_b_u, b_u);
    read_vector(file_b_y, b_y);

    N = N_;
    n_y = A_rows.size();

    //check if N and n_u are of compatible size
    if(b_u.size()%N == 0){
    n_u = b_u.size() / N;
    }
    else{
        cout << "dimensions are wrong" << endl;
        exit(1);
    }
    n_z = (N + 1) * n_y + N*n_u;

}

//destructor

MATRIXOP::~MATRIXOP() {
}

//todo skip comments in .mtx files
void MATRIXOP::read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) {
    //ifstream input(filename);
    ifstream ifs(filename.c_str(), ifstream::in);

    if (ifs.is_open()) {
        int n, i = 0;
        ifs >> n >> n >> n;
        rows.resize(n);
        cols.resize(n);
        vals.resize(n);

        while (ifs >> rows[i] >> cols[i] >> vals[i]) {
            ++i;
        }
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

//todo insert matrices Q, R
double MATRIXOP::eval_f(double* x, int n, double eps, double y_ref, double u_ref){

}
