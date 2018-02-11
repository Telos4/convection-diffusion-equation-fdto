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


using namespace std;

class MATRIXOP {
public:

    MATRIXOP();
    MATRIXOP(int N_, string file_A, string file_B_y, string file_B_w, string file_b_u, string file_b_y,
            double eps_, double y_ref_, double u_ref_, bool convection, bool closed_values, bool open_values);


    /** Default destructor */
    ~MATRIXOP();

    void read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals);
    void read_vector(string filename, valarray<double> &vals);
    void print_matrix(valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) const;
    void print_vector(valarray<double> &vals) const;
    void create_folder();

    double eval_f(valarray<double> &x);
    double vec_Q_vec(valarray<double> y, double y_ref);
    double vec_R_vec(valarray<double> y, double y_ref);
    double vec_W_vec(valarray<double> y);

    valarray<double> eval_grad_f(valarray<double> z);
    valarray<double> Q_vec(valarray<double> y);
    valarray<double> R_vec(valarray<double> u);
    valarray<double> W_vec(valarray<double> w);

    void A_eq();
    valarray<double> matrix_vector_mult(valarray<int> &rows,
        valarray<int> &cols, valarray<double> &vals, valarray<double> &vec);
    
    
    
    
    valarray<int> A_rows, A_cols, B_y_rows, B_y_cols, A_eq_rows, A_eq_cols, B_w_rows, B_w_cols;
    valarray<double> A_vals, B_y_vals, A_eq_vals, B_w_vals, b_u, b_y, y_old, u_old, w_old;
    int n_y, n_u, n_w, n_z, N, iter;
    double eps, y_ref, u_ref, closed_loop_cost;
    
    bool convection, closed_values, open_values;
    string foldername;
private:


};





#endif /* MATRIXOP_HPP_ */
