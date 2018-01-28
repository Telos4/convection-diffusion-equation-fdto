/*
 * hs071_nlp.hpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#ifndef MATRIXOP_HPP_
#define MATRIXOP_HPP_

#include <valarray>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>


using namespace std;

class MATRIXOP {
public:

    MATRIXOP();
    MATRIXOP(int N_, string file_A, string file_B, string file_b_u, string file_b_y,
            double eps_, double y_ref_, double u_ref_, bool closed_values, bool open_values);


    /** Default destructor */
    ~MATRIXOP();

    void read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals);
    void read_vector(string filename, valarray<double> &vals);
    void print_matrix(valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) const;
    void print_vector(valarray<double> &vals) const;

    double eval_f(valarray<double> &x);
    double vec_Q_vec(valarray<double> y, double y_ref);
    double vec_R_vec(valarray<double> y, double y_ref);

    valarray<double> eval_grad_f(valarray<double> z);
    valarray<double> Q_vec(valarray<double> y);
    valarray<double> R_vec(valarray<double> u);

    void A_eq();
    valarray<double> matrix_vektor_mult(valarray<int> &rows,
        valarray<int> &cols, valarray<double> &vals, valarray<double> &vec);
    
    
    
    
    valarray<int> A_rows, A_cols, B_rows, B_cols, A_eq_rows, A_eq_cols;
    valarray<double> A_vals, B_vals, A_eq_vals, b_u, b_y, y_old, u_old;
    int n_y, n_u, n_z, N, iter;
    double eps, y_ref, u_ref, closed_loop_cost;
    
    bool closed_values, open_values;
    //string closed_file, open_file;
private:


};





#endif /* MATRIXOP_HPP_ */
