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

class MATRIXOP
{
public:

    MATRIXOP();
    MATRIXOP(int N_, string file_A, string file_B, string file_b_u, string file_b_y);


  /** Default destructor */
    ~MATRIXOP();

    void read_matrix(string filename, valarray<int> &rows, valarray<int> &cols, valarray<double> &vals);
    void read_vector(string filename, valarray<double> &vals);
    void print_matrix(valarray<int> &rows, valarray<int> &cols, valarray<double> &vals) const;
    void print_vector(valarray<double> &vals) const;



    valarray<int> A_rows, A_cols, B_rows, B_cols;
    valarray<double> A_vals, B_vals, b_u, b_y;
    int n_y, n_u, n_z, N;
private:


};





#endif /* MATRIXOP_HPP_ */
