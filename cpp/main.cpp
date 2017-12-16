/*
 * main.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */
//Zum testen von matrixop

//#include "IpIpoptApplication.hpp"
#include "heat_nlp.hpp"
#include "matrixop.hpp"

using namespace Ipopt;

int main(int argv, char* argc[])
{
	valarray<int> a,b;
	valarray<double> c, d;
	MATRIXOP m;
	m.read_matrix("A_test.mtx", a, b, c);
	m.print_matrix(a ,b, c);

	m.read_vector("testvector", d);
	m.print_vector(d);

        MATRIXOP first(101, "A.mtx", "B_y.mtx", "b_u.txt", "b_y_out.txt");
        first.print_matrix(first.A_rows, first.A_cols, first.A_vals);
        first.print_matrix(first.B_rows, first.B_cols, first.B_vals);
        first.print_vector(first.b_u);
        first.print_vector(first.b_y);

	return 0;
}
