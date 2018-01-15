/*
 * main.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

//#include "IpIpoptApplication.hpp"
#include "heat_nlp.hpp"
//#include "matrixop.hpp"
#include "IpIpoptApplication.hpp"


using namespace Ipopt;

int main(int argv, char* argc[]) {

    int steps = 5;
    MATRIXOP data1(10, "A.mtx", "B_y.mtx", "b_u.txt", "b_y_out.txt", 10e-3, 0.5, 0.5);
    //MATRIXOP data2(1, "A_small.mtx", "B_y_small.mtx", "b_u_small.txt", "b_y_out_small.txt", 10e-3, 0.5, 0.5);

    //outer loop
    for (int i = 0; i < steps; ++i) {
        SmartPtr<TNLP> mynlp = new HEAT_NLP(data1);

        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
        app->RethrowNonIpoptException(true);

        // Change some options
        //app->Options()->SetStringValue("derivative_test", "second-order");
        //app->Options()->SetNumericValue("max_iter", 100);
        app->Options()->SetNumericValue("tol", 1e-5);
        //app->Options()->SetStringValue("output_file", "ipopt.out");
        app->Options()->SetNumericValue("print_level", (int)0);
        app->Options()->SetStringValue("jac_c_constant", "yes");
        app->Options()->SetStringValue("jac_d_constant", "yes");
        app->Options()->SetStringValue("hessian_constant", "yes");

        // The following overwrites the default name (ipopt.opt) of the
        // options file
        // app->Options()->SetStringValue("option_file_name", "hs071.opt");

        // Initialize the IpoptApplication and process the options
        ApplicationReturnStatus status;
        status = app->Initialize();
        if (status != Solve_Succeeded) {
            std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
            return (int) status;
        }

        // Ask Ipopt to solve the problem
        status = app->OptimizeTNLP(mynlp);

        if (status == Solve_Succeeded) {
            std::cout << "Step: " << i << " *** The problem solved!" << std::endl;
        }
        else {
            std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
            return (int) status;
        }
    }

    cout << endl << endl << "maybe this works" << endl;
    return 0;








    //function tests
    /*
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
     */
    /*
    MATRIXOP second(2, "A_small.mtx", "B_y_small.mtx", "b_u_small.txt", "b_y_out_small.txt", 1, 0, 0);
    valarray<double> ones(second.n_z);
    for(int i = 0; i < second.n_z; ++i){
        ones[i] = 1;
    }
    valarray<double> eval = second.eval_grad_f(ones);

    second.print_vector(eval);

    cout << "z: " << second.n_z << " n_u: " << second.n_u << " n_y: " << second.n_y << endl;

    double f = second.eval_f(ones);
    cout << "eval f: " << f << endl;
     */
    /*
       MATRIXOP third(2, "A_small.mtx", "B_y_small.mtx", "b_u_small.txt", "b_y_out_small.txt", 1, 0, 0);
       third.print_matrix(third.A_eq_rows, third.A_eq_cols, third.A_eq_vals);
     */
}
