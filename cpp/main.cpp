/*
 * main.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include <string>
#include <vector>
#include <time.h> 

#include "heat_nlp.hpp"
#include "IpIpoptApplication.hpp"
#include "args.hxx"


using namespace Ipopt;

int optimize(int steps, MATRIXOP & data);
int optimize(int steps, MATRIXOP & data, double & closed_loop_cost);
int optimize(int steps, MATRIXOP & data, string solver);
void closed_cost_vs_horizon_length(int max_MPC_horizon, int steps, string matrix_A, string matrix_B, string vec_b_u, string vec_b_y_out, double eps, double y_ref, double u_ref);
void closed_cost_vs_steps(int steps, int MPC_horizon, string matrix_A, string matrix_B, string vec_b_u, string vec_b_y_out, double eps, double y_ref, double u_ref);
void solver_test(MATRIXOP &data, int steps);

int main(int argc, char** argv) {

    bool closed_values = false;
    bool open_values = false;

    double eps = 10e-3;
    double y_ref = 0.5;
    double u_ref = 0.5;
    string matrix_A = "../A.mtx";
    string matrix_B = "../B_y.mtx";
    string vec_b_u = "../b_u.txt";
    string vec_b_y_out = "../b_y_out.txt";
    
    int steps = 200;
    int MPC_horizon = 10;

    int max_MPC_horizon = 10;


    args::ArgumentParser parser("This is a test program.", "This goes after the options.");
    
    args::Flag closed_values_(parser, "closed values", "save states of the closed loop",{'c', "closedvalues"});
    args::Flag open_values_(parser, "open values", "save states of all open loop",{'o', "openvalues"});

    args::ValueFlag<int> MPC_horizon_(parser, "mpc horizon", "MPC horizon",{'N', "mpc"});
    args::ValueFlag<int> steps_(parser, "steps", "closed loop step count",{'L', "steps"});
    args::ValueFlag<double> eps_(parser, "epsilon", "epsilon from cost functional",{'e', "eps"});
    args::ValueFlag<double> u_ref_(parser, "u_ref", "reference solution for control",{'u', "uref"});
    args::ValueFlag<double> y_ref_(parser, "y_ref", "reference solution for state",{'y', "yref"});
    args::ValueFlag<string> matrix_A_(parser, "Matrix A", "matrix A from PDE",{"matA"});
    args::ValueFlag<string> matrix_B_(parser, "Matrix B", "matrix B from PDE",{"matB"});
    args::ValueFlag<string> vec_b_u_(parser, "vec_b_u", "vector b_u from PDE",{"b_u"});
    args::ValueFlag<string> vec_b_y_out_(parser, "vec_b_y_out", "vector b_y_out from PDE",{"b_y_out"});

    try {
	parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
	std::cout << parser;
	return 0;
    }
    catch (args::ParseError e) {
	std::cerr << e.what() << std::endl;
	std::cerr << parser;
	return 1;
    }
    catch (args::ValidationError e) {
	std::cerr << e.what() << std::endl;
	std::cerr << parser;
	return 1;
    }

    closed_values = closed_values_;
    open_values = open_values_;

    if (MPC_horizon_) {
	MPC_horizon = args::get(MPC_horizon_);
    }
    if (steps_) {
	steps = args::get(steps_);
    }
    if (eps_) {
	eps = args::get(eps_);
    }
    if (u_ref_) {
	u_ref = args::get(u_ref_);
    }
    if (y_ref_) {
	y_ref = args::get(y_ref_);
    }
    if (matrix_A_) {
	matrix_A = args::get(matrix_A_);
    }
    if (matrix_B_) {
	matrix_B = args::get(matrix_B_);
    }
    if (vec_b_u_) {
	vec_b_u = args::get(vec_b_u_);
    }
    if (vec_b_y_out_) {
	vec_b_y_out = args::get(vec_b_y_out_);
    }

    //closed_cost_vs_horizon_length(max_MPC_horizon, steps, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref);

    //closed_cost_vs_steps(steps, MPC_horizon, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref);

    //MATRIXOP data(MPC_horizon, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref);
    //optimize(steps, data);



    MATRIXOP data(MPC_horizon, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref, closed_values, open_values);
    optimize(steps, data);


    return 0;
}

int optimize(int steps, MATRIXOP & data) {
    int last_status;
    //outer loop
    for (int i = 0; i < steps; ++i) {
	SmartPtr<TNLP> mynlp = new HEAT_NLP(data);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	//app->Options()->SetStringValue("derivative_test", "second-order");
	//app->Options()->SetNumericValue("max_iter", 100);
	app->Options()->SetNumericValue("tol", 1e-5);
	//app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetIntegerValue("print_level", 1);
	app->Options()->SetStringValue("jac_c_constant", "yes");
	app->Options()->SetStringValue("jac_d_constant", "yes");
	app->Options()->SetStringValue("hessian_constant", "yes");
	app->Options()->SetStringValue("linear_solver", "ma27");

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
	    //return status;

	}
	else {
	    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	    return (int) status;
	}

	last_status = (int) status;
    }
    return last_status;
}

int optimize(int steps, MATRIXOP & data, double & closed_loop_cost) {
    int last_status;
    //outer loop
    for (int i = 0; i < steps; ++i) {
	SmartPtr<TNLP> mynlp = new HEAT_NLP(data);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	//app->Options()->SetStringValue("derivative_test", "second-order");
	//app->Options()->SetNumericValue("max_iter", 100);
	app->Options()->SetNumericValue("tol", 1e-5);
	//app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetIntegerValue("print_level", 0);
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
	    //std::cout << "Step: " << i << " *** The problem solved!" << std::endl;
	    closed_loop_cost = data.closed_loop_cost;
	    //return status;

	}
	else {
	    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	    return (int) status;
	}

	last_status = (int) status;
    }
    return last_status;
}

int optimize(int steps, MATRIXOP & data, string solver) {
    int last_status;
    //outer loop
    for (int i = 0; i < steps; ++i) {
	SmartPtr<TNLP> mynlp = new HEAT_NLP(data);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	//app->Options()->SetStringValue("derivative_test", "second-order");
	//app->Options()->SetNumericValue("max_iter", 100);
	app->Options()->SetNumericValue("tol", 1e-5);
	//app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("jac_c_constant", "yes");
	app->Options()->SetStringValue("jac_d_constant", "yes");
	app->Options()->SetStringValue("hessian_constant", "yes");
	app->Options()->SetStringValue("linear_solver", solver);

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
	    //std::cout << "Step: " << i << " *** The problem solved!" << std::endl;
	    //return status;

	}
	else {
	    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	    return (int) status;
	}

	last_status = (int) status;
    }
    return last_status;
}

void closed_cost_vs_horizon_length(int max_MPC_horizon, int steps, string matrix_A, string matrix_B, string vec_b_u, string vec_b_y_out, double eps, double y_ref, double u_ref) {
    int status;
    bool closed_values = false;
    bool open_values = false;
    valarray<double> costs(max_MPC_horizon);
    for (int i = 0; i < max_MPC_horizon; ++i) {
	MATRIXOP data(i, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref, closed_values, open_values);
	status = optimize(steps, data, costs[i]);
	if (status == 0) {
	    std::cout << "Step: " << i << " *** success!" << std::endl;
	}
	else {
	    costs[i] = 0;
	    std::cout << "Step: " << i << " *** failure!" << std::endl;
	}
    }

    ofstream out("../results/closed_cost_vs_horizon_length.txt");
    for (int i = 0; i < max_MPC_horizon; ++i) {
	out << costs[i] << endl;
    }
    out.close();
}

void closed_cost_vs_steps(int steps, int MPC_horizon, string matrix_A, string matrix_B, string vec_b_u, string vec_b_y_out, double eps, double y_ref, double u_ref) {
    bool closed_values = false;
    bool open_values = false;
    
    valarray<double> costs(steps);
    MATRIXOP data(MPC_horizon, matrix_A, matrix_B, vec_b_u, vec_b_y_out, eps, y_ref, u_ref, closed_values, open_values);
    for (int i = 0; i < steps; ++i) {
	SmartPtr<TNLP> mynlp = new HEAT_NLP(data);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	app->Options()->SetNumericValue("tol", 1e-5);
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("jac_c_constant", "yes");
	app->Options()->SetStringValue("jac_d_constant", "yes");
	app->Options()->SetStringValue("hessian_constant", "yes");

	// Initialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
	    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
	}
	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

	if (status == Solve_Succeeded) {
	    std::cout << "Step: " << i << " *** The problem solved!" << std::endl;
	    costs[i] = data.closed_loop_cost;
	}
	else {
	    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	}
    }

    ofstream out("../results/closed_cost_vs_steps_10.txt");
    for (int i = 0; i < steps; ++i) {
	out << costs[i] << endl;
    }
    out.close();
}

void solver_test(MATRIXOP &data, int steps){
        vector<string> solver = {"ma27", "ma57", "ma77", "ma86", "ma97", "mumps"};
    
    for (unsigned i = 0; i < solver.size(); ++i) {
	clock_t t;
	t = clock();
	optimize(steps, data, solver[i]);
	t = clock() - t;

	cout << solver[i] << ": " << t << endl;
    }
}












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

}
 */