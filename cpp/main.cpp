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

int optimize(int steps, MATRIXOP & data, int outputlevel);
void closed_cost_vs_horizon_length(int max_MPC_horizon, int steps, string matrix_A, string matrix_B_y,
	string matrix_B_w, string vec_b_u, string vec_b_y_out, string dof_x, string dof_y, string pythonparam, double eps, double y_ref,
	double u_ref, double u_upper, double u_lower, double y_upper, double y_lower, double w_upper, double w_lower, bool dim2, bool convection, int outputlevel);
//void solver_test(MATRIXOP &data, int steps);
void write_parameters(int MPC_horizon, int steps, string matrix_A, string matrix_B_y, string matrix_B_w, string vec_b_u,
	string vec_b_y_out, string dof_x, string dof_y, string pythonparam_, double eps, double y_ref, double u_ref,
	double u_upper, double u_lower, double y_upper, double y_lower, double w_upper, double w_lower, bool dim2,
	bool convection, bool closed_values, bool open_values, bool free_init_value, bool cost_vs_horizon,
	string foldername);

int main(int argc, char** argv) {
    bool dim2 = false;
    bool convection = false;
    bool closed_values = false;
    bool open_values = false;
    bool cost_vs_horizon = false;
    bool free_init_value = false;

    double eps = 10e-3;
    double y_ref = 0.0;
    double u_ref = 0.0;
    double u_upper = 0.25;
    double u_lower = -0.25;
    double y_upper = 0.15;
    double y_lower = -0.15;
    double w_upper = 1e19;
    double w_lower = -1e19;

    string matrix_A = "../A.mtx";
    string matrix_B_y = "../B_y.mtx";
    string matrix_B_w = "../B_w.mtx";
    string vec_b_u = "../b_u.txt";
    string vec_b_y_out = "../b_y_out.txt";
    string dof_x = "../dof_x.txt";
    string dof_y = "../dof_y.txt";
    string result_folder = "../results";
    string result_folder_prefix = "";
    string pythonparam = "../python_parameters.txt";


    int steps = 1;
    int MPC_horizon = 100;

    int outputlevel = 0;

    args::ArgumentParser parser("convection diffusion equation 1d.", "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu",{'h', "help"});
    args::CompletionFlag completion(parser,{"complete"});

    args::Flag dim2_(parser, "2 dimensional", "turn dimension 2 on",{'d', "dim2"});
    args::Flag convection_(parser, "convection", "turn convection on",{'c', "convection"});
    args::Flag closed_values_(parser, "closed values", "save states of the closed loop",{"cv", "closedvalues"});
    args::Flag open_values_(parser, "open values", "save states of all open loop",{"ov", "openvalues"});
    args::Flag free_init_value_(parser, "free initial value", "indicates whether the intial state should be free",{"fi", "free-init-value"});
    args::Flag cost_vs_horizon_(parser, "cost vs mpc horizon plot",
	    "create data for cost vs mpc horizon plot, closed and open loop values will not be written, MPC_horizon functions as your maximal MPC horizon",{"cost_vs_horizon"});
    args::ValueFlag<int> MPC_horizon_(parser, "mpc horizon", "MPC horizon, N >= 1",{'N', "mpc"});
    args::ValueFlag<int> steps_(parser, "steps", "closed loop step count",{'L', "steps"});
    args::ValueFlag<int> outputlevel_(parser, "outputlevel", "outputlevel of IPopt",{"output"});
    args::ValueFlag<double> eps_(parser, "epsilon", "epsilon from cost functional",{'e', "eps"});
    args::ValueFlag<double> u_ref_(parser, "u_ref", "reference solution for control",{'u', "uref"});
    args::ValueFlag<double> y_ref_(parser, "y_ref", "reference solution for state",{'y', "yref"});
    args::ValueFlag<double> u_upper_(parser, "u_upper", "upper bound for control",{"u_upper"});
    args::ValueFlag<double> u_lower_(parser, "u_lower", "lower bound for control",{"u_lower"});
    args::ValueFlag<double> y_upper_(parser, "y_upper", "upper bound for state",{"y_upper"});
    args::ValueFlag<double> y_lower_(parser, "y_lower", "lower bound for state",{"y_lower"});
    args::ValueFlag<double> w_upper_(parser, "w_upper", "upper bound for convection",{"w_upper"});
    args::ValueFlag<double> w_lower_(parser, "w_lower", "lower bound for convection",{"w_lower"});
    args::ValueFlag<string> matrix_A_(parser, "Matrix A", "matrix A from PDE",{"matA"});
    args::ValueFlag<string> matrix_B_y_(parser, "Matrix B_y", "matrix B_y from PDE",{"matB_y"});
    args::ValueFlag<string> matrix_B_w_(parser, "Matrix B_w", "matrix B_w from PDE",{"matB_w"});
    args::ValueFlag<string> vec_b_u_(parser, "vec_b_u", "vector b_u from PDE",{"b_u"});
    args::ValueFlag<string> vec_b_y_out_(parser, "vec_b_y_out", "vector b_y_out from PDE",{"b_y_out"});
    args::ValueFlag<string> dof_x_(parser, "dof_x", "vector dof_x from PDE",{"dof_x"});
    args::ValueFlag<string> dof_y_(parser, "dof_y", "vector dof_y from PDE",{"dof_y"});
    args::ValueFlag<string> result_folder_(parser, "result_folder", "folder where results shall be stored",{"result_folder"});
    args::ValueFlag<string> result_folder_prefix_(parser, "result_folder_prefix", "prefix for created folder",{"result_folder_prefix"});
    args::ValueFlag<string> pythonparam_(parser, "pythonparam", "parameter file from PDE",{"pythonparam"});


    try {
	parser.ParseCLI(argc, argv);
    }
    catch (args::Completion e) {
	std::cout << e.what();
	return 0;
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

    dim2 = dim2_;
    convection = convection_;
    closed_values = closed_values_;
    open_values = open_values_;
    cost_vs_horizon = cost_vs_horizon_;
    free_init_value = free_init_value_;


    if (MPC_horizon_) {
	MPC_horizon = args::get(MPC_horizon_);
    }
    if (steps_) {
	steps = args::get(steps_);
    }
    if (outputlevel_) {
	outputlevel = args::get(outputlevel_);
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
    if (u_upper_) {
	u_upper = args::get(u_upper_);
    }
    if (u_lower_) {
	u_lower = args::get(u_lower_);
    }
    if (y_upper_) {
	y_upper = args::get(y_upper_);
    }
    if (y_lower_) {
	y_lower = args::get(y_lower_);
    }
    if (w_upper_) {
	w_upper = args::get(w_upper_);
    }
    if (w_lower_) {
	w_lower = args::get(w_lower_);
    }
    if (matrix_A_) {
	matrix_A = args::get(matrix_A_);
    }
    if (matrix_B_y_) {
	matrix_B_y = args::get(matrix_B_y_);
    }
    if (matrix_B_w_) {
	matrix_B_w = args::get(matrix_B_w_);
    }
    if (vec_b_u_) {
	vec_b_u = args::get(vec_b_u_);
    }
    if (vec_b_y_out_) {
	vec_b_y_out = args::get(vec_b_y_out_);
    }
    if (dof_x_) {
	dof_x = args::get(dof_x_);
    }
    if (dof_y_) {
	dof_y = args::get(dof_y_);
    }
    if (result_folder_) {
	result_folder = args::get(result_folder_);
    }
    if (result_folder_prefix_) {
	result_folder_prefix = args::get(result_folder_prefix_);
    }
    if (pythonparam_) {
	pythonparam = args::get(pythonparam_);
    }




    //start program
    if (cost_vs_horizon) {
	closed_cost_vs_horizon_length(MPC_horizon, steps, matrix_A, matrix_B_y, matrix_B_w, vec_b_u, vec_b_y_out, dof_x, dof_y, pythonparam, eps, y_ref, u_ref, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, 
		dim2, convection, outputlevel);
    }
    else {
	MATRIXOP data(MPC_horizon, matrix_A, matrix_B_y, matrix_B_w, vec_b_u, vec_b_y_out, dof_x, dof_y, eps, y_ref, u_ref, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, 
		dim2, convection, closed_values, open_values, free_init_value, result_folder, result_folder_prefix);
	int status = optimize(steps, data, outputlevel);

	if (status == 0) {
	    cout << "problem solved \n";
	}
	else {
	    cout << "IPopt could'nt solve the problem, use higher outputlevel to investigate \n";
	    return (1);
	}


	if (closed_values || open_values) {
	    write_parameters(MPC_horizon, steps, matrix_A, matrix_B_y, matrix_B_w, vec_b_u, vec_b_y_out, dof_x, dof_y, pythonparam, eps,
		    y_ref, u_ref, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, dim2, convection, closed_values, open_values, free_init_value, cost_vs_horizon,
		    data.foldername);
	}
    }

    return 0;
}

int optimize(int steps, MATRIXOP & data, int outputlevel) {
    int last_status;
    //outer loop
    for (int i = 0; i < steps; ++i) {
	SmartPtr<TNLP> mynlp = new HEAT_NLP(data);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	// Change some options
	//app->Options()->SetStringValue("derivative_test", "second-order");
	//app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetIntegerValue("print_level", outputlevel);

	app->Options()->SetNumericValue("tol", 1e-5);
	app->Options()->SetStringValue("jac_c_constant", "yes");
	app->Options()->SetStringValue("jac_d_constant", "yes");
	app->Options()->SetStringValue("hessian_constant", "yes");
	app->Options()->SetStringValue("linear_solver", "ma27");

	if (data.convection) {
	    app->Options()->SetStringValue("jac_c_constant", "no");
	    app->Options()->SetStringValue("jac_d_constant", "no");

	    //approximation of the hessian of the lagranage function by only using the hessian of the objective function
	    //error of ~10e-3
	    app->Options()->SetStringValue("hessian_constant", "yes");
	    //app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	}

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
	    //return (int) status;
	}
	else {
	    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
	    return (int) status;
	}

	last_status = (int) status;
    }
    return last_status;
}

void closed_cost_vs_horizon_length(int max_MPC_horizon, int steps, string matrix_A, string matrix_B_y,
	string matrix_B_w, string vec_b_u, string vec_b_y_out, string dof_x, string dof_y, string pythonparam,
	double eps, double y_ref, double u_ref, double u_upper, double u_lower, double y_upper, double w_upper, double w_lower, double y_lower, bool dim2, bool convection, int outputlevel) {
    int status;
    bool closed_values = false;
    bool open_values = false;
    bool free_init_value = false;
    valarray<double> costs(0.0, max_MPC_horizon);
    string foldername;
    string result_folder = "../results";
    string result_folder_prefix = "";


    //start at 1, 0 does not make sense
    for (int i = 1; i < max_MPC_horizon; ++i) {
	MATRIXOP data(i, matrix_A, matrix_B_y, matrix_B_w, vec_b_u, vec_b_y_out, dof_x, dof_y, eps, y_ref, u_ref, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, dim2,
		convection, closed_values, open_values, free_init_value, result_folder, result_folder_prefix);

	if (i == 1) {
	    data.create_folder("");
	    foldername = data.foldername;
	}

	status = optimize(steps, data, outputlevel);
	if (status == 0) {
	    std::cout << "Step: " << i << " *** success!" << std::endl;
	    costs[i] = data.closed_loop_cost;
	}
	    //set costs to 0 if no solution was found
	else {
	    costs[i] = 0;
	    std::cout << "Step: " << i << " *** failure!" << std::endl;
	}
    }

    ofstream out;
    out.open(foldername + "closed_cost_vs_horizon_length.txt", ofstream::app);
    if (out.is_open()) {
	for (int i = 0; i < max_MPC_horizon; ++i) {
	    out << costs[i] << endl;
	}
    }
    else {
	cout << "couldnt open closed_cost_vs_horizon_length.txt ... exiting" << endl;
	exit(1);
    }

    write_parameters(max_MPC_horizon, steps, matrix_A, matrix_B_y, matrix_B_w, vec_b_u, vec_b_y_out, dof_x, dof_y, pythonparam, eps,
	    y_ref, u_ref, u_upper, u_lower, y_upper, y_lower, w_upper, w_lower, dim2, convection, closed_values, open_values, free_init_value, true, foldername);

    out.close();
}

void write_parameters(int MPC_horizon, int steps, string matrix_A, string matrix_B_y, string matrix_B_w, string vec_b_u, string vec_b_y_out, string dof_x, string dof_y, string pythonparam_,
	double eps, double y_ref, double u_ref, double u_upper, double u_lower, double y_upper, double y_lower, double w_upper, double w_lower, bool dim2, bool convection, bool closed_values, bool open_values, bool free_init_value, bool cost_vs_horizon, string foldername) {
    double alpha = 0, beta = 0, gamma = 0;
    int n = 0;
    string trash;
    ifstream pythonparam(pythonparam_);
    if (pythonparam.is_open()) {
	pythonparam >> n >> trash >> alpha >> trash >> beta >> trash >> gamma;
	pythonparam.close();
    }
    else {
	cout << "can't open parameter file from python script" << endl;
	exit(1);
    }

    ofstream out(foldername + "parameters.txt");
    if (out.is_open()) {

	string function = "0.3 * sin(0.1 * i)";

	out << "cost_vs_horizon \t" << cost_vs_horizon << endl;
	out << "discretization_parameter \t" << n << endl;
	out << "alpha \t" << alpha << endl;
	out << "beta \t" << beta << endl;
	out << "gamma \t" << gamma << endl;
	out << "MPC_horizon \t" << MPC_horizon << endl;
	out << "steps \t" << steps << endl;
	out << "eps \t" << eps << endl;
	out << "y_ref \t" << y_ref << endl;
	out << "u_ref \t" << u_ref << endl;
	out << "u_upper \t" << u_upper << endl;
	out << "u_lower \t" << u_lower << endl;
	out << "y_upper \t" << y_upper << endl;
	out << "y_lower \t" << y_lower << endl;
	out << "w_upper \t" << w_upper << endl;
	out << "w_lower \t" << w_lower << endl;
	out << "dim2 \t" << dim2 << endl;
	out << "convection \t" << convection << endl;
	out << "closed_values \t" << closed_values << endl;
	out << "open_values \t" << open_values << endl;
	out << "free_init_value \t" << free_init_value << endl;
	out << "function \t" << function << endl;
	out << "matrix_A \t" << matrix_A << endl;
	out << "matrix_B_y \t" << matrix_B_y << endl;
	out << "matrix_B_w \t" << matrix_B_w << endl;
	out << "vec_b_u \t" << vec_b_u << endl;
	out << "vec_b_y_out \t" << vec_b_y_out << endl;
	out << "dof_x \t" << dof_x << endl;
	out << "dof_y \t" << dof_y << endl;

	out.close();
    }
    else {
	cout << "can't open paramater file to write to" << endl;
	exit(1);
    }
}



/*
void solver_test(MATRIXOP &data, int steps) {
    vector<string> solver = {"ma27", "ma57", "ma77", "ma86", "ma97", "mumps"};

    for (unsigned i = 0; i < solver.size(); ++i) {
	clock_t t;
	t = clock();
	optimize(steps, data, solver[i]);
	t = clock() - t;

	cout << solver[i] << ": " << t << endl;
    }
}
 */

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
