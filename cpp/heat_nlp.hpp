/*
 * hs071_nlp.hpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#ifndef HEAT_NLP_HPP_
#define HEAT_NLP_HPP_

#include "IpTNLP.hpp"
#include "matrixop.hpp"
#include <valarray>
#include <fstream>

using namespace Ipopt;

class HEAT_NLP : public TNLP {
public:

    //HEAT_NLP();
    //HEAT_NLP(int N_, string file_A, string file_B, string file_b_u, string file_b_y,
    //    double eps, double y_ref, double u_ref);
    HEAT_NLP(MATRIXOP &data_);    

    /** Default destructor */
    virtual ~HEAT_NLP();


    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
            Index& nnz_h_lag, IndexStyleEnum& index_style);


    /** overload this method to return the information about the bound
     *  on the variables and constraints. The value that indicates
     *  that a bound does not exist is specified in the parameters
     *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
     *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
     *  1e19. (see TNLPAdapter) */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
            Index m, Number* g_l, Number* g_u);


    /** overload this method to return the starting point. The bool
     *  variables indicate whether the algorithm wants you to
     *  initialize x, z_L/z_u, and lambda, respectively.  If, for some
     *  reason, the algorithm wants you to initialize these and you
     *  cannot, return false, which will cause Ipopt to stop.  You
     *  will have to run Ipopt with different options then.
     */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
            bool init_z, Number* z_L, Number* z_U,
            Index m, bool init_lambda,
            Number* lambda);


    /** overload this method to return the value of the objective function */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
            Number& obj_value);

    /** overload this method to return the vector of the gradient of
     *  the objective w.r.t. x */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
            Number* grad_f);

    /** overload this method to return the vector of constraint values */
    virtual bool eval_g(Index n, const Number* x, bool new_x,
            Index m, Number* g);

    /** overload this method to return the jacobian of the
     *  constraints. The vectors iRow and jCol only need to be set
     *  once. The first call is used to set the structure only (iRow
     *  and jCol will be non-NULL, and values will be NULL) For
     *  subsequent calls, iRow and jCol will be NULL. */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
            Index m, Index nele_jac, Index* iRow,
            Index *jCol, Number* values);

    /** overload this method to return the hessian of the
     *  lagrangian. The vectors iRow and jCol only need to be set once
     *  (during the first call). The first call is used to set the
     *  structure only (iRow and jCol will be non-NULL, and values
     *  will be NULL) For subsequent calls, iRow and jCol will be
     *  NULL. This matrix is symmetric - specify the lower diagonal
     *  only.  A default implementation is provided, in case the user
     *  wants to use quasi-Newton approximations to estimate the second
     *  derivatives and doesn't not neet to implement this method. */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
            Number obj_factor, Index m, const Number* lambda,
            bool new_lambda, Index nele_hess,
            Index* iRow, Index* jCol, Number* values);

    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,
            Index n, const Number* x, const Number* z_L, const Number* z_U,
            Index m, const Number* g, const Number* lambda,
            Number obj_value,
            const IpoptData* ip_data,
            IpoptCalculatedQuantities* ip_cq);

    //@}

private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    //TNLP();

    /** Copy Constructor */
    HEAT_NLP(const HEAT_NLP&);

    /** Overloaded Equals Operator */
    void operator=(const TNLP&);
    //@}


    MATRIXOP & data;
};





#endif /* HS071_NLP_HPP_ */
