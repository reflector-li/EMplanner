/*
Created by LinkX on 2022/3/31.
This code is using ipopt library to solve QP problem in frenet smoothing.
The detail usage of ipopt can be found at https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CPP
*/

#ifndef SMOOTHFRENET_SMOOTHSOLVER_H
#define SMOOTHFRENET_SMOOTHSOLVER_H

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "frenetUtils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>


using namespace Ipopt;

typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> T;

class Trajectory_NLP: public TNLP
{

private:
    bool printiterate_;
    double boundBuff_;
    Eigen::VectorXd Xr;
    int pointNumber = 0;
    double w1 = 2, w2= 2, w3 = 3; // w1->smooth weight, w2->sameSize weight, w3->length weight
    spMat Hessian;


public:
    std::vector<waypoint> Xoptimal;
    Trajectory_NLP( std::vector<waypoint> &Trajectory, bool printiterate,double boundBuff);

    ~Trajectory_NLP() override;


    virtual bool get_nlp_info(
            Index&          n,
            Index&          m,
            Index&          nnz_jac_g,
            Index&          nnz_h_lag,
            IndexStyleEnum& index_style
    );


    virtual bool get_bounds_info(
            Index   n,
            Number* x_l,
            Number* x_u,
            Index   m,
            Number* g_l,
            Number* g_u
    );

    virtual bool get_starting_point(
            Index   n,
            bool    init_x,
            Number* x,
            bool    init_z,
            Number* z_L,
            Number* z_U,
            Index   m,
            bool    init_lambda,
            Number* lambda
    );


    virtual bool eval_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number&       obj_value
    );


    virtual bool eval_grad_f(
            Index         n,
            const Number* x,
            bool          new_x,
            Number*       grad_f
    );


    virtual bool eval_g(
            Index         n,
            const Number* x,
            bool          new_x,
            Index         m,
            Number*       g
    );


    virtual bool eval_jac_g(
            Index         n,
            const Number* x,
            bool          new_x,
            Index         m,
            Index         nele_jac,
            Index*        iRow,
            Index*        jCol,
            Number*       values
    );


    virtual bool eval_h(
            Index         n,
            const Number* x,
            bool          new_x,
            Number        obj_factor,
            Index         m,
            const Number* lambda,
            bool          new_lambda,
            Index         nele_hess,
            Index*        iRow,
            Index*        jCol,
            Number*       values
    );


    virtual void finalize_solution(
            SolverReturn               status,
            Index                      n,
            const Number*              x,
            const Number*              z_L,
            const Number*              z_U,
            Index                      m,
            const Number*              g,
            const Number*              lambda,
            Number                     obj_value,
            const IpoptData*           ip_data,
            IpoptCalculatedQuantities* ip_cq
    );
    //@}

    bool intermediate_callback(
            AlgorithmMode              mode,
            Index                      iter,
            Number                     obj_value,
            Number                     inf_pr,
            Number                     inf_du,
            Number                     mu,
            Number                     d_norm,
            Number                     regularization_size,
            Number                     alpha_du,
            Number                     alpha_pr,
            Index                      ls_trials,
            const IpoptData*           ip_data,
            IpoptCalculatedQuantities* ip_cq
    );

};

void trajectorySolver(const SmartPtr<Trajectory_NLP>& mynlp,std::vector<waypoint> &optimal_traject);

#endif //SMOOTHFRENET_SMOOTHSOLVER_H
