/*
Created by LinkX on 2022/3/31.
This code is using ipopt library to solve QP problem in frenet smoothing.
The detail usage of ipopt can be found at https://coin-or.github.io/Ipopt/INTERFACES.html#INTERFACE_CPP
*/

#include "smoothSolver.h"

Trajectory_NLP::Trajectory_NLP(std::vector<waypoint> &Trajectory, bool printiterate,double boundBuff):printiterate_(printiterate),boundBuff_(boundBuff)
{
    pointNumber = static_cast<int>(Trajectory.size());
    Xr.resize(2*pointNumber);
    Xoptimal.reserve(pointNumber);
    for(int i = 0;i<pointNumber;++i) {
        Xr[2 * i] = Trajectory.at(i).x;
        Xr[2 * i + 1] = Trajectory.at(i).y;
    }

    std::vector<T> triVector;

    // create H1
    spMat H1(2*pointNumber,2*(pointNumber-2));
    for(int i = 0;i<pointNumber-2;++i)
    {
        std::vector<T> tempVec = {
                {2*i,2*i,1},
                {2*i+1,2*i+1,1},
                {2*i+2,2*i,-2},
                {2*i+3,2*i+1,-2},
                {2*i+4,2*i,1},
                {2*i+5,2*i+1,1}
        };
        triVector.insert(triVector.end(),tempVec.begin(),tempVec.end());
    }
    H1.setFromTriplets(triVector.begin(),triVector.end());

    // create H2
    triVector.clear();
    for(int i = 0;i<pointNumber;++i)
    {
        triVector.emplace_back(i,i,1);
    }
    spMat H2(2*pointNumber,2*pointNumber);
    H2.setFromTriplets(triVector.begin(),triVector.end());

    // create H3
    triVector.clear();
    for(int i = 0;i<pointNumber-1;i++){
        std::vector<T> tempVec = {
                {2*i,2*i,1},
                {2*i+1,2*i+1,1},
                {2*i+2,2*i,-1},
                {2*i+3,2*i+1,-1},
        };
        triVector.insert(triVector.end(),tempVec.begin(),tempVec.end());
    }
    spMat H3(2*pointNumber,2*(pointNumber-1));
    H3.setFromTriplets(triVector.begin(),triVector.end());

    spMat H = w1*H1*H1.transpose()+w2*H2+w3*H3*H3.transpose();

    Hessian = H.triangularView<Eigen::Lower>();

}

Trajectory_NLP::~Trajectory_NLP()
= default;

bool Trajectory_NLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style)
{
   n = 2*pointNumber;
   m = 0;
   nnz_jac_g = 0;
   nnz_h_lag = static_cast<int>(Hessian.nonZeros());
   index_style = TNLP::C_STYLE;
   return true;
}

bool Trajectory_NLP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u)
{
    assert(n == 2*pointNumber);
    assert(m == 0);
    for(int i =0;i<n;++i)
    {
        x_l[i] = Xr[i] - boundBuff_;
        x_u[i] = Xr[i] + boundBuff_;
    }

    return true;

}

bool Trajectory_NLP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                                        bool init_lambda, Number *lambda)
{
    assert(n == 2*pointNumber);
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    for(int i = 0;i<n;++i)
    {
        x[i] = Xr[i];
    }
    return true;
}

bool Trajectory_NLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
{
    assert(n == 2*pointNumber);
    Eigen::VectorXd temp_x(n);
    for(int i = 0;i<n;++i)
    {
        temp_x[i] = x[i];
    }
    double firstDegree = 2*w2*Xr.transpose()*temp_x;
    obj_value = temp_x.transpose()*Hessian*temp_x-firstDegree;
    return true;
}

bool Trajectory_NLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
{
    assert(n == 2*pointNumber);
    Eigen::VectorXd temp_x(n);
    Eigen::VectorXd grad(n);
    for(int i = 0;i<n;++i)
    {
        temp_x[i] = x[i];
    }
    grad = 2*temp_x.transpose()*Hessian-2*w2*Xr.transpose();
    for(int i =0;i<n;++i)
    {
        grad_f[i] = grad[i];
    }
    return true;
}

bool Trajectory_NLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g)
{
    assert(n == 2*pointNumber);
    assert(m ==0);
    return true;
}

bool Trajectory_NLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                                Number *values)
{
    assert(n == 2*pointNumber);
    assert(m ==0);
    return true;
}

bool Trajectory_NLP::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m, const Number *lambda,
                            bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values)
{
    assert(n == 2*pointNumber);
    assert(m ==0);
    if( values == nullptr)
    {
        Index idx = 0;
        for(int i = 0;i<Hessian.outerSize();++i)
     {
         for(spMat::InnerIterator it(Hessian,i);it;++it){
             iRow[idx] = static_cast<int>(it.row());
             jCol[idx] = static_cast<int>(it.col());
             idx++;
         }
     }
        assert(idx == nele_hess);
    }
    else{
        Index idx = 0;
        for(int i = 0;i<Hessian.outerSize();++i)
        {
            for(spMat::InnerIterator it(Hessian,i);it;++it){
                values[idx++] = obj_factor*it.value();
            }
        }
    }

    return true;
}

void Trajectory_NLP::finalize_solution(
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
)
{
    std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
    waypoint tempPoint;
    for( Index i = 0; i < n; i++ )
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
        if(i%2 == 0)
            tempPoint.x = x[i];
        else {
            tempPoint.y = x[i];
            Xoptimal.push_back(tempPoint);
        }
    }

    std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
    for( Index i = 0; i < n; i++ )
    {
        std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
    }
    for( Index i = 0; i < n; i++ )
    {
        std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
    }

    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;
}

bool Trajectory_NLP::intermediate_callback(
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
)
{
    return true;
}

void trajectorySolver(const SmartPtr<Trajectory_NLP>& mynlp,std::vector<waypoint> &optimal_traject)
{
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetNumericValue("tol", 3.82e-6);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded )
    {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if( status == Solve_Succeeded )
    {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else
    {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    optimal_traject = mynlp->Xoptimal;
}
