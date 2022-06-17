//
// Created by 12971 on 2022/4/14.
//

#ifndef SMOOTHFRENET_OSQPSOLVER_H
#define SMOOTHFRENET_OSQPSOLVER_H

#include "OsqpEigen/OsqpEigen.h"
#include "frenetUtils.h"

typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> T;

typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> T;

class refLineSmoother{
private:
    spMat hessian;
    Eigen::VectorXd gradient;
    spMat linearMatrix;
    Eigen::VectorXd lowerBound,upperBound;
    int NumberOfVariables,NumberOfConstraints;
    OsqpEigen::Solver solver;

public:
    double w_cost_smooth, w_cost_ref, w_cost_length;
    double low_bound,upper_bound;
    Eigen::VectorXd Xinit;
    int size;

    refLineSmoother(double w_cost_smooth_,double w_cost_ref_,double w_cost_length_,double low_bound_,double upper_bound_,int size_);
    void updateRefLine(std::vector<waypoint> &trajectory);
    void getNewRefLine(std::vector<waypoint> &optimal_trajectory);

};

class qpPathSolver{
private:
    spMat hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double,Eigen::RowMajor> linearMatrix;
    Eigen::VectorXd lowerBound,upperBound;
    int NumberOfVariables,NumberOfConstraints;
    OsqpEigen::Solver solver;

public:
    int size;
    qpPlanConfigure config;
    frenet plan_start;
    double ds = 1;

    qpPathSolver(const qpPlanConfigure &config_,const frenet& start_point);
    void updateBound(const Eigen::VectorXd &low_conv_space, const Eigen::VectorXd &upper_conv_space,const frenet &start_point);
    void updateBoundWithDpPath(const Eigen::VectorXd &low_conv_space, const Eigen::VectorXd &upper_conv_space,const std::vector<frenet> &dp_path,const frenet &start_point);
    void getNewPath(std::vector<frenet> &optimal_path,const frenet &start_point);
};

#endif //SMOOTHFRENET_OSQPSOLVER_H
