//
// Created by 12971 on 2022/4/14.
//
#include "osqpSolver.h"

refLineSmoother::refLineSmoother(double w_cost_smooth_, double w_cost_ref_, double w_cost_length_, double low_bound_,
                                 double upper_bound_,int size_)
        :w_cost_smooth(w_cost_smooth_),w_cost_ref(w_cost_ref_),w_cost_length(w_cost_length_),
         low_bound(low_bound_),upper_bound(upper_bound_),size(size_)
{
    Xinit = Eigen::VectorXd::Zero(2*size);
    NumberOfVariables = 2*size;
    NumberOfConstraints = 2*size;

    /* create gradient vector*/
    gradient = -2*w_cost_ref*Xinit;

    /* create hessian matrix */
    // smooth cost matrix
    std::vector<T> triVector;
    spMat H1(2*size,2*(size-2));
    for(int i = 0;i<size-2;++i)
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
    triVector.clear();

    // close to original line cost matrix
    spMat H2(2*size,2*size);
    for(int i = 0;i<2*size;++i)
    {
        triVector.emplace_back(i,i,1);
    }
    H2.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    //length cost matrix
    spMat H3(2*size,2*(size-1));
    for(int i = 0;i<size-1;i++){
        std::vector<T> tempVec = {
                {2*i,2*i,1},
                {2*i+1,2*i+1,1},
                {2*i+2,2*i,-1},
                {2*i+3,2*i+1,-1},
        };
        triVector.insert(triVector.end(),tempVec.begin(),tempVec.end());
    }
    H3.setFromTriplets(triVector.begin(),triVector.end());
    hessian = 2*(w_cost_smooth*H1*H1.transpose()+w_cost_ref*H2+w_cost_length*H3*H3.transpose());


    /* create linearMatrix */
    triVector.clear();
    linearMatrix.resize(2*size,2*size);
    for(int i=0;i<2*size;++i)
    {
        triVector.emplace_back(i,i,1);
    }
    linearMatrix.setFromTriplets(triVector.begin(),triVector.end());

    /* create bound*/
    lowerBound = Xinit + low_bound*Eigen::VectorXd::Ones(2*size);
    upperBound = Xinit + upper_bound*Eigen::VectorXd::Ones(2*size);


    /* init solver*/
    solver.settings()->setWarmStart(true);


    solver.data()->setNumberOfVariables(NumberOfVariables);
    solver.data()->setNumberOfConstraints(NumberOfConstraints);
    if(!solver.data()->setHessianMatrix(hessian)) return;
    if(!solver.data()->setGradient(gradient)) return;
    if(!solver.data()->setLinearConstraintsMatrix(linearMatrix)) return;
    if(!solver.data()->setLowerBound(lowerBound)) return;
    if(!solver.data()->setUpperBound(upperBound)) return;

    /* initialize solver */
    if(!solver.initSolver()) return;
}

void refLineSmoother::updateRefLine(std::vector<waypoint> &trajectory)
{
    if(NumberOfVariables != 2*trajectory.size())
    {
        std::cout<<"path size not match the init size!"<<std::endl;
        return;
    }
    for(int i = 0;i<size;++i) {
        Xinit[2 * i] = trajectory.at(i).x;
        Xinit[2 * i + 1] = trajectory.at(i).y;
    }
    gradient = -2*w_cost_ref*Xinit;
    upperBound = Xinit + upper_bound*Eigen::VectorXd::Ones(2*size);
    lowerBound = Xinit + low_bound*Eigen::VectorXd::Ones(2*size);
    solver.updateBounds(lowerBound,upperBound);
    solver.updateGradient(gradient);
}

void refLineSmoother::getNewRefLine(std::vector<waypoint> &optimal_trajectory)
{
    if(solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) return ;
    Eigen::VectorXd Xr;
    Xr = solver.getSolution();
    optimal_trajectory.reserve(size);
    waypoint tempPoint;
    for( int i = 0; i < NumberOfVariables; i++ )
    {
        if(i%2 == 0)
            tempPoint.x = Xr[i];
        else {
            tempPoint.y = Xr[i];
            optimal_trajectory.push_back(tempPoint);
        }
    }

}



qpPathSolver::qpPathSolver(const qpPlanConfigure &config_,const frenet& start_point)
        :size(config_.size), config(config_), plan_start(start_point)
{
    NumberOfVariables = 3*size;
    NumberOfConstraints = 11*size-4;

    /*create hessian matrix*/
    // dl matrix
    spMat H_dl(NumberOfVariables,NumberOfVariables);
    std::vector<T> triVector;
    for(int i=0;i<size;++i)
        triVector.emplace_back(3*i+1,3*i+1,1);
    H_dl.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    // ddl matrix
    spMat H_ddl(NumberOfVariables,NumberOfVariables);
    for(int i =0;i<size;++i)
        triVector.emplace_back(3*i+2,3*i+2,1);
    H_ddl.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    // dddl matrix = P*P.transpose()
    spMat P_dddl(NumberOfVariables,size-1);
    for(int i=0;i<size-1;++i)
    {
        triVector.emplace_back(3*i+2,i,-1);
        triVector.emplace_back(3*i+5,i,1);
    }
    P_dddl.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    // refline matrix
    spMat H_ref(NumberOfVariables,NumberOfVariables);
    for(int i = 0;i<size;++i)
        triVector.emplace_back(3*i,3*i,1);
    H_ref.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    // center line matrix
    spMat H_center = H_ref;

    /* create end point constraints */
    spMat H_end_l(3*size,3*size);
    triVector.emplace_back(3*size -3,3*size-3,1);
    H_end_l.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    spMat H_end_dl(3*size,3*size);
    triVector.emplace_back(3*size-2,3*size-2,1);
    H_end_dl.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    spMat H_end_ddl(3*size,3*size);
    triVector.emplace_back(3*size-1,3*size-1,1);
    H_end_ddl.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    hessian = 2*(config.w_cost_dl*H_dl + config.w_cost_ddl*H_ddl + config.w_cost_dddl*P_dddl*P_dddl.transpose()
                 + config.w_cost_ref*H_ref + config.w_cost_center*H_center
                 + config.w_cost_end_l*H_end_l + config.w_cost_end_dl*H_end_dl + config.w_cost_end_ddl*H_end_ddl);


    Eigen::VectorXd low_conv_space= -6*Eigen::VectorXd::Ones(size);
    Eigen::VectorXd upper_conv_space = 6*Eigen::VectorXd::Ones(size);

    /* create gradient vector */
    double end_l_desire = 0;
    double end_dl_desire = 0;
    double end_ddl_desire = 0;
    gradient = Eigen::VectorXd::Zero(3*size);



    /*create equation constraints matrix Aeq */
    spMat Aeq(2*size-2,3*size);
    for(int i = 0;i<size-1;++i)
    {
        std::vector<T> tem_vec = {
                {2*i,3*i,1},{2*i,3*i+1,ds},{2*i,3*i+2,ds*ds/3},{2*i,3*i+3,-1},{2*i,3*i+5,ds*ds/6},
                {2*i+1,3*i+1,1},{2*i+1,3*i+2,ds/2},{2*i+1,3*i+4,-1},{2*i+1,3*i+5,ds/2}
        };
        triVector.insert(triVector.end(),tem_vec.begin(),tem_vec.end());
    }
    Aeq.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();
    Eigen::VectorXd low_eq = Eigen::VectorXd::Zero(2*size-2);
    Eigen::VectorXd upper_eq = Eigen::VectorXd::Zero(2*size-2);


    /* create inequation constraints matrix A_ieq */
    double d1 = config.host_d1;
    double d2 = config.host_d2;
    double w = config.host_w;
    spMat A_ieq(4*size,3*size);
    // neglect the start point
    for(int i =1;i<size;++i)
    {
        std::vector<T> tem_vec = {
                {4*i,3*i,1},{4*i,3*i+1,d1},
                {4*i+1,3*i,1},{4*i+1,3*i+1,d1},
                {4*i+2,3*i,1},{4*i+2,3*i+1,-d2},
                {4*i+3,3*i,1},{4*i+3,3*i+1,-d2}
        };
        triVector.insert(triVector.end(),tem_vec.begin(),tem_vec.end());
    }
    A_ieq.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();


    /* create inequation constraints matrix low_ieq and upper_ieq */
    int vehicle_front_index = std::ceil(d1/ds), vehicle_back_index = std::ceil(d2/ds);
    Eigen::VectorXd low_ieq = Eigen::VectorXd::Zero(4*size);
    Eigen::VectorXd upper_ieq = Eigen::VectorXd::Zero(4*size);




    /* create start point constraints */
    spMat A_start(3*size,3*size);
    for(int i =0;i<3*size;++i)
    {
        triVector.emplace_back(i,i,1);
    }
    A_start.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();


    Eigen::VectorXd low_start = -1e8*Eigen::VectorXd::Ones(3*size);
    Eigen::VectorXd upper_start = 1e8*Eigen::VectorXd::Ones(3*size);
    low_start[0] = plan_start.l;
    low_start[1] = plan_start.l_diff;
    low_start[2] = plan_start.l_2diff;
    upper_start[0] = plan_start.l;
    upper_start[1] = plan_start.l_diff;
    upper_start[2] = plan_start.l_2diff;



    /* create dl and ddl difference minus constraints */
    spMat A_dl_minus(size-1,3*size);
    for(int i =0;i<size-1;++i)
    {
        std::vector<T> temp_vec = {{i,3*i+1,-1},{i,3*i+4,1}};
        triVector.insert(triVector.end(),temp_vec.begin(),temp_vec.end());
    }
    A_dl_minus.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    Eigen::VectorXd low_dl_minus = -config.delta_dl_max*Eigen::VectorXd::Ones(size-1);
    Eigen::VectorXd upper_dl_minus = config.delta_dl_max*Eigen::VectorXd::Ones(size-1);

    spMat A_ddl_minus(size-1,3*size);
    for(int i =0;i<size-1;++i)
    {
        std::vector<T> temp_vec = {{i,3*i+2,-1},{i,3*i+5,1}};
        triVector.insert(triVector.end(),temp_vec.begin(),temp_vec.end());
    }
    A_ddl_minus.setFromTriplets(triVector.begin(),triVector.end());
    triVector.clear();

    Eigen::VectorXd low_ddl_minus = -config.delta_ddl_max*Eigen::VectorXd::Ones(size-1);
    Eigen::VectorXd upper_ddl_minus = config.delta_ddl_max*Eigen::VectorXd::Ones(size-1);



    /* get final matrixes */
    linearMatrix.resize(11*size-4,3*size);
    linearMatrix.topRows(2*size-2) = Aeq;
    linearMatrix.middleRows(2*size-2,4*size) = A_ieq;
    linearMatrix.middleRows(6*size-2,3*size) = A_start;
    linearMatrix.middleRows(9*size-2,size-1) = A_dl_minus;
    linearMatrix.bottomRows(size-1) = A_ddl_minus;

    lowerBound.resize(11*size-4);
    upperBound.resize(11*size-4);

    //set for lowBound
    lowerBound.head(2*size-2) = low_eq;
    lowerBound.segment(2*size-2,4*size) = low_ieq;
    lowerBound.segment(6*size-2,3*size) = low_start;
    lowerBound.segment(9*size-2,size-1) = low_dl_minus;
    lowerBound.tail(size-1) = low_ddl_minus;
    // std::cout<<"lowerBound_ori: \n"<<lowerBound<<std::endl;

    //set for upperBound
    upperBound.head(2*size-2) = upper_eq;
    upperBound.segment(2*size-2,4*size) = upper_ieq;
    upperBound.segment(6*size-2,3*size) = upper_start;
    upperBound.segment(9*size-2,size-1) = upper_dl_minus;
    upperBound.tail(size-1) = upper_ddl_minus;


    /* init solver*/
    solver.settings()->setWarmStart(true);

    solver.data()->setNumberOfVariables(NumberOfVariables);
    solver.data()->setNumberOfConstraints(NumberOfConstraints);
    if(!solver.data()->setHessianMatrix(hessian)) return;
    if(!solver.data()->setGradient(gradient)) return;
    if(!solver.data()->setLinearConstraintsMatrix(linearMatrix)) return;
    if(!solver.data()->setLowerBound(lowerBound)) return;
    if(!solver.data()->setUpperBound(upperBound)) return;

    /* initialize solver */
    if(!solver.initSolver()) return;
}

void qpPathSolver::updateBound(const Eigen::VectorXd &low_conv_space, const Eigen::VectorXd &upper_conv_space,const frenet &start_point)
{
    if(low_conv_space.size() != size || upper_conv_space.size()!=size)
    {
        std::cout<<"bound size not match init size!"<<std::endl;
        return;
    }
    // low_eq and upper_eq
    Eigen::VectorXd low_eq = Eigen::VectorXd::Zero(2*size-2);
    Eigen::VectorXd upper_eq = Eigen::VectorXd::Zero(2*size-2);

    // low_ieq and upper_ieq
    double d1 = config.host_d1;
    double d2 = config.host_d2;
    double w = config.host_w;
    int vehicle_front_index = std::ceil(d1/ds), vehicle_back_index = std::ceil(d2/ds);
    Eigen::VectorXd low_ieq = Eigen::VectorXd::Zero(4*size);
    Eigen::VectorXd upper_ieq = Eigen::VectorXd::Zero(4*size);
    // neglect the start point
    for(int i =1;i<size;i++)
    {
        int front_index = std::min(i+vehicle_front_index,size-1);
        int back_index = std::max(i-vehicle_back_index,0);
        low_ieq[4*i] = low_conv_space[front_index] - w/2;
        low_ieq[4*i+1] = low_conv_space[front_index] + w/2;
        low_ieq[4*i+2] = low_conv_space[back_index] - w/2;
        low_ieq[4*i+3] = low_conv_space[back_index] + w/2;

        upper_ieq[4*i] = upper_conv_space[front_index] - w/2;
        upper_ieq[4*i+1] = upper_conv_space[front_index] + w/2;
        upper_ieq[4*i+2] = upper_conv_space[back_index] - w/2;
        upper_ieq[4*i+3] = upper_conv_space[back_index] + w/2;
    }

    // low_start and upper_start
    plan_start = start_point;
    Eigen::VectorXd low_start = -1e8*Eigen::VectorXd::Ones(3*size);
    Eigen::VectorXd upper_start = 1e8*Eigen::VectorXd::Ones(3*size);
    low_start[0] = plan_start.l;
    low_start[1] = plan_start.l_diff;
    low_start[2] = plan_start.l_2diff;
    upper_start[0] = plan_start.l;
    upper_start[1] = plan_start.l_diff;
    upper_start[2] = plan_start.l_2diff;
    // add dl and ddl constraints
    for(int i=1;i<size;++i)
    {
        low_start[3*i+1] = -config.dl_max;
        low_start[3*i+2] = -config.ddl_max;
        upper_start[3*i+1] = config.dl_max;
        upper_start[3*i+2] = config.ddl_max;
    }

    /* update gradient vector */
    double end_l_desire = 0;
    double end_dl_desire = 0;
    double end_ddl_desire = 0;
    for(int i =0;i<size;++i)
    {
        // double x = (upper_conv_space[i]+low_conv_space[i])/2;
        double x = low_conv_space[i] + 1.5;
        if(abs(x)>0.3)
            gradient[3*i] = -2*config.w_cost_center*x;
        else
            gradient[3*i] = -2*x;

    }
    gradient[3*size-3] = gradient[3*size-3] - 2*config.w_cost_end_l*end_l_desire;
    gradient[3*size-2] = gradient[3*size-2] - 2*config.w_cost_end_dl*end_dl_desire;
    gradient[3*size-1] = gradient[3*size-1] - 2*config.w_cost_end_ddl*end_ddl_desire;


    lowerBound.segment(2*size-2,4*size) = low_ieq;
    lowerBound.segment(6*size-2,3*size) = low_start;
    // std::cout<<"------------"<<std::endl;
    // std::cout<<"upper_bound: \n"<<upper_conv_space.head(6)<<std::endl;
    // std::cout<<"low_bound: \n"<<low_conv_space.head(6)<<std::endl;
    // std::cout<<"gradient:\n"<<gradient.head(9)<<std::endl;
    // std::cout<<"low_ieq: \n"<<low_ieq.head(16)<<std::endl;
    // std::cout<<"high_ieq: \n"<<upper_ieq.head(16)<<std::endl;
    //std::cout<<"low_start: \n"<<low_start<<std::endl;
    // std::cout<<"high_start: \n"<<upper_start.head(12)<<std::endl;



    //set for upperBound
    upperBound.segment(2*size-2,4*size) = upper_ieq;
    upperBound.segment(6*size-2,3*size) = upper_start;
    
    /* update bound and gradient */
    solver.updateBounds(lowerBound,upperBound);
    solver.updateGradient(gradient);

}

void qpPathSolver::getNewPath(std::vector<frenet> &optimal_path,const frenet &start_point)
{
    if(solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError) return ;
    Eigen::VectorXd Xr;
    Xr = solver.getSolution();
    optimal_path.reserve(size);
    frenet tempPoint;
    int count = 0;
    for( int i = 0; i < NumberOfVariables; i++ )
    {
        if(i%3 == 0)
            tempPoint.l = Xr[i];
        else if(i%3 == 1)
        {
            tempPoint.l_diff = Xr[i];
        }
        else
        {
            tempPoint.l_2diff = Xr[i];
            tempPoint.s = start_point.s+count*ds;
            optimal_path.push_back(tempPoint);
            count++;
        }
    }
}

void qpPathSolver::updateBoundWithDpPath(const Eigen::VectorXd &low_conv_space, const Eigen::VectorXd &upper_conv_space,const std::vector<frenet> &dp_path,const frenet &start_point)
{
    if(low_conv_space.size() != size || upper_conv_space.size()!=size || dp_path.size() != size)
    {
        std::cout<<"bound size not match init size!"<<std::endl;
        return;
    }
    // low_eq and upper_eq
    Eigen::VectorXd low_eq = Eigen::VectorXd::Zero(2*size-2);
    Eigen::VectorXd upper_eq = Eigen::VectorXd::Zero(2*size-2);

    // low_ieq and upper_ieq
    double d1 = config.host_d1;
    double d2 = config.host_d2;
    double w = config.host_w;
    int vehicle_front_index = std::ceil(d1/ds), vehicle_back_index = std::ceil(d2/ds);
    Eigen::VectorXd low_ieq(4*size),upper_ieq(4*size);
    for(int i =0;i<size;i++)
    {
        int front_index = std::min(i+vehicle_front_index,size-1);
        int back_index = std::max(i-vehicle_back_index,0);
        low_ieq[4*i] = low_conv_space[front_index] - w/2;
        low_ieq[4*i+1] = low_conv_space[front_index] + w/2;
        low_ieq[4*i+2] = low_conv_space[back_index] - w/2;
        low_ieq[4*i+3] = low_conv_space[back_index] + w/2;

        upper_ieq[4*i] = upper_conv_space[front_index] - w/2;
        upper_ieq[4*i+1] = upper_conv_space[front_index] + w/2;
        upper_ieq[4*i+2] = upper_conv_space[back_index] - w/2;
        upper_ieq[4*i+3] = upper_conv_space[back_index] + w/2;
    }

    // low_start and upper_start
    plan_start = start_point;
    Eigen::VectorXd low_start = -1e8*Eigen::VectorXd::Ones(3*size);
    Eigen::VectorXd upper_start = 1e8*Eigen::VectorXd::Ones(3*size);
    low_start[0] = plan_start.l;
    low_start[1] = plan_start.l_diff;
    low_start[2] = plan_start.l_2diff;
    upper_start[0] = plan_start.l;
    upper_start[1] = plan_start.l_diff;
    upper_start[2] = plan_start.l_2diff;

    /* update gradient vector */
    double end_l_desire = 0;
    double end_dl_desire = 0;
    double end_ddl_desire = 0;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(3*size);
    for(int i =0;i<size;++i)
//        f[3*i] = low_conv_space[i]+2;
        f[3*i] = dp_path.at(i).l;
    gradient = -2*config.w_cost_center*f;
    gradient[3*size-3] = gradient[3*size-3] - 2*config.w_cost_end_l*end_l_desire;
    gradient[3*size-2] = gradient[3*size-2] - 2*config.w_cost_end_dl*end_dl_desire;
    gradient[3*size-1] = gradient[3*size-1] - 2*config.w_cost_end_ddl*end_ddl_desire;


    lowerBound.head(2*size-2) = low_eq;
    lowerBound.segment(2*size-2,4*size) = low_ieq;
    lowerBound.tail(3*size) = low_start;

    //set for upperBound
    upperBound.head(2*size-2) = upper_eq;
    upperBound.segment(2*size-2,4*size) = upper_ieq;
    upperBound.tail(3*size) = upper_start;


    /* update bound and gradient */
    solver.updateBounds(lowerBound,upperBound);
    solver.updateGradient(gradient);
}