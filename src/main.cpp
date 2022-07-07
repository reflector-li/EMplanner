#include <iostream>
#include "frenetUtils.h"
#include "smoothSolver.h"
#include "osqpSolver.h"
#include "matplotlibcpp.h"

#include <vector>

namespace plt = matplotlibcpp;
using namespace std;
using namespace Ipopt;




int main() {
    
    // dp algorithm parameter set
    initPlanConfigure dp_config;
    dp_config.w_cost_obs = 1e6;
    dp_config.w_cost_dl = 300;
    dp_config.w_cost_ddl = 2000;
    dp_config.w_cost_dddl = 10000;
    dp_config.w_cost_ref = 20;
    dp_config.rows = 9;
    dp_config.cols = 6;
    dp_config.sample_s = 10;
    dp_config.sample_l = 1;
    
    
    // qp algorithm parameter set
    qpPlanConfigure qp_config;
    qp_config.size = 60; //60
    qp_config.w_cost_dl = 100000; //25000 for original
    qp_config.w_cost_ddl = 50;
    qp_config.w_cost_dddl = 20;
    qp_config.w_cost_ref = 2;
    qp_config.w_cost_end_l = 15;
    qp_config.w_cost_end_dl = 15;
    qp_config.w_cost_end_ddl = 15;
    qp_config.w_cost_center = 1200; //1200
    qp_config.host_d1 = 3;
    qp_config.host_d2 = 3;
    qp_config.host_w = 1.63;
    qp_config.obs_width = 2;
    qp_config.obs_length = 5;
    qp_config.delta_dl_max = 10;
    qp_config.delta_ddl_max = 20;
    qp_config.dl_max = 10;
    qp_config.ddl_max = 3;

    // init reference line smooth, reference line must have 181 points
    refLineSmoother smoother(2,2,3,-0.2,0.2,181);
    frenet init_point;

    //qp solver init
    qpPathSolver final_path_solver(qp_config,init_point);

    //global obstacle set
    vector<planPointInfo> obs_array;
    readObs(obs_array,"../data/obstacles.txt");

    //vehicle init position
    planPointInfo vehicle;
    vehicle.x = 55;
    vehicle.y = 448; //460
    vehicle.time = 0;

    
    //read global path and resample it with 1m 
    vector<waypoint> global_path,left_edge,right_edge;
    readTraje(global_path,"../data/waypoints.txt");
    resampleTraje(1,global_path);
    for(int i =0;i<global_path.size();++i)
    {
        getDirAndK(global_path, i);
    }
    createEdge(global_path,left_edge,right_edge);


    int control = 0;

    vector<waypoint> pre_path,stitch_path;
    int ref_preMatch_index = -1; // used to get reference line
    while (control<global_path.size())
    {

        int next_index = 0;

        // get reference line and smoother it 
        vector<waypoint> ori_ref_path,ref_path;
        int ref_match_index = getMatchPoint(vehicle,ref_preMatch_index,global_path);
        ref_preMatch_index = ref_match_index;
        getRefTrajectory(ref_match_index,global_path,ori_ref_path);
        smoother.updateRefLine(ori_ref_path);
        smoother.getNewRefLine(ref_path);
        for(int i =0;i<ref_path.size();++i)
        {
            getDirAndK(ref_path, i);
        }
    

        // get start plan point
        planPointInfo start_plan_point;
        getStartPlanPoint(pre_path,vehicle,start_plan_point,stitch_path);


        // select obstacles
        vector<planPointInfo> new_obs_array;
        obsSelect(vehicle,obs_array,new_obs_array);


        // covert start point and obstacles to frenet coordinate
        int fre_match_index = getMatchPoint(vehicle,-1,ref_path);
        waypoint veh_match_point = ref_path.at(fre_match_index);
        waypoint veh_project_point;
        getProjectPoint(vehicle,veh_match_point,veh_project_point);
        Eigen::VectorXd index_S;
        index2S(fre_match_index,veh_project_point,ref_path,index_S);
        //    cout<<index_S<<endl;


        int spp_match_index = getMatchPoint(start_plan_point,-1,ref_path);
        waypoint spp_match_point = ref_path.at(spp_match_index);
        waypoint spp_project_point;
        getProjectPoint(start_plan_point,spp_match_point,spp_project_point);
        frenet spp_frenet;
        cardesian2Frenet(start_plan_point,spp_match_index,spp_project_point,ref_path,index_S,spp_frenet);

        vector<frenet> obs_array_frenet;
        for(auto &obs:new_obs_array)
        {
            int obs_match_index = getMatchPoint(obs,-1,ref_path);
            waypoint obs_match_point = ref_path.at(obs_match_index);
            waypoint obs_project_point;
            getProjectPoint(obs,obs_match_point,obs_project_point);
            frenet obs_frenet;
            cardesian2Frenet(obs,obs_match_index,obs_project_point,ref_path,index_S,obs_frenet);
            obs_array_frenet.push_back(obs_frenet);
        }

        // get init plan value
        vector<frenet> init_path;
        getInitPlanValue(obs_array_frenet,spp_frenet,dp_config,init_path);
        trajectoryInterp(init_path);
        vector<waypoint> no_stitch_path;
        frenet2Cardesian(init_path, ref_path,index_S,no_stitch_path);

        // qp plan
        Eigen::VectorXd low_bound,upper_bound;
        getCovexBound(init_path,obs_array_frenet,qp_config,low_bound,upper_bound);


        final_path_solver.updateBound(low_bound,upper_bound,spp_frenet);
        // final_path_solver.updateBoundWithDpPath(low_bound,upper_bound,init_path,spp_frenet);
        vector<frenet> qp_frenet_path;
        final_path_solver.getNewPath(qp_frenet_path,spp_frenet);
        vector<waypoint> qp_final_path;
        frenet2Cardesian(qp_frenet_path,ref_path,index_S,qp_final_path);
        timeAndVelocity(qp_final_path,vehicle.time+0.1);


        // stitch the pre_plan path
        vector<waypoint> final_path;
        copy(stitch_path.begin(),stitch_path.end(), back_inserter(final_path));
        final_path.insert(final_path.end(),qp_final_path.begin(),qp_final_path.end());

        plt::cla();
        plt::plotTrajectory(global_path);
        plt::plotTrajectory(left_edge,"green");
        plt::plotTrajectory(right_edge,"green");
        plt::plotTrajectory(ref_path,"black");
        plt::plotTrajectory(qp_final_path,"purple");
        plt::plotTrajectory(no_stitch_path,"y");
//        plt::plotTrajectory(final_path,"r");
        plt::plotTrajectory(new_obs_array,".r");
        plt::plot(vector<double> {vehicle.x},vector<double> {vehicle.y},"vc");
//        plt::plot(vector<double>{no_stitch_path.at(0).x},vector<double>{no_stitch_path.at(0).y},"vg");

        plt::xlim(vehicle.x - 40,vehicle.x + 80);
        plt::ylim(vehicle.y - 30,vehicle.y + 80);
        plt::title("plot frenet_path");
        plt::grid(true);
        plt::pause(0.5);

        control = ref_match_index+150;
        vehicle = covertFromWaypoint(qp_final_path.at(0));
        pre_path = final_path;
    }


//    int ref_match_index = -1,ref_preMatch_index = -1; // used to get reference line
//    int next_index = 0;
//
//    vector<waypoint> ori_ref_path,ref_path;
//    ref_match_index = getMatchPoint(vehicle,ref_preMatch_index,global_path);
//    getRefTrajectory(ref_match_index,global_path,ori_ref_path);
//    smoother.updateRefLine(ori_ref_path);
//    smoother.getNewRefLine(ref_path);
//    for(int i =0;i<ref_path.size();++i)
//    {
//        getDirAndK(ref_path, i);
//    }
//
//    // get start plan point
//    planPointInfo start_plan_point;
//    getStartPlanPoint(pre_path,vehicle,start_plan_point,stitch_path);
//
//
//    // select obstacles
//    vector<planPointInfo> new_obs_array;
//    obsSelect(vehicle,obs_array,new_obs_array);
//
//
//    // covert start point and obstacles to frenet coordinate
//    int fre_match_index = getMatchPoint(vehicle,-1,ref_path);
//    waypoint veh_match_point = ref_path.at(fre_match_index);
//    waypoint veh_project_point;
//    getProjectPoint(vehicle,veh_match_point,veh_project_point);
//    Eigen::VectorXd index_S;
//    index2S(fre_match_index,veh_project_point,ref_path,index_S);
////    cout<<index_S<<endl;
//
//
//    int spp_match_index = getMatchPoint(start_plan_point,-1,ref_path);
//    waypoint spp_match_point = ref_path.at(spp_match_index);
//    waypoint spp_project_point;
//    getProjectPoint(start_plan_point,spp_match_point,spp_project_point);
//    frenet spp_frenet;
//    cardesian2Frenet(start_plan_point,spp_match_index,spp_project_point,ref_path,index_S,spp_frenet);
//
//    vector<frenet> obs_array_frenet;
//    for(auto &obs:new_obs_array)
//    {
//        int obs_match_index = getMatchPoint(obs,-1,ref_path);
//        waypoint obs_match_point = ref_path.at(obs_match_index);
//        waypoint obs_project_point;
//        getProjectPoint(obs,obs_match_point,obs_project_point);
//        frenet obs_frenet;
//        cardesian2Frenet(obs,obs_match_index,obs_project_point,ref_path,index_S,obs_frenet);
//        //cout<<obs_frenet.s<<" "<<obs_frenet.l<<endl;
//        obs_array_frenet.push_back(obs_frenet);
//    }
//
//    // get init plan value
//    vector<frenet> init_path;
//    getInitPlanValue(obs_array_frenet,spp_frenet,dp_config,init_path);
//    trajectoryInterp(init_path);
//    vector<waypoint> no_stitch_path;
//    frenet2Cardesian(init_path, ref_path,index_S,no_stitch_path);
//    timeAndVelocity(no_stitch_path,vehicle.time+0.1);
//
//    // qp plan
//    Eigen::VectorXd low_bound,upper_bound;
//    getCovexBound(init_path,obs_array_frenet,qp_config,low_bound,upper_bound);
////    writeEigen(low_bound,"../data/lowBound.txt");
////    writeEigen(upper_bound,"../data/uppperBound.txt");
//
////    final_path_solver.updateBound(low_bound,upper_bound,spp_frenet);
//    final_path_solver.updateBoundWithDpPath(low_bound,upper_bound,init_path,spp_frenet);
//
//    vector<frenet> qp_frenet_path;
//    final_path_solver.getNewPath(qp_frenet_path,spp_frenet);
//    vector<waypoint> qp_final_path;
//    frenet2Cardesian(qp_frenet_path,ref_path,index_S,qp_final_path);
//    timeAndVelocity(qp_final_path,vehicle.time+0.1);
//
//
//    // stitch the pre_plan path
//    vector<waypoint> final_path;
//    copy(stitch_path.begin(),stitch_path.end(), back_inserter(final_path));
//    final_path.insert(final_path.end(),no_stitch_path.begin(),no_stitch_path.end());
//
//    plt::plotTrajectory(global_path);
//    plt::plotTrajectory(left_edge,"green");
//    plt::plotTrajectory(right_edge,"green");
//    plt::plotTrajectory(ref_path,"black");
//    plt::plotTrajectory(final_path,"red");
//    plt::plotTrajectory(qp_final_path,"purple");
//    plt::plotTrajectory(new_obs_array,".r");
//    plt::plot(vector<double> {vehicle.x},vector<double> {vehicle.y},"vc");
//    plt::xlim(vehicle.x - 20,vehicle.x + 50);
//    plt::ylim(vehicle.y - 30,vehicle.y + 50);
//    plt::title("plot frenet_path");
//    plt::grid(true);
//    plt::show();


    return 0;
}




