#include <iostream>
#include "frenetUtils.h"
#include "smoothSolver.h"
#include "osqpSolver.h"
#include <vector>
#include "visualization.h"
#include <mcap/writer.hpp>
#include <mcap/reader.hpp>
#include "BuildFileDescriptorSet.h"

using namespace std;
using namespace Ipopt;

#define NS_PER_MS 1000000
#define NS_PER_S 1000000000

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <output.mcap>" << std::endl;
        return 1;
    }
    const char* outputFilename = argv[1];

    mcap::McapWriter writer;
    {
        auto options = mcap::McapWriterOptions("");
        const auto res = writer.open(outputFilename, options);
        if (!res.ok()) {
        std::cerr << "Failed to open " << outputFilename << " for writing: " << res.message
                    << std::endl;
        return 1;
        }
    }


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


    mcap::ChannelId scene_channel;
    {
        mcap::Schema schema(
            "foxglove.SceneUpdate", "protobuf",
            foxglove::BuildFileDescriptorSet(foxglove::SceneUpdate::descriptor()).SerializeAsString());
        writer.addSchema(schema);
        mcap::Channel channel("road_msg", "protobuf", schema.id);
        writer.addChannel(channel);
        scene_channel = channel.id;
    }

    mcap::ChannelId data_channel;
    {
        mcap::Schema schema(
            "foxglove.Point2", "protobuf",
            foxglove::BuildFileDescriptorSet(foxglove::Point2::descriptor()).SerializeAsString());
        writer.addSchema(schema);
        mcap::Channel channel("along_edge", "protobuf", schema.id);
        writer.addChannel(channel);
        data_channel = channel.id;
    }

    mcap::ChannelId tf_channel;
    {
        mcap::Schema schema(
            "foxglove.FrameTransform", "protobuf",
            foxglove::BuildFileDescriptorSet(foxglove::FrameTransform::descriptor()).SerializeAsString());
        writer.addSchema(schema);
        mcap::Channel channel("root2bese_link", "protobuf", schema.id);
        writer.addChannel(channel);
        tf_channel = channel.id;
    }
    int control = 0;
    vector<waypoint> pre_path,stitch_path;
    int ref_preMatch_index = -1; // used to get reference line

    mcap::Timestamp startTime = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                    std::chrono::system_clock::now().time_since_epoch())
                                    .count();
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
            obs.dirAngle = obs_project_point.dirAngle;
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


        // tf message
        mcap::Timestamp curr_time = startTime + (static_cast<uint64_t>(control) * 100 * NS_PER_MS);
        
        foxglove::FrameTransform tf;
        auto tf_timestamp = tf.mutable_timestamp();
        tf_timestamp->set_seconds(static_cast<int64_t>(curr_time) / NS_PER_S);
        tf_timestamp->set_nanos(static_cast<int>(curr_time % NS_PER_S));
        RootToEgoTF(tf,"root","base_link",vehicle);
        std::string tf_serializedMsg = tf.SerializeAsString();
        
        mcap::Message tf_msg;
        tf_msg.channelId = tf_channel;
        tf_msg.publishTime = curr_time;
        tf_msg.logTime = curr_time;
        tf_msg.data = reinterpret_cast<const std::byte*>(tf_serializedMsg.data());
        tf_msg.dataSize = tf_serializedMsg.size();
        writer.write(tf_msg);

        // scene message
        foxglove::SceneUpdate scene;
        auto *entity = scene.add_entities();
        auto scene_timestamp = entity->mutable_timestamp();
        scene_timestamp->set_seconds(static_cast<int64_t>(curr_time) / NS_PER_S);
        scene_timestamp->set_nanos(static_cast<int>(curr_time % NS_PER_S));
        entity->set_frame_id("root");
        visualizeLine(entity,global_path,{.thickness = 1.0,
                                          .type = 2,
                                          .color_r = 1.0,
                                          .color_g = 0.64,
                                          .color_b = 0.0,
                                          .color_a = 1.0});
        visualizeLine(entity,left_edge,{.thickness = 1.0,
                                          .type = 0,
                                          .color_r = 0.0,
                                          .color_g = 0.0,
                                          .color_b = 0.0,
                                          .color_a = 1.0});
        visualizeLine(entity,right_edge,{.thickness = 1.0,
                                          .type = 0,
                                          .color_r = 0.0,
                                          .color_g = 0.0,
                                          .color_b = 0.0,
                                          .color_a = 1.0});
        visualizeLine(entity,ref_path,{.thickness = 1.0,
                                          .type = 0,
                                          .color_r = 1.0,
                                          .color_g = 0.0,
                                          .color_b = 0.0,
                                          .color_a = 1.0});
        visualizeLine(entity,qp_final_path,{.thickness = 1.0,
                                          .type = 0,
                                          .color_r = 0.54,
                                          .color_g = 0.17,
                                          .color_b = 0.88,
                                          .color_a = 1.0});

        
        for(auto &obs:new_obs_array){
            visualizeAgent(entity,obs,{.width = 2,
                                       .length = 5,
                                       .color_r = 0.48,
                                       .color_g = 0.41,
                                       .color_b = 0.93,
                                       .color_a = 1.0
                                      });
        }

        visualizeAgent(entity,vehicle,{.width = 1.63,
                                       .length = 6,
                                       .color_r = 1.0,
                                       .color_g = 0.49,
                                       .color_b = 0.31,
                                       .color_a = 1.0
                                      });
        std::string scene_serializedMsg = scene.SerializeAsString();
        
        mcap::Message scene_msg;
        scene_msg.channelId = scene_channel;
        scene_msg.publishTime = curr_time;
        scene_msg.logTime = curr_time;
        scene_msg.data = reinterpret_cast<const std::byte*>(scene_serializedMsg.data());
        scene_msg.dataSize = scene_serializedMsg.size();
        writer.write(scene_msg);


        // get current vehicle Lateral distance
        int test_match_index = getMatchPoint(vehicle,-1,ref_path);
        waypoint test_match_point = ref_path.at(test_match_index);
        waypoint test_project_point;
        getProjectPoint(vehicle,test_match_point,test_project_point);
        frenet test_frenet;
        cardesian2Frenet(vehicle,test_match_index,test_project_point,ref_path,index_S,test_frenet);

        control = ref_match_index+150;
        vehicle = covertFromWaypoint(qp_final_path.at(0));
        pre_path = final_path;

        foxglove::Point2 vis_data_point;
        vis_data_point.set_x(static_cast<double>(control));
        vis_data_point.set_y(test_frenet.l);
        std::string data_serializedMsg = vis_data_point.SerializeAsString();
        
        mcap::Message data_msg;
        data_msg.channelId = data_channel;
        data_msg.publishTime = curr_time;
        data_msg.logTime = curr_time;
        data_msg.data = reinterpret_cast<const std::byte*>(data_serializedMsg.data());
        data_msg.dataSize = data_serializedMsg.size();
        writer.write(data_msg);
    }

    // After the while loop, close the writer
    writer.close();

    return 0;
}



