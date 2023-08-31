#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include "foxglove/SceneUpdate.pb.h"
#include "foxglove/FrameTransform.pb.h"
#include "server_config.h"
#include "frenetUtils.h"
#include <string>
#include <vector>
#include <chrono>
#include <cmath>

struct LineVisConfig
{
    double thickness;
    int type; // 0 LINE_STRIP, 1 LINE_LOOP, 2 LINE_LIST
    double color_r;
    double color_g;
    double color_b;
    double color_a;
};

struct AgentVisConfig{
    double width;
    double length;
    double color_r;
    double color_g;
    double color_b;
    double color_a;
};


uint64_t nanosecondsSinceEpoch();
void EulerToQuaternion(foxglove::Quaternion* q, double yaw, double pitch, double roll);
void setAxisAngle(foxglove::Quaternion* q, double x, double y, double z, double angle);
void visualizeLine(foxglove::SceneEntity *entity,const std::vector<waypoint> &path,const LineVisConfig &config);
void visualizeAgent(foxglove::SceneEntity *entity,const planPointInfo &pose, const AgentVisConfig &config);
void RootToEgoTF(foxglove::FrameTransform &tf, std::string parent_frame,std::string child_frame,const planPointInfo &pose);

#endif