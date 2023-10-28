#include "visualization.h"

uint64_t nanosecondsSinceEpoch(){
    return uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::system_clock::now().time_since_epoch())
                    .count());
}

void EulerToQuaternion(foxglove::Quaternion* q, double yaw, double pitch, double roll){
    double cy = cos(yaw * 0.5);
    double sy = sin(yaw * 0.5);
    double cp = cos(pitch * 0.5);
    double sp = sin(pitch * 0.5);
    double cr = cos(roll * 0.5);
    double sr = sin(roll * 0.5);
 
    q->set_w(cy * cp * cr + sy * sp * sr);
    q->set_x(cy * cp * sr - sy * sp * cr);
    q->set_y(sy * cp * sr + cy * sp * cr);
    q->set_z(sy * cp * cr - cy * sp * sr);
}

void setAxisAngle(foxglove::Quaternion* q, double x, double y, double z, double angle) {
  double s = std::sin(angle / 2);
  q->set_x(x * s);
  q->set_y(y * s);
  q->set_z(z * s);
  q->set_w(std::cos(angle / 2));
}


void visualizeLine(foxglove::SceneEntity *entity,const std::vector<waypoint> &path,const LineVisConfig &config){
    auto *line = entity->add_lines();
    switch (config.type)
    {
    case 0:
        line->set_type(foxglove::LinePrimitive_Type::LinePrimitive_Type_LINE_STRIP);
        break;
    case 1:
        line->set_type(foxglove::LinePrimitive_Type::LinePrimitive_Type_LINE_LOOP);
        break;
    case 2:
        line->set_type(foxglove::LinePrimitive_Type::LinePrimitive_Type_LINE_LIST);
        break;
    default:
        line->set_type(foxglove::LinePrimitive_Type::LinePrimitive_Type_LINE_STRIP);
        break;
    }
    line->set_thickness(config.thickness);
    line->set_scale_invariant(true);
    auto color = line->mutable_color();
    color->set_r(config.color_r);
    color->set_g(config.color_g);
    color->set_b(config.color_b);
    color->set_a(config.color_a);
    for(auto point:path){
        auto *foxglove_point = line->add_points();
        foxglove_point->set_x(point.x);
        foxglove_point->set_y(point.y);
        foxglove_point->set_z(0);
    }
}

void visualizeAgent(foxglove::SceneEntity *entity,const planPointInfo &pose, const AgentVisConfig &config){
    auto *cube = entity->add_cubes();
    auto *size = cube->mutable_size();
    size->set_x(config.length);
    size->set_y(config.width);
    size->set_z(1); 
    auto* position = cube->mutable_pose()->mutable_position();
    position->set_x(pose.x);
    position->set_y(pose.y);
    position->set_z(0);
    auto *orientation = cube->mutable_pose()->mutable_orientation();
    EulerToQuaternion(orientation,pose.dirAngle,0,0);
    auto* color = cube->mutable_color();
    color->set_r(config.color_r);
    color->set_g(config.color_g);
    color->set_b(config.color_b);
    color->set_a(config.color_a);
}


void RootToEgoTF(foxglove::FrameTransform &tf,std::string parent_frame,std::string child_frame,const planPointInfo &pose){
    tf.set_parent_frame_id(parent_frame);
    tf.set_child_frame_id(child_frame);
    auto *position = tf.mutable_translation();
    position->set_x(pose.x);
    position->set_y(pose.y);
    position->set_z(0);
    auto *orientation = tf.mutable_rotation();
    EulerToQuaternion(orientation,pose.dirAngle,0,0);
}

