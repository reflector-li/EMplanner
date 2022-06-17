//
// Created by 12971 on 2022/3/31.
//

#include <iomanip>
#include "frenetUtils.h"
#include <algorithm>

double NormalizeAngle(const double angle) {
    double a = std::fmod(angle + M_PI, 2.0 * M_PI);
    if (a < 0) {
        a += M_PI;
    } else {
        a -= M_PI;
    }
    return a;
}

void readTraje(std::vector<waypoint> &vecTraj,const std::string &fileName){
    std::ifstream infile(fileName,std::ios::in);
    if(!infile.is_open())
    {
        std::cout<<"can not open fine"<<std::endl;
        return;
    }
    waypoint Point;
    std::string line;
    std::stringstream ss;
    while(getline(infile,line))
    {
        ss<<line;
        ss>>Point.x>>Point.y>>Point.dirAngle>>Point.k;
        vecTraj.push_back(Point);
        ss.clear();
    }
    infile.close();
}

void readObs(std::vector<planPointInfo> &obs_array,const std::string &fileName)
{
    std::ifstream infile(fileName,std::ios::in);
    if(!infile.is_open())
    {
        std::cout<<"can not open fine"<<std::endl;
        return;
    }
    planPointInfo Point;
    std::string line;
    std::stringstream ss;
    while(getline(infile,line))
    {
        ss<<line;
        ss>>Point.x>>Point.y>>Point.dirAngle>>Point.k;
        obs_array.push_back(Point);
        ss.clear();
    }
    infile.close();
}

void writeTraje(std::vector<waypoint> &vecTraj, const std::string &fileName){
    std::ofstream out_file(fileName,std::ios::trunc);
    for(auto el:vecTraj)
    {
        out_file<<std::setprecision(8)<<el.x<<' '<<el.y<<' '<<el.dirAngle<<' '<<el.k<<std::endl;
    }
    out_file.close();
}

void writeEigen(const Eigen::VectorXd &vec,const std::string &file_name)
{
    std::ofstream out_file(file_name,std::ios::trunc);
    out_file<<vec;
    out_file.close();
}



void getDirAndK(std::vector<waypoint> &vecTraj, int index){
    if(index == 0 || index == vecTraj.size()-1)
    {
        return;
    }
    waypoint &p0 = vecTraj.at(index-1), &p1 = vecTraj.at(index),&p2 = vecTraj.at(index+1);
    double d =2*((p0.y-p2.y)*(p0.x-p1.x)-(p0.y-p1.y)*(p0.x-p2.x));
    if(fabs(d)<1e-8)
    {
        vecTraj.at(index).k = 0;
        double module = (p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y);
        std::vector<double> Line = {(p2.x - p1.x)/module,(p2.y - p1.y)/module};
        vecTraj.at(index).dirAngle = atan2(Line.at(1),Line.at(0));
        return;

    }
    double a = p0.y*p0.y - p1.y*p1.y + p0.x*p0.x-p1.x*p1.x;
    double b = p0.y*p0.y - p2.y*p2.y + p0.x*p0.x-p2.x*p2.x;

    double cx = ((p0.y - p2.y)*a - (p0.y - p1.y)*b)/d;
    double cy = ((p0.x - p2.x)*a - (p0.x - p1.x)*b)/(-d);

    double dx = cx - p1.x;
    double dy = cy - p1.y;
    double R = sqrt(dx*dx+dy*dy);
    vecTraj.at(index).k = 1/R;

    std::vector<double> differVector = {p2.x - p1.x,p2.y - p1.y};
    std::vector<double> dirVec = {-dy/R,dx/R};
    if((differVector.at(0)*dirVec.at(0)+ differVector.at(1)*dirVec.at(1))<0)
    {
        dirVec.at(0) = -dirVec.at(0);
        dirVec.at(1) = -dirVec.at(1);
    }


    vecTraj.at(index).dirAngle = atan2(dirVec.at(1),dirVec.at(0));
    if(index == 1)
    {
        dx = cx - p0.x;
        dy = cy - p0.y;
        R = sqrt(dx*dx+dy*dy);
        vecTraj.at(index-1).k = 1/R;
        differVector = {p1.x - p0.x,p1.y - p0.y};
        dirVec = {-dy/R,dx/R};
        if((differVector.at(0)*dirVec.at(0)+ differVector.at(1)*dirVec.at(1))<0)
        {
            dirVec.at(0) = -dirVec.at(0);
            dirVec.at(1) = -dirVec.at(1);
        }
        vecTraj.at(index-1).dirAngle = atan2(dirVec.at(1),dirVec.at(0));
    }
    if(index == vecTraj.size()-2)
    {
        dx = cx - p2.x;
        dy = cy - p2.y;
        R = sqrt(dx*dx+dy*dy);
        vecTraj.at(index+1).k = 1/R;
        differVector = {p2.x - p1.x,p2.y - p1.y};
        dirVec = {-dy/R,dx/R};
        if((differVector.at(0)*dirVec.at(0)+ differVector.at(1)*dirVec.at(1))<0)
        {
            dirVec.at(0) = -dirVec.at(0);
            dirVec.at(1) = -dirVec.at(1);
        }
        vecTraj.at(index+1).dirAngle = atan2(dirVec.at(1),dirVec.at(0));
    }
}

int getMatchPoint(planPointInfo hostPoint, int prePointIdx,const std::vector<waypoint> &vecTraj)
{
    if(prePointIdx == -1)
    {
        std::vector<double> distVec;
        for(int i = 0;i<vecTraj.size()-1;++i)
        {
            double tempDist = pow(hostPoint.x - vecTraj.at(i).x,2)+ pow(hostPoint.y - vecTraj.at(i).y,2);
            distVec.push_back(tempDist);
        }
        auto itr = std::min_element(distVec.begin(),distVec.end());
        return static_cast<int>(std::distance(distVec.begin(),itr));
    }
    else{

        double preDist = 0,dist = 0,nextDist = 0;
        int process = 0;
        int minIdx = 0;
        int dir = 0;

        // determine direction
        std::vector<double> vec = {hostPoint.x - vecTraj.at(prePointIdx).x, hostPoint.y - vecTraj.at(prePointIdx).y};
        std::vector<double> matchDir = {cos(vecTraj.at(prePointIdx).dirAngle), sin(vecTraj.at(prePointIdx).dirAngle)};
        double innerProd = vec.at(0)*matchDir.at(0)+vec.at(1)*vec.at(1);
        dir = innerProd > 0? 1:(innerProd < 0)?-1:0;

        // find match point forward or backward
        if(dir != 0)
        {
            int idx = prePointIdx;
            preDist = 1e10;
            dist = pow(hostPoint.x - vecTraj.at(prePointIdx).x,2)+ pow(hostPoint.y - vecTraj.at(prePointIdx).y,2);
            nextDist = pow(hostPoint.x - vecTraj.at(prePointIdx+dir).x,2)+ pow(hostPoint.y - vecTraj.at(prePointIdx+dir).y,2);
            while (idx >= 0 && idx <vecTraj.size())
            {
                preDist = dist;
                dist = pow(hostPoint.x - vecTraj.at(idx).x,2)+ pow(hostPoint.y - vecTraj.at(idx).y,2);
                if( dist<preDist)
                {
                    process = 0;
                    minIdx = idx;
                }
                else if(dist > preDist)
                {
                    process++;
                }
                else{}
                if (process >= 5) break;
                idx = idx + dir;
            }
        } else{
            minIdx = prePointIdx;
        }
        return minIdx;
    }
}

void getProjectPoint(planPointInfo hostPoint, const waypoint &matchPoint, waypoint &projectPoint)
{
    Eigen::Vector2d hostPos = {hostPoint.x,hostPoint.y},matchPointPos = {matchPoint.x,matchPoint.y};
    Eigen::Vector2d d = hostPos - matchPointPos;
    Eigen::Vector2d tao = {cos(matchPoint.dirAngle), sin(matchPoint.dirAngle)};
    Eigen::Vector2d projectPose = matchPointPos+tao*(d.transpose()*tao);
    projectPoint.k = matchPoint.k;
    projectPoint.dirAngle = matchPoint.dirAngle+matchPoint.k*static_cast<double>((d.transpose()*tao));
    projectPoint.x = projectPose[0];
    projectPoint.y = projectPose[1];
}

void getProjectPoint(const frenet &frenet_point,  const std::vector<waypoint> &ref_trajectory, const Eigen::VectorXd &index_S,waypoint &projectPoint)
{
    int match_index = -1;
    for(int i = 0;i<index_S.size()-1;++i)
    {
        if(index_S[i]<frenet_point.s && index_S[i+1]>frenet_point.s)
            match_index = i;
    }
    if(match_index == -1)
        return;
    waypoint match_point = ref_trajectory.at(match_index);
    Eigen::Vector2d r_n = {match_point.x,match_point.y};
    Eigen::Vector2d tao_n = {cos(match_point.dirAngle),sin(match_point.dirAngle)};
    double ds = frenet_point.s - index_S[match_index];

    Eigen::Vector2d r_r = r_n + ds*tao_n;
    projectPoint.x = r_r[0];
    projectPoint.y = r_r[1];
    projectPoint.dirAngle = match_point.dirAngle + ds*match_point.k;
    projectPoint.k = (match_point.k + ref_trajectory.at(match_index+1).k)/2;
}





void resampleTraje(double resample_interval, std::vector<waypoint> &vecTraj)
{
    std::vector<waypoint> original_trajectory(vecTraj);
    vecTraj.clear();
    vecTraj.push_back(original_trajectory.at(0));
    vecTraj.reserve(ceil(1.5*calcPathLength(original_trajectory)/resample_interval));
    for(int i = 1;i<original_trajectory.size();++i)
    {
        std::vector<waypoint> curve_points = {vecTraj.back(),
                                              original_trajectory.at(i),
                                              i<original_trajectory.size()-1? original_trajectory.at(i+1):original_trajectory.back()};
        calcCurvePara(curve_points);
        if(curve_points.at(1).k == 0)
        {
            resampleOnStraight(vecTraj,curve_points,resample_interval);
        }
        else{
            resampleOnCurve(vecTraj,curve_points,resample_interval);
        }
    }

}

double calcPathLength(const std::vector<waypoint> &vecTraj)
{
    double length = 0;
    for(int i = 0;i<vecTraj.size()-1;++i)
    {
        Eigen::Vector2d point_vector = {vecTraj.at(i+1).x - vecTraj.at(i).x,vecTraj.at(i+1).y - vecTraj.at(i).y};
        length += point_vector.norm();
    }
    return length;
}

void calcCurvePara(std::vector<waypoint> &curve_points){
    if(curve_points.size() != 3)
        return;
    waypoint &p0 = curve_points.at(0), &p1 = curve_points.at(1),&p2 = curve_points.at(2);
    double d =2*((p0.y-p2.y)*(p0.x-p1.x)-(p0.y-p1.y)*(p0.x-p2.x));
    if(fabs(d)<1e-8)
    {
        p1.k = 0;
        std::vector<double> Line = {(p1.x - p0.x),(p1.y - p0.y)};
        p1.dirAngle = atan2(Line.at(1),Line.at(0));
        return;

    }
    double a = p0.y*p0.y - p1.y*p1.y + p0.x*p0.x-p1.x*p1.x;
    double b = p0.y*p0.y - p2.y*p2.y + p0.x*p0.x-p2.x*p2.x;

    double cx = ((p0.y - p2.y)*a - (p0.y - p1.y)*b)/d;
    double cy = ((p0.x - p2.x)*a - (p0.x - p1.x)*b)/(-d);

    double dx = cx - p1.x;
    double dy = cy - p1.y;
    double R = sqrt(dx*dx+dy*dy);
    p1.k = 1/R;
    p0.k = 1/R;

    std::vector<double> differVector = {p2.x - p1.x,p2.y - p1.y};
    std::vector<double> dirVec = {-dy/R,dx/R};
    if((differVector.at(0)*dirVec.at(0)+ differVector.at(1)*dirVec.at(1))<0)
    {
        dirVec.at(0) = -dirVec.at(0);
        dirVec.at(1) = -dirVec.at(1);
    }
    p1.dirAngle = atan2(dirVec.at(1),dirVec.at(0));

    // calculate k and dir for p0
    dx = cx - p0.x;
    dy = cy - p0.y;
    differVector = {p1.x - p0.x,p1.y - p0.y};
    dirVec = {-dy/R,dx/R};
    if((differVector.at(0)*dirVec.at(0)+ differVector.at(1)*dirVec.at(1))<0)
    {
        dirVec.at(0) = -dirVec.at(0);
        dirVec.at(1) = -dirVec.at(1);
    }
    p0.dirAngle = atan2(dirVec.at(1),dirVec.at(0));
}

void resampleOnStraight(std::vector<waypoint> &vecTraj,std::vector<waypoint> &curve_points,double resample_interval){
    if(curve_points.size() != 3)
        return;
    waypoint prePoint(vecTraj.back());
    Eigen::Vector2d vecDiff = {curve_points.at(1).x - prePoint.x,curve_points.at(1).y - prePoint.y};
    double dist = vecDiff.norm();
    double coeff = resample_interval/dist;
    vecDiff[0] *= coeff;
    vecDiff[1] *= coeff;
    for(;dist>resample_interval;dist -= resample_interval)
    {
        prePoint.x += vecDiff[0];
        prePoint.y += vecDiff[1];
        vecTraj.push_back(prePoint);
    }
}

void resampleOnCurve(std::vector<waypoint> &vecTraj,std::vector<waypoint> &curve_points,double resample_interval)
{
    if(curve_points.size() != 3)
        return;
    waypoint prePoint(curve_points.at(0));
    double R = 1/curve_points.at(1).k;
    int dir = 0;

    // judge whither clockwise or anticlockwise
    Eigen::Vector2d p0 = {cos(prePoint.dirAngle), sin(prePoint.dirAngle)};
    Eigen::Vector2d p1 = {cos(curve_points.at(1).dirAngle), sin(curve_points.at(1).dirAngle)};

    double cross = p0[0]*p1[1]- p0[1]*p1[0];
    double dot = p0.dot(p1);
    double theta = acos(dot/(p0.norm()*p1.norm()));
    if(cross>0)
    {
        dir = 1;
    }
    else
    {
        dir = -1;//clockwise
    }

    double dist = fabs(theta)*R;
    double theta_diff = resample_interval*curve_points.at(1).k;
    for(;dist>resample_interval;dist -= resample_interval){

        if(vecTraj.size() == vecTraj.capacity())
            break;

        Eigen::Vector2d vec = {cos(prePoint.dirAngle+dir*theta_diff/2), sin(prePoint.dirAngle+dir*theta_diff/2)};
        vec = 2*R* sin(theta_diff/2)*vec;
        prePoint.dirAngle += dir*theta_diff;
        prePoint.x += vec[0];
        prePoint.y += vec[1];
        vecTraj.push_back(prePoint);
    }
}

// use 30 points before matchPoint and 150 points after it, so total frenet_trajectory size is 181
void getRefTrajectory(int matchPoint_index, const std::vector<waypoint> &vecTraj,std::vector<waypoint> &frenet_trajectory)
{

    frenet_trajectory.clear();
    frenet_trajectory.reserve(181);
    if (vecTraj.size()<181)
    {
        std::cout<<"the global waypoints is too short!"<<std::endl;
        return;
    }
    int startIndex = 0;
    if(matchPoint_index<29)
        startIndex = 0;
    else if(vecTraj.size() - matchPoint_index<151)
        startIndex = static_cast<int>(vecTraj.size()) - 181;
    else
        startIndex = matchPoint_index - 30;

    // use std::copy, it copys [start,end) to new vector.
    std::copy(vecTraj.begin()+startIndex,vecTraj.begin()+startIndex+181,std::back_inserter(frenet_trajectory));
}

double getDistance(const waypoint &p1,const waypoint &p2)
{
    return sqrt(pow(p1.x - p2.x,2)+ pow(p1.y-p2.y,2));
}

double calcSFromIndex2S(int matchPoint_index_local,const waypoint &project_point,const std::vector<waypoint> &frenet_trajectory,const Eigen::VectorXd &index_S)
{
    double heading = frenet_trajectory.at(matchPoint_index_local).dirAngle;
    Eigen::Vector2d tangent = {cos(heading), sin(heading)};
    Eigen::Vector2d dir = {project_point.x - frenet_trajectory.at(matchPoint_index_local).x,project_point.y - frenet_trajectory.at(matchPoint_index_local).y};
    double s_temp = getDistance(frenet_trajectory.at(matchPoint_index_local),project_point);
    double s0 = dir.dot(tangent) >0?index_S[matchPoint_index_local]+s_temp:index_S[matchPoint_index_local]-s_temp;
    return s0;
}

void index2S(int matchPoint_index_local,const waypoint &project_point,const std::vector<waypoint> &frenet_trajectory,Eigen::VectorXd &index_S)
{
    int size = static_cast<int>(frenet_trajectory.size());
    index_S.resize(size);
    index_S[0] = 0;
    for(int i = 1;i<frenet_trajectory.size();++i)
    {
        index_S[i] = getDistance(frenet_trajectory.at(i),frenet_trajectory.at(i-1))+index_S[i-1];
    }
    double s0 = calcSFromIndex2S(matchPoint_index_local,project_point,frenet_trajectory,index_S);
    index_S = index_S - s0*Eigen::VectorXd::Ones(index_S.size());
}


// if Longitudinal error is bigger than 1.5, or lateral error is bigger than 0.5, thus contorl ability is believed worse.
bool controlWeak(const planPointInfo &current_vehicle,const std::vector<waypoint> &pre_trajectory)
{
    int match_index = -1;
    for(int i = 0;i<pre_trajectory.size();i++)
    {
        if(pre_trajectory.at(i).time <= current_vehicle.time && pre_trajectory.at(i+1).time > current_vehicle.time ) {
            match_index = i;
            break;
        }
    }
    double heading = pre_trajectory.at(match_index).dirAngle;
    Eigen::Vector2d tangential = {cos(heading), sin(heading)};
    Eigen::Vector2d normal = {-sin(heading),cos(heading)};
    Eigen::Vector2d error_vec = {current_vehicle.x - pre_trajectory.at(match_index).x,current_vehicle.y - pre_trajectory.at(match_index).y};
    double lon_error = error_vec.dot(tangential); // Longitudinal error
    double lat_error = error_vec.dot(normal); // lateral error
    if(lon_error > 2.5 || lat_error > 0.5)
        return true;
    else
        return false;
}

// stitch 20 points of pre_trajectory to new trajectory except the start point.
void getStartPlanPoint(const std::vector<waypoint> &pre_trajectory,const planPointInfo &current_vehicle,planPointInfo &start_plan_point,std::vector<waypoint> &stitch_trajectory)
{
    double dt = 0.1;
    stitch_trajectory.clear();
    stitch_trajectory.reserve(20);
    if(pre_trajectory.empty())
    {
        start_plan_point.x = current_vehicle.x;
        start_plan_point.y = current_vehicle.y;
        start_plan_point.velocity_x = 0;
        start_plan_point.velocity_y = 0;
        start_plan_point.dirAngle = current_vehicle.dirAngle;
        start_plan_point.accel_x = 0;
        start_plan_point.accel_y = 0;
        start_plan_point.k  = 0;
        start_plan_point.time = current_vehicle.time+dt;
        return;
    }
    else
    {
        if(controlWeak(current_vehicle,pre_trajectory))
        {
            double vx_global = current_vehicle.velocity_x*cos(current_vehicle.dirAngle) - current_vehicle.velocity_y*sin(current_vehicle.dirAngle);
            double vy_global = current_vehicle.velocity_x*sin(current_vehicle.dirAngle) + current_vehicle.velocity_y*cos(current_vehicle.dirAngle);
            double ax_global = current_vehicle.accel_x*cos(current_vehicle.dirAngle) - current_vehicle.accel_y*sin(current_vehicle.dirAngle);
            double ay_global = current_vehicle.accel_x*sin(current_vehicle.dirAngle) + current_vehicle.accel_y*cos(current_vehicle.dirAngle);

            start_plan_point.x = current_vehicle.x + vx_global*dt+0.5*ax_global*dt*dt;
            start_plan_point.y = current_vehicle.y + vy_global*dt+0.5*ay_global*dt*dt;
            start_plan_point.velocity_x = vx_global+ax_global*dt;
            start_plan_point.velocity_y = vy_global+ay_global*dt;

            start_plan_point.dirAngle = atan2(start_plan_point.velocity_y,start_plan_point.velocity_x);

            start_plan_point.accel_x = current_vehicle.accel_x;
            start_plan_point.accel_y = current_vehicle.accel_y;
            start_plan_point.k  = 0;
            start_plan_point.time = current_vehicle.time+dt;
            return;

        } else
        {
            int expert_Index = -1;
            for(int i = 0;i<pre_trajectory.size();++i)
            {
                if(pre_trajectory.at(i).time<=(current_vehicle.time+dt) && pre_trajectory.at(i+1).time>(current_vehicle.time+dt))
                {
                    expert_Index = i;
                    break;
                }

            }
            start_plan_point.x = pre_trajectory.at(expert_Index).x;
            start_plan_point.y = pre_trajectory.at(expert_Index).y;
            start_plan_point.dirAngle = pre_trajectory.at(expert_Index).dirAngle;
            start_plan_point.velocity_x = pre_trajectory.at(expert_Index).velocity*cos(start_plan_point.dirAngle);
            start_plan_point.velocity_y = pre_trajectory.at(expert_Index).velocity*sin(start_plan_point.dirAngle);
            start_plan_point.k  = pre_trajectory.at(expert_Index).k;

            //converting waypoint acceleration to ax_global and ay_global needs tangential acceleration and normal acceleration
            Eigen::Vector2d tangentail = {cos(start_plan_point.dirAngle), sin(start_plan_point.dirAngle)};
            Eigen::Vector2d normal = {-sin(start_plan_point.dirAngle),cos(start_plan_point.dirAngle)};
            Eigen::Vector2d tangential_accel = pre_trajectory.at(expert_Index).accel*tangentail;
            Eigen::Vector2d normal_accel = pow(pre_trajectory.at(expert_Index).velocity,2)*start_plan_point.k*normal;

            start_plan_point.accel_x = (tangential_accel+normal_accel)[0];
            start_plan_point.accel_y = (tangential_accel+normal_accel)[1];

            start_plan_point.time = pre_trajectory.at(expert_Index).time;

            //stitch 20 points of pre_trajectory to current period plan
            if (expert_Index>=20)
                std::copy(pre_trajectory.begin()+expert_Index-20,pre_trajectory.begin()+expert_Index,std::back_inserter(stitch_trajectory));
            else
                std::copy(pre_trajectory.begin(),pre_trajectory.begin()+expert_Index,std::back_inserter(stitch_trajectory));
            return;
        }
    }
    
    
}

void cardesian2Frenet(const planPointInfo &cardesian_point, int match_point_index, const waypoint &project_point, const std::vector<waypoint> &ref_trajectory,const Eigen::VectorXd &index_S,frenet &frenet_pos,cardesian2FrentMode mode)
{
    if(ref_trajectory.empty())
    {
        std::cout<<"reference line is empty!"<<std::endl;
        return;
    }
    Eigen::Vector2d r_h = {cardesian_point.x,cardesian_point.y};
    Eigen::Vector2d r_r = {project_point.x,project_point.y};
    Eigen::Vector2d tao_r = {cos(project_point.dirAngle),sin(project_point.dirAngle)};
    Eigen::Vector2d n_r = {-sin(project_point.dirAngle), cos(project_point.dirAngle)};
    Eigen::Vector2d n_h = {-sin(cardesian_point.dirAngle), cos(cardesian_point.dirAngle)};


    frenet_pos.s = calcSFromIndex2S(match_point_index,project_point,ref_trajectory,index_S);

    frenet_pos.l = (r_h - r_r).dot(n_r);

    if(mode == firstInfo || mode == secondInfo)
    {
        Eigen::Vector2d v = {cardesian_point.velocity_x,cardesian_point.velocity_y};

        frenet_pos.l_dot = v.dot(n_r);
        frenet_pos.s_dot = (1/(1-frenet_pos.l*project_point.k))* v.dot(tao_r);
        if(fabs(frenet_pos.s_dot)<1e-8)
            frenet_pos.l_diff = 0;
        else
            frenet_pos.l_diff = frenet_pos.l_dot/ frenet_pos.s_dot;
        if(mode == secondInfo)
        {
            Eigen::Vector2d a = {cardesian_point.accel_x,cardesian_point.accel_y};

            frenet_pos.l_2dot = a.dot(n_r) - project_point.k*(1-project_point.k*frenet_pos.l)*pow(frenet_pos.s_dot,2);
            // because d_k/ds is small, use 0 replace
            frenet_pos.s_2dot =(1/(1-frenet_pos.l*project_point.k))*(a.dot(tao_r)+2*project_point.k*frenet_pos.l_diff*pow(frenet_pos.s_dot,2));
            if(fabs(frenet_pos.s_dot)<1e-8)
                frenet_pos.l_2diff = 0;
            else
                frenet_pos.l_2diff = (1/pow(frenet_pos.s_dot,2))*(frenet_pos.l_2dot - frenet_pos.l_diff*frenet_pos.s_2dot);
        }
        return;
    }
}

// only consider obstacles with lontitude error within [-10,60], lateral error within [-10,10]
void obsSelect(const planPointInfo &currentvehicle,const std::vector<planPointInfo> &obstacle_array_original, std::vector<planPointInfo> &obstacle_array_selected)
{
    obstacle_array_selected.clear();
    obstacle_array_selected.reserve(32);
    Eigen::Vector2d tor = {cos(currentvehicle.dirAngle),sin(currentvehicle.dirAngle)};
    Eigen::Vector2d nor = {-sin(currentvehicle.dirAngle),cos(currentvehicle.dirAngle)};
    for(const auto &obstacle : obstacle_array_original)
    {
        Eigen::Vector2d pose_error = {obstacle.x - currentvehicle.x,obstacle.y - currentvehicle.y};
        double lon_error = pose_error.dot(tor);
        double lat_error = pose_error.dot(nor);
        if((lon_error > -10 && lon_error<60) && (lat_error > -10 && lat_error<10))
            obstacle_array_selected.push_back(obstacle);
    }
}

// using dp algorithm to calculate the init path.

double calcNeighborCost(const std::vector<frenet> &obstacle_array,const frenet &start_point, const frenet &end_point, const initPlanConfigure &config)
{
    Eigen::VectorXd coeff(6);
    calcQuinticCoeff(start_point,end_point,coeff);
    
    //use 10 points to calculate the cost
    double cost_sum = 0;
    for(int i = 0;i<10;++i)
    {
        double ds = start_point.s + i*config.sample_s/10;
        double ds_2 = ds*ds, ds_3 = ds_2*ds, ds_4 = ds_3*ds, ds_5 = ds_4*ds;
        double l = coeff[0] + coeff[1]*ds + coeff[2]*ds_2 + coeff[3]*ds_3 + coeff[4]*ds_4 + coeff[5]*ds_5;
        double l_diff = coeff[1] + 2*coeff[2]*ds + 3*coeff[3]*ds_2 + 4*coeff[4]*ds_3 + 5*coeff[5]*ds_4;
        double l_2diff = 2*coeff[2] + 6*coeff[3]*ds + 12*coeff[4]*ds_2 + 20*coeff[5]*ds_3;
        double l_3diff = 6*coeff[3] + 24*coeff[4]*ds + 60*coeff[5]*ds_2;
        
        double cost_smooth = config.w_cost_dl*pow(l_diff,2) + config.w_cost_ddl*pow(l_2diff,2)+config.w_cost_dddl*pow(l_3diff,2);
        double cost_ref = config.w_cost_ref* pow(l,2);
        double cost_collision = 0;
        for(const auto &obstacle:obstacle_array)
        {
            Eigen::Vector2d pose_error = {obstacle.s - ds, obstacle.l - l};
            double square_dist = pose_error.dot(pose_error); // Not Euclidean distance, just approximate
            double cost_collision_once = calcObsCost(config.w_cost_obs,square_dist);
            cost_collision += cost_collision_once;
        }

        cost_sum = cost_sum + cost_smooth + cost_ref + cost_collision;
    }

}


double calcObsCost(double w_cost_obs, double square_dist)
{
    // if dist>4 cost = 0; if dist in [3,4], cost = 1000/square_dist; if dist<3, cost = w_cost_obs
    double cost = 0;
    if(square_dist >= 9 && square_dist <= 16)
        cost = 1000/square_dist;
    else if(square_dist < 9)
        cost = w_cost_obs;
    else
        cost = 0;
    return cost;

}

void calcQuinticCoeff(const frenet &start_point, const frenet &end_point, Eigen::VectorXd &coeff)
{
    // l = a0 + a1*s +a2*s^2 + a3*s^3 + a4*s^4 + a5*s^5;
    // l' = a1 + 2*a2*s + 3*a3*s^2 + 4*a4*s^3 + 5*a5*s^4
    // l'' = 2*a2 + 6*a3*s + 12*a4*s^2 + 20*a5*s^3

    Eigen::Matrix<double,6,6> A;
    Eigen::VectorXd B(6);
    double S_2 = start_point.s * start_point.s;
    double S_3 = S_2 * start_point.s;
    double S_4 = S_3 * start_point.s;
    double S_5 = S_4 * start_point.s;

    double E_2 = end_point.s * end_point.s;
    double E_3 = E_2 * end_point.s;
    double E_4 = E_3 * end_point.s;
    double E_5 = E_4 * end_point.s;

    A<<1, start_point.s,   S_2,   S_3, S_4, S_5,
       0,  1, 2*start_point.s, 3*S_2,  4*S_3, 5*S_4,
       0, 0, 2, 6*start_point.s, 12*S_2, 20*S_3,
       1, end_point.s,   E_2,   E_3, E_4, E_5,
       0, 1, 2*end_point.s, 3*E_2,  4*E_3, 5*E_4,
       0, 0, 2, 6*end_point.s, 12*E_2, 20*E_3;
    
    B<<start_point.l, start_point.l_diff,start_point.l_2diff,end_point.l,end_point.l_diff,end_point.l_2diff;
    coeff = A.inverse()*B;
}


void getInitPlanValue(const std::vector<frenet> &obstacle_array,const frenet &start_point, const initPlanConfigure &config, std::vector<frenet> &init_path)
{
    init_path.clear();
    std::vector<std::vector<costElement> > cost_matrix;
    cost_matrix.resize(config.rows);
    for(int i = 0;i<config.rows;++i)
    {
        cost_matrix[i].resize(config.cols);
    }
    for(int i = 0;i<config.rows;i++)
    {
        double l = (config.rows/2 - i)*config.sample_l;
        for(int j = 0;j<config.cols;j++)
        {
            cost_matrix[i][j].pose.l = l;
            cost_matrix[i][j].pose.s = start_point.s+(j+1)*config.sample_s;
            cost_matrix[i][j].pose.l_diff = 0;
            cost_matrix[i][j].pose.l_2diff = 0;
        }
    }

    // calculate the first column cost
    for(int i = 0;i<config.rows;++i)
    {
        double cost = calcNeighborCost(obstacle_array,start_point,cost_matrix[i][0].pose,config);
        cost_matrix[i][0].miniCost = cost;
    }

    for(int j = 1;j<config.cols;++j)
    {
        for(int i = 0;i<config.rows;++i)
        {
            for(int k = 0;k<config.rows;++k)
            {
                double cost_temp = calcNeighborCost(obstacle_array,cost_matrix[k][j-1].pose,cost_matrix[i][j].pose,config);
                double pre_mini_cost = cost_matrix[k][j-1].miniCost;
                double cost_cur = pre_mini_cost+cost_temp;
                if(cost_cur < cost_matrix[i][j].miniCost)
                {
                    cost_matrix[i][j].miniCost = cost_cur;
                    cost_matrix[i][j].preRow = k;
                    cost_matrix[i][j].preCol = j-1;
                }
            }
        }
    }

    int mini_cost_row = 0;
    int mini_cost_col = config.cols -1;
    double mini_cost_last = 1e8;
    for(int i = 0;i<config.rows;++i)
    {

        if(cost_matrix[i][mini_cost_col].miniCost < mini_cost_last)
        {
            mini_cost_last = cost_matrix[i][mini_cost_col].miniCost;
            mini_cost_row = i;
        }
    }

    init_path.push_back(cost_matrix[mini_cost_row][mini_cost_col].pose);
    for(int i = 0;i<config.cols-1;++i)
    {
        mini_cost_row = cost_matrix[mini_cost_row][mini_cost_col].preRow;
        mini_cost_col = cost_matrix[mini_cost_row][mini_cost_col].preCol;
        init_path.push_back(cost_matrix[mini_cost_row][mini_cost_col].pose);
    }
    init_path.push_back(start_point);
    std::reverse(init_path.begin(),init_path.end());
}

// interpolate the init_path, point is taken every 1m of s, maximum 60m
void trajectoryInterp(std::vector<frenet> &init_path)
{
    if (init_path.empty())
        return;
    std::vector<frenet> ori_path = init_path;
    init_path.clear();
    init_path.reserve(61);
    double ds = 1;
    init_path.push_back(ori_path.at(0));
    double s_cur = ori_path.at(0).s + ds;
    int count = 1;
    for(int i = 0;i<ori_path.size()-1;++i)
    {
        Eigen::VectorXd coeff(6);
        frenet p1 = ori_path.at(i), p2 = ori_path.at(i+1);
        calcQuinticCoeff(p1,p2,coeff);
        while(s_cur<ori_path.at(i+1).s)
        {
            frenet point_curr;
            point_curr.s = s_cur;
            double s_2 = s_cur*s_cur , s_3 = s_2*s_cur , s_4 = s_3*s_cur , s_5 = s_4*s_cur ;
            point_curr.l= coeff[0] + coeff[1]*s_cur + coeff[2]*s_2 + coeff[3]*s_3 + coeff[4]*s_4 + coeff[5]*s_5;
            point_curr.l_diff = coeff[1] + 2*coeff[2]*s_cur + 3*coeff[3]*s_2 + 4*coeff[4]*s_3 + 5*coeff[5]*s_4;
            point_curr.l_2diff = 2*coeff[2] + 6*coeff[3]*s_cur + 12*coeff[4]*s_2 + 20*coeff[5]*s_3;
            init_path.push_back(point_curr);
            s_cur += ds;
            count++;
            if(count>61)
                return;
        }
    }
}

// dk_r/ds is thought to be 0
void frenet2Cardesian(const std::vector<frenet> &dense_init_path, const std::vector<waypoint> &ref_trajectory,const Eigen::VectorXd &index_S, std::vector<waypoint> &init_path_cardesian)
{
    if(dense_init_path.empty())
        return;
    init_path_cardesian.clear();
    init_path_cardesian.reserve(dense_init_path.size());
    for(const auto &frenet_point:dense_init_path)
    {
        waypoint project_point,cur_point;
        getProjectPoint(frenet_point,ref_trajectory,index_S,project_point);
        Eigen::Vector2d r_r = {project_point.x,project_point.y};
        Eigen::Vector2d tau_r = {cos(project_point.dirAngle),sin(project_point.dirAngle)};
        Eigen::Vector2d n_r = {-sin(project_point.dirAngle), cos(project_point.dirAngle)};

        Eigen::Vector2d r_h = r_r + frenet_point.l*n_r;
        cur_point.x = r_h[0];
        cur_point.y = r_h[1];
        cur_point.dirAngle = project_point.dirAngle+ atan2(frenet_point.l_diff,(1-frenet_point.l*project_point.k));
        cur_point.dirAngle = NormalizeAngle(cur_point.dirAngle);

        double d_theta = cur_point.dirAngle - project_point.dirAngle;
        cur_point.k = ((frenet_point.l_2diff+project_point.k*frenet_point.l_diff* tan(d_theta))*
                      pow(cos(d_theta),2)*(1/(1 - project_point.k*frenet_point.l))+project_point.k)*
                              cos(d_theta)*(1/(1 - project_point.k*frenet_point.l));
        init_path_cardesian.push_back(cur_point);
    }
}

void timeAndVelocity(std::vector<waypoint> &init_path_cardesian, double startPoint_time)
{
    if(init_path_cardesian.empty())
        return;
    int count = 0;
    for(auto &waypoint_:init_path_cardesian)
    {
        waypoint_.time = startPoint_time + count*0.02;
        count++;
    }
}

waypoint covertFromPlanPointInfo(const planPointInfo &point)
{
    waypoint result;
    result.x = point.x;
    result.y = point.y;
    result.dirAngle = point.dirAngle;
    result.k = point.k;
    result.time = point.time;
    return result;
}

planPointInfo covertFromWaypoint(const waypoint &point)
{
    planPointInfo result;
    result.x = point.x;
    result.y = point.y;
    result.dirAngle = point.dirAngle;
    result.k = point.k;
    result.time = point.time;
    return result;
}

// obstacle_length = 5, obstacle_width = 2
void getCovexBound(const std::vector<frenet> &init_path,const std::vector<frenet> &obstacle_array,const qpPlanConfigure &config,Eigen::VectorXd &low_bound, Eigen::VectorXd &upper_bound)
{
    if(init_path.empty())
        return;
    low_bound = -6*Eigen::VectorXd::Ones(60); //right bound of vehicle
    upper_bound = 6*Eigen::VectorXd::Ones(60); //left bound of vehicle
    for(auto &obs:obstacle_array)
    {
        double obs_s_min = obs.s - config.obs_length/2.0;
        double obs_s_max = obs.s + config.obs_length/2.0;
        int index_s_min = getMatchSIndex(init_path,obs_s_min);
        int index_s_max = getMatchSIndex(init_path,obs_s_max);
        int index_obs = getMatchSIndex(init_path,obs.s);
        if(index_s_min==0 || index_s_max == 0)
            continue;
        if(init_path.at(index_obs).l>obs.l) // dp-path on left of obstacles
        {
            for(int i = index_s_min;i<=index_s_max;++i)
                low_bound[i] =std::max(obs.l + config.obs_width/2.0,low_bound[i]);
        }
        else // dp-path on right of obstacles
        {
            for(int i = index_s_min;i<=index_s_max;++i)
                upper_bound[i] = std::min(obs.l - config.obs_width/2.0,upper_bound[i]);
        }
    }
}

int getMatchSIndex(const std::vector<frenet> &init_path,double s)
{
    int count = 0,index = -1;
    if(s<init_path.at(0).s)
        return index = 0;
    else if(s>init_path.back().s)
        return index = 59;
    else
    {
        for(auto &point:init_path)
        {
            if(point.s>s)
                break;
            count++;
        }
        if((init_path.at(count).s - s) > (s - init_path.at(count-1).s))
            return index = count;
        else
            return index = count-1;
    }
}



void createEdge(const std::vector<waypoint> &global_path, std::vector<waypoint> &upper_edge, std::vector<waypoint> &low_edge)
{
    upper_edge.reserve(global_path.size());
    low_edge.reserve(global_path.size());

    for(const auto & point : global_path)
    {
        waypoint temp_point;
        Eigen::Vector2d pos = {point.x,point.y};
        Eigen::Vector2d tor = {-sin(point.dirAngle), cos(point.dirAngle)};
        Eigen::Vector2d upper_point = pos+6*tor;
        temp_point.x = upper_point[0];
        temp_point.y = upper_point[1];
        temp_point.dirAngle = point.dirAngle;
        temp_point.k = point.k;
        upper_edge.push_back(temp_point);

        upper_point = pos-6*tor;
        temp_point.x = upper_point[0];
        temp_point.y = upper_point[1];
        temp_point.dirAngle = point.dirAngle;
        temp_point.k = point.k;
        low_edge.push_back(temp_point);
    }
}

