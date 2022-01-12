#include <ros/ros.h>
#include <iostream>
#include <nav_msgs/OccupancyGrid.h>
#include <vector>
#include <string>
#include <utility>
#include <geometry_msgs/PointStamped.h>
#include <ros_viz_tools/ros_viz_tools.h>
#include "reference_line_smoother.h"
#include "Eigen/Eigen"

using std::cout;
using std::endl;
using namespace Eigen;

struct Waypoint{
    Waypoint(const double x, const double y):
    x_(x),
    y_(y){}
    
    Waypoint() = default;

    double x_{};
    double y_{};
};

std::vector<Waypoint> way_point_list;   //原始路径

void waypointCb(const geometry_msgs::PointStampedConstPtr &p) {
    Waypoint way_point;
    way_point.x_ = p->point.x;
    way_point.y_ = p->point.y;
    way_point_list.emplace_back(way_point);
    cout << "way point received." << endl;
    // for(const auto& p: way_point_list){
    //     cout << "x: " << p.x_ << ", y: " << p.y_ << endl;
    // } 
}

int main(int argc, char **argv){

    ros::init(argc, argv, "map_pub_node");

    ros::NodeHandle nh;

    ros::Publisher map_pub = nh.advertise<nav_msgs::OccupancyGrid>("map_pub",1,true);

    ros::Subscriber way_point_sub = nh.subscribe("/clicked_point", 1, waypointCb);

    ros::Rate rate(10);

    nav_msgs::OccupancyGrid map_data;
    //fill map data
    map_data.header.stamp = ros::Time::now();
    map_data.header.frame_id = "map";

    map_data.info.map_load_time = ros::Time(0);
    map_data.info.resolution = 0.5;
    map_data.info.width = 200;
    map_data.info.height = 100;
    map_data.info.origin.position.x = -50;
    map_data.info.origin.position.y = -50;
    map_data.info.origin.position.z = 0;
    map_data.info.origin.orientation.x = 0;
    map_data.info.origin.orientation.y = 0;
    map_data.info.origin.orientation.z = 0;

    int data_temp[map_data.info.width * map_data.info.height] = {0};
    //test point
    data_temp[9] = 100;  
    int map_size = map_data.info.width * map_data.info.height;
    std::vector<signed char> a(data_temp, data_temp + map_size);
    map_data.data = a;
    
    ros_viz_tools::RosVizTools markers(nh, "markers");
    std::string frame_id = "/map";
    while(ros::ok()){ 
        int id = 0; 
        markers.clear();
        //publish markers
        visualization_msgs::Marker marker = markers.newSphereList(0.1, "reference point", id++, ros_viz_tools::RED, frame_id);
        for(const auto& p_tmp : way_point_list){
            geometry_msgs::Point p;
            p.x = p_tmp.x_;
            p.y = p_tmp.y_;
            p.z = 1.0;
            marker.points.push_back(p);   
        }
        markers.append(marker);

        //way_point_list转化为smoother输入
        std::vector<std::pair<double, double>> ref_points;
        for(const auto & waypoints : way_point_list){          
            ref_points.emplace_back(std::make_pair(waypoints.x_, waypoints.y_));
        }
        //设置bound
        std::vector<double> ref_bounds_x;
        std::vector<double> ref_bounds_y;
        for(int i = 0; i < ref_points.size(); ++i){
            ref_bounds_x.emplace_back(0.2);
            ref_bounds_y.emplace_back(0.2);
        }

        ReferenceLineSmoother smoother_test;
        smoother_test.set_ref_point(ref_points);
        smoother_test.set_ref_bound_x(ref_bounds_x);
        smoother_test.set_ref_bound_y(ref_bounds_y);
        smoother_test.Solve();

        std::vector<double> result_path_x = smoother_test.opt_x();
        std::vector<double> result_path_y = smoother_test.opt_y();
        //填充结果并可视化
        std::vector<Waypoint> result_path;      //平滑后的路径
        for(int i = 0; i < result_path_x.size(); ++i){
            result_path.emplace_back(Waypoint{result_path_x.at(i), result_path_y.at(i)});
        } 

        visualization_msgs::Marker result_path_marker =
            markers.newSphereList(0.1, "result path point", id++, ros_viz_tools::PINK, frame_id);
        for(const auto& p_tmp : result_path){
            geometry_msgs::Point p;
            p.x = p_tmp.x_;
            p.y = p_tmp.y_;
            p.z = 1.0;
            result_path_marker.points.push_back(p);   
        }
        markers.append(result_path_marker);

        visualization_msgs::Marker result_path_line =
            markers.newLineStrip(0.03, "result path line", id++, ros_viz_tools::YELLOW, frame_id);
        for(const auto& p_tmp : result_path){
            geometry_msgs::Point p;
            p.x = p_tmp.x_;
            p.y = p_tmp.y_;
            p.z = 1.0;
            result_path_line.points.push_back(p);   
        }
        markers.append(result_path_line);
        //publish markers
        markers.publish();

        //publish map
        map_pub.publish(map_data);

        ros::spinOnce();

        rate.sleep();

    }

}