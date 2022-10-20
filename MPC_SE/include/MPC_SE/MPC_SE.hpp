#pragma once

#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <cmath>
#include <list>

#include <eigen3/Eigen/Dense>

#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/twist.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "tf2_msgs/msg/tf_message.hpp"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

#include <vector>

#include <omp.h> // multi-core

#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>

static size_t N = 6;
static double dt = 0.1;

static int node_number=2;

using namespace std::chrono_literals;
USING_NAMESPACE_ACADO

class MPC_SE : public rclcpp::Node
{
public:
    MPC_SE(): Node("MPC_SE")
    {
        node_name = this->get_name();
        node_name = "Vehicle1";
        robot_num = node_name[7];

        for(int i = 0; i < node_number; i++)
        {
            std::string ads = node_name;
            ads[7] = '1' + i;
            odom_sub_list[i] = this->create_subscription<nav_msgs::msg::Odometry>(ads+"/odom/unfiltered", 10, [this, i, ads](const nav_msgs::msg::Odometry::SharedPtr msg){odom[i] = *msg; if(ads == this->node_name) sub_odom_callback(odom[i]);});
        }

        this->cmd_vel = this->create_publisher<geometry_msgs::msg::Twist>(node_name+"/cmd_vel", 10);
        this->timer = this->create_wall_timer(100ms, std::bind(&MPC_SE::timer_callback, this));

        f << dot(x) == cos(theta)*vx - sin(theta)*vy;
        f << dot(y) == sin(theta)*vx + cos(theta)*vy;
        f << dot(theta) == w;

        h << x-0.5;
        h << y-0.5;
        h << theta;

        Q.resize(3, 3);
        Q.setIdentity();

        r.resize(3);
        r.setAll(0.0);
    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    std::string node_name;
    char robot_num;

    DifferentialState x, y, theta;
    Control vx, vy, w;
    DifferentialEquation f;
    Function h;

    DMatrix Q;
    DVector r;

    double x0, y0, theta0;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x0 = odom1.pose.pose.position.x;
        y0 = odom1.pose.pose.position.y;
        theta0 = Orientation2Elur(odom1);
    }
    double Orientation2Elur(nav_msgs::msg::Odometry ori_odom)
    {
        double roll, pitch, yaw;
        tf2::Quaternion imu(ori_odom.pose.pose.orientation.x, ori_odom.pose.pose.orientation.y, ori_odom.pose.pose.orientation.z, ori_odom.pose.pose.orientation.w);
        tf2::Matrix3x3 m(imu);
        m.getRPY(roll, pitch, yaw);
        return yaw;
    }
    void timer_callback()
    {
        OCP ocp;

        ocp.minimizeLSQ(Q, h, r);
        ocp.subjectTo(f);
        ocp.subjectTo(AT_START, x == x0);
        ocp.subjectTo(AT_START, y == y0);
        ocp.subjectTo(AT_START, theta == theta0);
        ocp.subjectTo(-1.0<=vx<=1.0);
        ocp.subjectTo(-0.0<=vy<=0.0);
        ocp.subjectTo(-0.3<=w<=0.3);

        OptimizationAlgorithm algorithm(ocp);
        algorithm.solve();

        VariablesGrid controls;

        algorithm.getControls(controls);

        DMatrix U = controls.getMatrix(0);

        geometry_msgs::msg::Twist twist;
        twist.linear.x = U(0);
        twist.linear.y = U(1);
        twist.angular.z = U(2);

        cmd_vel->publish(twist);
    }
    /*std::vector<double> Solve()
    {

    }*/
};
