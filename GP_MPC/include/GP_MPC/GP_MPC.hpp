#pragma once

#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <cmath>

#include "gurobi_c++.h"
#include <eigen3/Eigen/Dense>

#include "gp.h"
#include "gp_utils.h"

#include "GP_MPC/gnuplot_i.hpp"

#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/twist.hpp"
#include "geometry_msgs/msg/quaternion.hpp"
#include "nav_msgs/msg/odometry.hpp"
#include "tf2_msgs/msg/tf_message.hpp"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.h>

#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <Eigen/Core>
#include <vector>

#include <manif/manif.h>

static size_t N = 6;
static double dt = 0.1;

static int node_number=1;

using namespace std::chrono_literals;
using CppAD::AD;
#if 0
class FG_eval // a small test
{
public:
    FG_eval(double x=0.0, double y=0.0)
    {
        x_ = x;
        y_ = y;
    }
    void setXY(double x, double y)
    {
        x_ = x;
        y_ = y;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < N-1; i++)
        {
            fg[0] += 10.0*(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2));
        }
        for (int i = 0; i < N-1; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2));
        }
        fg[0] += 100.0*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2) + CppAD::pow(state[3*N-1]-1.57, 2)); // The finnal point...

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];

        for (int i = 0; i < N-1; i++) // constraints
        {
            AD<double> x1 = state[i+1];
            AD<double> y1 = state[N+i+1];
            AD<double> theta1 = state[2*N+i+1];

            AD<double> x0 = state[i];
            AD<double> y0 = state[N+i];
            AD<double> theta0 = state[2*N+i];
            AD<double> v0 = state[3*N+i];
            AD<double> w0 = state[4*N+i];

            AD<double> kx1 = v0 * cos(theta0);
            AD<double> ky1 = v0 * sin(theta0);
            AD<double> kw1 = w0;

            AD<double> kx2 = v0 * cos(theta0 + dt/2*kw1);
            AD<double> ky2 = v0 * sin(theta0 + dt/2*kw1);
            AD<double> kw2 = w0;

            AD<double> kx3 = v0 * cos(theta0 + dt/2*kw2);
            AD<double> ky3 = v0 * sin(theta0 + dt/2*kw2);
            AD<double> kw3 = w0;

            AD<double> kx4 = v0 * cos(theta0 + dt*kw3);
            AD<double> ky4 = v0 * sin(theta0 + dt*kw3);
            AD<double> kw4 = w0;

            fg[2 + i]    = x1 - (x0 + dt*(kx1+2*kx2+2*kx3+kx4)/6);
            fg[2 + N + i]    = y1 - (y0 + dt*(ky1+2*ky2+2*ky3+ky4)/6);
            fg[2 + 2*N + i]  = theta1 - (theta0 + dt*(kw1+2*kw2+2*kw3+kw4)/6);

            fg[2 + 3*N + i] = (CppAD::pow(x1+2, 2) + CppAD::pow(y1+2.25, 2) - 1.5*1.5) - (1-1)*(CppAD::pow(x0+2, 2) + CppAD::pow(y0+2.25, 2) - 1.5*1.5);

        }
        //fg[5*N+1] = CppAD::pow(state[N-1]-5, 2) + CppAD::pow(state[2*N-1]-5, 2);

    }
private:
    double x_, y_;
};

class GP_MPC : public rclcpp::Node
{
public:
    GP_MPC(): Node("GP_MPC")
    {
        this->odom_sub = this->create_subscription<nav_msgs::msg::Odometry>("Vehicle1/odomg", 10, std::bind(&GP_MPC::sub_odom_callback, this, std::placeholders::_1));
        this->cmd_vel = this->create_publisher<geometry_msgs::msg::Twist>("Vehicle1/cmd_vel", 10);
        this->timer = this->create_wall_timer(100ms, std::bind(&GP_MPC::timer_callback, this));

        train_times = 0;
    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub;
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom;
    rclcpp::TimerBase::SharedPtr timer;

    FG_eval fg_eval;

    std::vector<double> tss, x_train, y_train, x_test, y_test;

    double x, y, z, theta, v, w;

    int train_times;



    void sub_odom_callback(nav_msgs::msg::Odometry::SharedPtr _odom)
    {
        this->odom = *_odom;
        x = odom.pose.pose.position.x;
        y = odom.pose.pose.position.y;
        z = odom.pose.pose.position.z;
        theta = Orientation2Elur(odom);
        v = odom.twist.twist.linear.x;
        w = odom.twist.twist.angular.z;
    }
    void timer_callback()
    {
        std::vector<double> cmd = Solve();
        geometry_msgs::msg::Twist twist;
        twist.linear.x = cmd[0];
        twist.angular.z = cmd[1];
        cmd_vel->publish(twist);
    }
    double Orientation2Elur(nav_msgs::msg::Odometry ori_odom)
    {
        double roll, pitch, yaw;
        tf2::Quaternion imu(ori_odom.pose.pose.orientation.x, ori_odom.pose.pose.orientation.y, ori_odom.pose.pose.orientation.z, ori_odom.pose.pose.orientation.w);
        tf2::Matrix3x3 m(imu);
        m.getRPY(roll, pitch, yaw);
        return yaw;
    }
    std::vector<double> Solve()
    {


        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*6;

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);

        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = x;
        state[N] = y;
        state[2*N] = theta;

        for (int i = 0; i < N; i++)
        {
            state_down[i] = -1.0e19;
            state_up[i] = 1.0e19;

            state_down[N+i] = -1.0e19;
            state_up[N+i] = 1.0e19;

            state_down[2*N+i] = -1.0e19;
            state_up[2*N+i] = 1.0e19;

            state_down[3*N+i] = -0.5;
            state_up[3*N+i] = 1.0;

            state_down[4*N+i] = -0.3;
            state_up[4*N+i] = 0.3;
        }

        Dvector const_up(n_const), const_down(n_const);

        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = x;
        const_up[0] = x;

        const_down[N] = y;
        const_up[N] = y;

        const_down[2*N] = theta;
        const_up[2*N] = theta;

        for (int i = 0; i < N-1; i++)
        {
            const_down[3*N+i+1] = 0;
            const_up[3*N+i+1] = 1.0e19;
        }

        //const_down[5*N] = 0;
        //const_up[5*N] = fmax((x-20)*(x-20)+(y-20)*(y-20)-0.05, 1);

        std::string options;
        options += "Integer print_level  0\n";
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        options += "Numeric max_cpu_time 0.5\n";

        CppAD::ipopt::solve_result<Dvector> solution;

        CppAD::ipopt::solve<Dvector, FG_eval>(options, state, state_down, state_up, const_down, const_up, fg_eval, solution);

        double V = solution.x[3*N];
        double W = solution.x[4*N];
        std::vector<double> output = {V, W};

        return output;

    }
};


class FG_eval // the 麦伦
{
public:
    FG_eval(double x=20.0, double y=20.0)
    {
        x_ = x;
        y_ = y;
        r = 1;
        lambda = 0.1;
        alpha = 0.1;
    }
    void setXY(double x, double y)
    {
        x_ = x;
        y_ = y;
    }
    void setOXY(double x, double y, int i)
    {
        xo_[i] = x;
        yo_[i] = y;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < N-1; i++)
        {
            fg[0] += 0.1*(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2));
            //fg[0] += 10.0*(CppAD::pow(state[3*N+i+1]-state[3*N+i],2) + CppAD::pow(state[4*N+i+1]-state[4*N+i], 2) + CppAD::pow(state[5*N+i+1]-state[5*N+i], 2));
            //fg[0] += 5.0*(CppAD::pow(state[3*N+i]-0.5, 2) + CppAD::pow(state[4*N+i]-0.5, 2) + CppAD::pow(state[5*N+i]-0, 2));
        }
        for (int i = 0; i < N; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2) + CppAD::pow(state[5*N+i], 2));
        }
        fg[0] += 0.1*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2)); // The finnal point...

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];

        for (int i = 0; i < N-1; i++) // constraints
        {
            AD<double> x1 = state[i+1];
            AD<double> y1 = state[N+i+1];
            AD<double> theta1 = state[2*N+i+1];

            AD<double> x0 = state[i];
            AD<double> y0 = state[N+i];
            AD<double> theta0 = state[2*N+i];
            AD<double> vx0 = state[3*N+i];
            AD<double> vy0 = state[4*N+i];
            AD<double> w0 = state[5*N+i];

            fg[2 + i]    = x1 - (x0 + (cos(theta0) * vx0 - sin(theta0) * vy0) * dt); // note that use a trans point to overcome nonholonomic model
            fg[2 + N + i]    = y1 - (y0 + (sin(theta0) * vx0 + cos(theta0) * vy0) * dt);
            fg[2 + 2*N + i] = theta1 - (theta0 + w0 * dt);

            fg[2 + 3*N + i] = (CppAD::pow(x1-10, 2) + CppAD::pow(y1-10, 2) - 5*5) - (1-lambda)*(CppAD::pow(x0-10, 2) + CppAD::pow(y0-10, 2) - 5*5); // barrier cbf
            for (int j = 0; j < node_number-1; j++)
                fg[2+(3+j+1)*N+i] = (CppAD::pow(x1-xo_[j], 2) + CppAD::pow(y1-yo_[j], 2) - r*r) - (1-lambda)*(CppAD::pow(x0-xo_[j], 2) + CppAD::pow(y0-yo_[j], 2) - r*r); // agents cbf
            fg[2+(3+node_number)*N+i] = (CppAD::pow(vx0, 2) + CppAD::pow(vy0, 2));

        }
    }
private:
    double x_, y_, r, lambda, alpha;
    double xo_[20], yo_[20];
};

class GP_MPC : public rclcpp::Node
{
public:
    GP_MPC(): Node("GP_MPC")
    {
        node_name = this->get_name();

        for(int i = 0; i < node_number; i++)
        {
            std::string ads = node_name;
            ads[7] = '1' + i;
            odom_sub_list[i] = this->create_subscription<nav_msgs::msg::Odometry>(ads+"/odom/unfiltered", 10, [this, i, ads](const nav_msgs::msg::Odometry::SharedPtr msg){odom[i] = *msg; if(ads == this->node_name) sub_odom_callback(odom[i]);});
        }

        this->cmd_vel = this->create_publisher<geometry_msgs::msg::Twist>(node_name+"/cmd_vel", 10);
        this->timer = this->create_wall_timer(100ms, std::bind(&GP_MPC::timer_callback, this));

        times = 0;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    FG_eval fg_eval;

    std::string node_name;

    double x, y, theta, vx, vy, w, times;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x = odom1.pose.pose.position.x;
        y = odom1.pose.pose.position.y;
        theta = Orientation2Elur(odom1);
        vx = odom1.twist.twist.linear.x;
        vy = odom1.twist.twist.linear.y;
        w = odom1.twist.twist.angular.z;
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
        times++;
        std::vector<double> cmd = Solve();
        geometry_msgs::msg::Twist twist;
        twist.linear.x = cmd[0];
        twist.linear.y = cmd[1];
        twist.angular.z = cmd[2];
        cmd_vel->publish(twist);
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*(4+node_number);

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);

        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = x;
        state[N] = y;
        state[2*N] = theta;

        for (int i = 0; i < N; i++)
        {
            state_down[i] = -1.0e19;
            state_up[i] = 1.0e19;

            state_down[N+i] = -1.0e19;
            state_up[N+i] = 1.0e19;

            state_down[2*N+i] = -1.0e19;
            state_up[2*N+i] = 1.0e19;

            state_down[3*N+i] = -1.0;
            state_up[3*N+i] = 1.0;

            state_down[4*N+i] = -1.0;
            state_up[4*N+i] = 1.0;

            state_down[5*N+i] = -0.3;
            state_up[5*N+i] = 0.3;
        }

        Dvector const_up(n_const), const_down(n_const);

        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = x;
        const_up[0] = x;

        const_down[N] = y;
        const_up[N] = y;

        const_down[2*N] = theta;
        const_up[2*N] = theta;

        for (int i = 0; i < N-1; i++)
        {
            const_down[3*N+i+1] = 0;
            const_up[3*N+i+1] = 1.0e19;
            for (int j = 1; j < node_number; j++)
            {
                const_down[(3+j)*N+i+1] = 0.0;
                const_up[(3+j)*N+i+1] = 1.0e19;
            }
            const_down[(3+node_number)*N+i+1] = -0.4;
            const_up[(3+node_number)*N+i+1] = 0.4;
        }

        /*if (node_name[7] == '1')
            fg_eval.setOXY(odom[1].pose.pose.position.x, odom[1].pose.pose.position.y);
        else
            fg_eval.setOXY(odom[0].pose.pose.position.x, odom[0].pose.pose.position.y);*/

        for (int i = 0, j = 0; i < node_number; i++)
        {
            if (node_name[7] == '1'+i)
            {

                //fg_eval.setXY(10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                //RCLCPP_INFO(this->get_logger(), "i=%d, x=%lf, y=%lf\n", i, 10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                continue;
            }
            fg_eval.setOXY(odom[i].pose.pose.position.x, odom[i].pose.pose.position.y, j);
            j++;
        }

        //const_down[5*N] = 0;
        //const_up[5*N] = fmax((x-20)*(x-20)+(y-20)*(y-20)-0.05, 1);

        std::string options;
        options += "Integer print_level  0\n";
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        options += "Numeric max_cpu_time 0.5\n";

        CppAD::ipopt::solve_result<Dvector> solution;

        CppAD::ipopt::solve<Dvector, FG_eval>(options, state, state_down, state_up, const_down, const_up, fg_eval, solution);

        double VX = solution.x[3*N];
        double VY = solution.x[4*N];
        double W = solution.x[5*N];
        std::vector<double> output = {VX, VY, W};

        return output;

    }
};


class FG_eval // the project
{
public:
    FG_eval(double x=0.0, double y=0.0, double yaw=0.0)
    {
        x_ = x;
        y_ = y;
        yaw_ = yaw;
    }
    void setXY(double x, double y, double yaw)
    {
        x_ = x;
        y_ = y;
        yaw_ = yaw;
    }
    void setOXY(double x, double y, int i)
    {
        xo_[i] = x;
        yo_[i] = y;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < N-1; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2) + CppAD::pow(state[2*N+i]-yaw_, 2));
            //fg[0] += 10.0*(CppAD::pow(state[3*N+i+1]-state[3*N+i],2) + CppAD::pow(state[4*N+i+1]-state[4*N+i], 2) + CppAD::pow(state[5*N+i+1]-state[5*N+i], 2));
            //fg[0] += 5.0*(CppAD::pow(state[3*N+i]-0.5, 2) + CppAD::pow(state[4*N+i]-0.5, 2) + CppAD::pow(state[5*N+i]-0, 2));
            for (int j = 0; j < node_number-1; j++)
            {
                fg[0] +=20.0/(CppAD::pow(state[i]-xo_[j], 2) + CppAD::pow(state[N+i]-yo_[j], 2));
            }
        }
        for (int i = 0; i < N; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2) + CppAD::pow(state[5*N+i], 2));
        }
        fg[0] += 10.0*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2)); // The finnal point...

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];

        for (int i = 0; i < N-1; i++) // constraints
        {
            AD<double> x1 = state[i+1];
            AD<double> y1 = state[N+i+1];
            AD<double> theta1 = state[2*N+i+1];

            AD<double> x0 = state[i];
            AD<double> y0 = state[N+i];
            AD<double> theta0 = state[2*N+i];
            AD<double> vx0 = state[3*N+i];
            AD<double> vy0 = state[4*N+i];
            AD<double> w0 = state[5*N+i];

            fg[2 + i]    = x1 - (x0 + (cos(theta0) * vx0 - sin(theta0) * vy0) * dt); // note that use a trans point to overcome nonholonomic model
            fg[2 + N + i]    = y1 - (y0 + (sin(theta0) * vx0 + cos(theta0) * vy0) * dt);
            fg[2 + 2*N + i] = theta1 - (theta0 + w0 * dt);

            //for (int j = 0; j < node_number-1; j++)
                //fg[2+(3+j+1)*N+i] = (CppAD::pow(x1-xo_[j], 2) + CppAD::pow(y1-yo_[j], 2) - 2*2) - (1-0.1)*(CppAD::pow(x0-xo_[j], 2) + CppAD::pow(y0-yo_[j], 2) - 2*2); // agents cbf


        }
    }
private:
    double x_, y_, yaw_, r, lambda, alpha;
    double xo_[20], yo_[20];
};

class GP_MPC : public rclcpp::Node
{
public:
    GP_MPC(): Node("GP_MPC")
    {
        node_name = this->get_name();

        for(int i = 0; i < node_number; i++)
        {
            std::string ads = node_name;
            ads[7] = '1' + i;
            odom_sub_list[i] = this->create_subscription<nav_msgs::msg::Odometry>(ads+"/odom/unfiltered", 10, [this, i, ads](const nav_msgs::msg::Odometry::SharedPtr msg){odom[i] = *msg; if(ads == this->node_name) sub_odom_callback(odom[i]);});
        }

        this->cmd_vel = this->create_publisher<geometry_msgs::msg::Twist>(node_name+"/cmd_vel", 10);
        this->timer = this->create_wall_timer(100ms, std::bind(&GP_MPC::timer_callback, this));

        this->sub_formation = this->create_subscription<geometry_msgs::msg::Pose>("/hn/formation", 10, [this](const geometry_msgs::msg::Pose::SharedPtr msg){formation = *msg;});

        times = 0;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;
    rclcpp::Subscription<geometry_msgs::msg::Pose>::SharedPtr sub_formation;
    geometry_msgs::msg::Pose formation;

    FG_eval fg_eval;

    std::string node_name;

    double x, y, theta, vx, vy, w, times;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x = odom1.pose.pose.position.x;
        y = odom1.pose.pose.position.y;
        theta = Orientation2Elur(odom1);
        vx = odom1.twist.twist.linear.x;
        vy = odom1.twist.twist.linear.y;
        w = odom1.twist.twist.angular.z;
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
        times++;
        std::vector<double> cmd = Solve();
        geometry_msgs::msg::Twist twist;
        twist.linear.x = cmd[0];
        twist.linear.y = cmd[1];
        twist.angular.z = cmd[2];
        cmd_vel->publish(twist);
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*(4+node_number);

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);

        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = x;
        state[N] = y;
        state[2*N] = theta;

        for (int i = 0; i < N; i++)
        {
            state_down[i] = -1.0e19;
            state_up[i] = 1.0e19;

            state_down[N+i] = -1.0e19;
            state_up[N+i] = 1.0e19;

            state_down[2*N+i] = -1.0e19;
            state_up[2*N+i] = 1.0e19;

            state_down[3*N+i] = -0.25;
            state_up[3*N+i] = 0.25;

            state_down[4*N+i] = -0.0;
            state_up[4*N+i] = 0.0;

            state_down[5*N+i] = -0.15;
            state_up[5*N+i] = 0.15;
        }

        Dvector const_up(n_const), const_down(n_const);

        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = x;
        const_up[0] = x;

        const_down[N] = y;
        const_up[N] = y;

        const_down[2*N] = theta;
        const_up[2*N] = theta;

        if (node_name[7] == '1')
            fg_eval.setXY(0.1+x, y, theta);
        else if (node_name[7] == '2')
            fg_eval.setXY(odom[0].pose.pose.position.x-formation.orientation.x, odom[0].pose.pose.position.y-formation.orientation.y, Orientation2Elur(odom[0]));
        else
            fg_eval.setXY(odom[0].pose.pose.position.x-formation.orientation.z, odom[0].pose.pose.position.y-formation.orientation.w, Orientation2Elur(odom[0]));

        for (int i = 0; i < N-1; i++)
        {
            const_down[3*N+i+1] = 0;
            const_up[3*N+i+1] = 1.0e19;
            for (int j = 1; j < node_number; j++)
            {
                const_down[(3+j)*N+i+1] = 0.0;
                const_up[(3+j)*N+i+1] = 1.0e19;
            }
            const_down[(3+node_number)*N+i+1] = -0.4;
            const_up[(3+node_number)*N+i+1] = 0.4;
        }

        for (int i = 0, j = 0; i < node_number; i++)
        {
            if (node_name[7] == '1'+i)
            {

                //fg_eval.setXY(10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                //RCLCPP_INFO(this->get_logger(), "i=%d, x=%lf, y=%lf\n", i, 10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                continue;
            }
            fg_eval.setOXY(odom[i].pose.pose.position.x, odom[i].pose.pose.position.y, j);
            j++;
        }

        std::string options;
        options += "Integer print_level  0\n";
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        options += "Numeric max_cpu_time 0.5\n";

        CppAD::ipopt::solve_result<Dvector> solution;

        CppAD::ipopt::solve<Dvector, FG_eval>(options, state, state_down, state_up, const_down, const_up, fg_eval, solution);

        double VX = solution.x[3*N];
        double VY = solution.x[4*N];
        double W = solution.x[5*N];
        std::vector<double> output = {VX, VY, W};

        return output;

    }
};
#endif

class FG_eval // the SE(3)
{
public:
    FG_eval(double x=20.0, double y=20.0)
    {
        x_ = x;
        y_ = y;
        r = 1;
        lambda = 0.1;
        alpha = 0.1;
    }
    void setXY(double x, double y)
    {
        x_ = x;
        y_ = y;
    }
    void setOXY(double x, double y, int i)
    {
        xo_[i] = x;
        yo_[i] = y;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < N-1; i++)
        {
            AD<double> d1 = CppAD::sqrt(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2));
            fg[0] += 2.0*(CppAD::sqrt((2+d1*d1+d1*CppAD::sqrt(d1*d1+4))/2));
        }
        for (int i = 0; i < N; i++)
        {
            fg[0] += 0.5*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2) + CppAD::pow(state[5*N+i], 2));
        }
        AD<double> d = CppAD::sqrt(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2));
        fg[0] += 2.0*(CppAD::sqrt((2+d*d+d*CppAD::sqrt(d*d+4))/2));
        //fg[0] += 0.1*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2)); // The finnal point...

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];

        for (int i = 0; i < N-1; i++) // constraints
        {
            AD<double> x1 = state[i+1];
            AD<double> y1 = state[N+i+1];
            AD<double> theta1 = state[2*N+i+1];

            AD<double> x0 = state[i];
            AD<double> y0 = state[N+i];
            AD<double> theta0 = state[2*N+i];
            AD<double> vx0 = state[3*N+i];
            AD<double> vy0 = state[4*N+i];
            AD<double> w0 = state[5*N+i];

            AD<double> d1 = CppAD::sqrt(CppAD::pow(x1-10, 2) + CppAD::pow(y1-10, 2));
            AD<double> d0 = CppAD::sqrt(CppAD::pow(x0-10, 2) + CppAD::pow(y0-10, 2));

            fg[2 + i]    = x1 - (x0 + (cos(theta0) * vx0 - sin(theta0) * vy0) * dt); // note that use a trans point to overcome nonholonomic model
            fg[2 + N + i]    = y1 - (y0 + (sin(theta0) * vx0 + cos(theta0) * vy0) * dt);
            fg[2 + 2*N + i] = theta1 - (theta0 + w0 * dt);

            fg[2 + 3*N + i] = (CppAD::sqrt((2+d1*d1+d1*CppAD::sqrt(d1*d1+4))/2)-5.2) - (1-lambda)*(CppAD::sqrt((2+d0*d0+d0*CppAD::sqrt(d0*d0+4))/2)-5.2); // barrier cbf

            fg[2+(3+node_number)*N+i] = (CppAD::pow(vx0, 2) + CppAD::pow(vy0, 2));

        }
    }
private:
    double x_, y_, r, lambda, alpha;
    double xo_[20], yo_[20];
};

class GP_MPC : public rclcpp::Node
{
public:
    GP_MPC(): Node("GP_MPC")
    {
        node_name = this->get_name();

        for(int i = 0; i < node_number; i++)
        {
            std::string ads = node_name;
            ads[7] = '1' + i;
            odom_sub_list[i] = this->create_subscription<nav_msgs::msg::Odometry>(ads+"/odom/unfiltered", 10, [this, i, ads](const nav_msgs::msg::Odometry::SharedPtr msg){odom[i] = *msg; if(ads == this->node_name) sub_odom_callback(odom[i]);});
        }

        this->cmd_vel = this->create_publisher<geometry_msgs::msg::Twist>(node_name+"/cmd_vel", 10);
        this->timer = this->create_wall_timer(100ms, std::bind(&GP_MPC::timer_callback, this));

        times = 0;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    FG_eval fg_eval;

    std::string node_name;

    double x, y, theta, vx, vy, w, times;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x = odom1.pose.pose.position.x;
        y = odom1.pose.pose.position.y;
        theta = Orientation2Elur(odom1);
        vx = odom1.twist.twist.linear.x;
        vy = odom1.twist.twist.linear.y;
        w = odom1.twist.twist.angular.z;
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
        times++;
        std::vector<double> cmd = Solve();
        geometry_msgs::msg::Twist twist;
        twist.linear.x = cmd[0];
        twist.linear.y = cmd[1];
        twist.angular.z = cmd[2];
        cmd_vel->publish(twist);
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*(4+node_number);

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);

        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = x;
        state[N] = y;
        state[2*N] = theta;

        for (int i = 0; i < N; i++)
        {
            state_down[i] = -1.0e19;
            state_up[i] = 1.0e19;

            state_down[N+i] = -1.0e19;
            state_up[N+i] = 1.0e19;

            state_down[2*N+i] = -1.0e19;
            state_up[2*N+i] = 1.0e19;

            state_down[3*N+i] = -1.0;
            state_up[3*N+i] = 1.0;

            state_down[4*N+i] = -1.0;
            state_up[4*N+i] = 1.0;

            state_down[5*N+i] = -0.3;
            state_up[5*N+i] = 0.3;
        }

        Dvector const_up(n_const), const_down(n_const);

        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = x;
        const_up[0] = x;

        const_down[N] = y;
        const_up[N] = y;

        const_down[2*N] = theta;
        const_up[2*N] = theta;

        for (int i = 0; i < N-1; i++)
        {
            const_down[3*N+i+1] = 0;
            const_up[3*N+i+1] = 1.0e19;
            for (int j = 1; j < node_number; j++)
            {
                const_down[(3+j)*N+i+1] = 0.0;
                const_up[(3+j)*N+i+1] = 1.0e19;
            }
            const_down[(3+node_number)*N+i+1] = -0.4;
            const_up[(3+node_number)*N+i+1] = 0.4;
        }

        /*if (node_name[7] == '1')
            fg_eval.setOXY(odom[1].pose.pose.position.x, odom[1].pose.pose.position.y);
        else
            fg_eval.setOXY(odom[0].pose.pose.position.x, odom[0].pose.pose.position.y);*/

        for (int i = 0, j = 0; i < node_number; i++)
        {
            if (node_name[7] == '1'+i)
            {

                //fg_eval.setXY(10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                //RCLCPP_INFO(this->get_logger(), "i=%d, x=%lf, y=%lf\n", i, 10-10*cos(2.0*double(i)*3.1415926/double(node_number)), 10-10*sin(2.0*double(i)*3.1415926/double(node_number)));
                continue;
            }
            fg_eval.setOXY(odom[i].pose.pose.position.x, odom[i].pose.pose.position.y, j);
            j++;
        }

        //const_down[5*N] = 0;
        //const_up[5*N] = fmax((x-20)*(x-20)+(y-20)*(y-20)-0.05, 1);

        std::string options;
        options += "Integer print_level  0\n";
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        options += "Numeric max_cpu_time 0.5\n";

        CppAD::ipopt::solve_result<Dvector> solution;

        CppAD::ipopt::solve<Dvector, FG_eval>(options, state, state_down, state_up, const_down, const_up, fg_eval, solution);

        double VX = solution.x[3*N];
        double VY = solution.x[4*N];
        double W = solution.x[5*N];
        std::vector<double> output = {VX, VY, W};

        return output;

    }
};
