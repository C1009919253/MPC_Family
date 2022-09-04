#pragma once

#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <cmath>
#include <list>

#include "gurobi_c++.h"
#include <eigen3/Eigen/Dense>

#include "Robust_MPC/gnuplot_i.hpp"
#include "Robust_MPC/DLQR.hpp"
#include "Robust_MPC/Data_Driven_Model.hpp"

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

#include <gurobi_c++.h>

#include <omp.h> // multi-core

static size_t N = 6;
static double dt = 0.1;

static int node_number=2;

using namespace std::chrono_literals;
using CppAD::AD;

typedef struct state_input
{
    double x, y, theta, vx, vy, w;
}state_input;

class FG_eval // the SE(2)
{
public:
    FG_eval(double x=20.0, double y=0.0)
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
    void sets_i(std::list<state_input> si) // 步骤：首先使用s_i得到阿尔法，通过state构建汉科尔矩阵乘阿尔法得到模型输出。
    {
        s_i = si;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < N-1; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2) + CppAD::pow(state[2*N+i], 2));
            //fg[0] += 10.0*(CppAD::pow(state[3*N+i]-state[3*N+i+1],2) + CppAD::pow(state[4*N+i]-state[4*N+i+1], 2) + CppAD::pow(state[5*N+i]-state[5*N+i+1], 2));
        }
        for (int i = 0; i < N; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2) + CppAD::pow(state[5*N+i], 2));
        }
        fg[0] += 1.0*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2) + CppAD::pow(state[3*N-1], 2));

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];
#pragma omp parallel for
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


            fg[2 + 3*N + i] = (CppAD::pow(x1-10.0, 2) + CppAD::pow(y1-1.0, 2)-25.0) - (1-lambda)*(CppAD::pow(x0-10.0, 2) + CppAD::pow(y0-1.0, 2)-25.0); // barrier cbf

            fg[2+(3+node_number)*N+i] = (CppAD::pow(vx0, 2) + CppAD::pow(vy0, 2));

        }
    }
private:
    double x_, y_, r, lambda, alpha;
    double xo_[20], yo_[20];
    std::list<state_input> s_i;
};

class Robust_MPC : public rclcpp::Node
{
public:
    Robust_MPC(): Node("Robust_MPC")
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
        this->timer = this->create_wall_timer(100ms, std::bind(&Robust_MPC::timer_callback, this));

        times = 0;

        first = 0;

        rand_noisy = 0.0;

        L_u = N - 1;

        N_u = 50;

        n_u = 3;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    FG_eval fg_eval;

    std::string node_name;

    double x, y, theta, vx, vy, w;
    double zx, zy, ztheta, zvx, zvy, zw;

    double rand_noisy;

    int first, times;

    char robot_num;

    Eigen::MatrixXd Hu, Hx, w;

    int n_u, N_u, L_u;

    std::list<state_input> s_i;

    //DLQR *dlqr;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x = odom1.pose.pose.position.x;
        y = odom1.pose.pose.position.y;
        theta = Orientation2Elur(odom1);
        vx = odom1.twist.twist.linear.x;
        vy = odom1.twist.twist.linear.y;
        w = odom1.twist.twist.angular.z;

        if (first == 0)
        {
            zx = odom1.pose.pose.position.x;
            zy = odom1.pose.pose.position.y;
            ztheta = Orientation2Elur(odom1);
            zvx = odom1.twist.twist.linear.x;
            zvy = odom1.twist.twist.linear.y;
            zw = odom1.twist.twist.angular.z;
            first++;
        }

        state_input kkkj;
        kkkj.x = x;
        kkkj.y = y;
        kkkj.theta = theta;
        kkkj.vx = vx;
        kkkj.vy = vy;
        kkkj.w = w;

        if (s_i.size() == 50)
            s_i.pop_front();
        s_i.push_back(kkkj);

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

        // to estimate system
        zx = zx + (cos(ztheta) * cmd[0] - sin(ztheta) * cmd[1]) * dt;
        zy = zy + (sin(ztheta) * cmd[0] + cos(ztheta) * cmd[1]) * dt;
        ztheta = ztheta + cmd[2] * dt;
        zvx = cmd[0];
        zvy = cmd[1];
        zw = cmd[2];

        /*Eigen::MatrixXd x_cur;
        x_cur.resize(3, 1);
        x_cur << x,
                 y,
                 theta;
        Eigen::MatrixXd A;
        A.resize(3, 3);
        A << 1.0, 0, 0,
             0, 1.0, 0,
             0, 0, 1.0;
        Eigen::MatrixXd B;
        B.resize(3, 3);
        B << cos(theta)*dt, -sin(theta)*dt, 0,
             sin(theta)*dt,  cos(theta)*dt, 0,
                   0,            0,    1.0*dt;
        DLQR *dlqr = new DLQR(A, B, 3, 3, dt, DISCRETE);

        dlqr->DLQRInit();
        dlqr->DLQRRun();

        auto k = dlqr->K;
        Eigen::MatrixXd x_now;
        x_now.resize(3,1);
        x_now<<x-zx,y-zy,theta-ztheta;
        Eigen::MatrixXd u = -k*x_now;*/


        Eigen::MatrixXd x_now, k;
        x_now.resize(3,1);
        k.resize(3, 3);
        double que = 10.0;
        k <<  cos(theta)*que, sin(theta)*que, 0,
             -sin(theta)*que, cos(theta)*que, 0,
                    0,            0,    1.0*que;
        x_now<<x-zx,y-zy,theta-ztheta;
        Eigen::MatrixXd u = -k*x_now;

        if (times % 20 == 0)
        {
            //rand_noisy = drand48() * drand48() * 2.0 - 1.0;
        }
        if (true)
        {
            rand_noisy = rand_noisy*drand48()+drand48()*0.5-0.25;
        }
        rand_noisy=0;
        geometry_msgs::msg::Twist twist;
        if (robot_num == '2')
        {
            twist.linear.x = 1.0*cmd[0]+1.0*u(0);
            twist.linear.y = 1.0*cmd[1]+1.0*u(1)+rand_noisy;
            twist.angular.z = 1.0*cmd[2]+1.0*u(2);
        }
        else
        {
            twist.linear.x = 1.0*cmd[0]+0.0*u(0);
            twist.linear.y = 1.0*cmd[1]+0.0*u(1)+rand_noisy;
            twist.angular.z = 1.0*cmd[2]+0.0*u(2);
        }
        cmd_vel->publish(twist);
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*(4+node_number);

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);
#pragma omp parallel for
        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = zx;
        state[N] = zy;
        state[2*N] = ztheta;
#pragma omp parallel for
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

            state_down[5*N+i] = -0.1;
            state_up[5*N+i] = 0.1;
        }

        Dvector const_up(n_const), const_down(n_const);
#pragma omp parallel for
        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = zx;
        const_up[0] = zx;

        const_down[N] = zy;
        const_up[N] = zy;

        const_down[2*N] = ztheta;
        const_up[2*N] = ztheta;
#pragma omp parallel for
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
#pragma omp parallel for
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
        if (robot_num == '2')
            fg_eval.setXY(20, -20);

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

#if 0
class FG_eval // the SE(2)
{
public:
    FG_eval(double x=20.0, double y=0.0)
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
            fg[0] += 1.0*(CppAD::pow(state[i]-x_, 2) + CppAD::pow(state[N+i]-y_, 2) + CppAD::pow(state[2*N+i], 2));
            //fg[0] += 10.0*(CppAD::pow(state[3*N+i]-state[3*N+i+1],2) + CppAD::pow(state[4*N+i]-state[4*N+i+1], 2) + CppAD::pow(state[5*N+i]-state[5*N+i+1], 2));
        }
        for (int i = 0; i < N; i++)
        {
            fg[0] += 1.0*(CppAD::pow(state[3*N+i],2) + CppAD::pow(state[4*N+i], 2) + CppAD::pow(state[5*N+i], 2));
        }
        fg[0] += 1.0*(CppAD::pow(state[N-1]-x_, 2) + CppAD::pow(state[2*N-1]-y_, 2) + CppAD::pow(state[3*N-1], 2));

        fg[1] = state[0];
        fg[N+1] = state[N];
        fg[2*N+1] = state[2*N];
#pragma omp parallel for
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


            fg[2 + 3*N + i] = (CppAD::pow(x1-10.0, 2) + CppAD::pow(y1-1.0, 2)-25.0) - (1-lambda)*(CppAD::pow(x0-10.0, 2) + CppAD::pow(y0-1.0, 2)-25.0); // barrier cbf

            fg[2+(3+node_number)*N+i] = (CppAD::pow(vx0, 2) + CppAD::pow(vy0, 2));

        }
    }
private:
    double x_, y_, r, lambda, alpha;
    double xo_[20], yo_[20];
};

class Robust_MPC : public rclcpp::Node
{
public:
    Robust_MPC(): Node("Robust_MPC")
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
        this->timer = this->create_wall_timer(100ms, std::bind(&Robust_MPC::timer_callback, this));

        times = 0;

        first = 0;

        rand_noisy = 0.0;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    FG_eval fg_eval;

    std::string node_name;

    double x, y, theta, vx, vy, w;
    double zx, zy, ztheta, zvx, zvy, zw;

    double rand_noisy;

    int first, times;

    char robot_num;

    //DLQR *dlqr;

    void sub_odom_callback(nav_msgs::msg::Odometry odom1)
    {
        x = odom1.pose.pose.position.x;
        y = odom1.pose.pose.position.y;
        theta = Orientation2Elur(odom1);
        vx = odom1.twist.twist.linear.x;
        vy = odom1.twist.twist.linear.y;
        w = odom1.twist.twist.angular.z;

        if (first == 0)
        {
            zx = odom1.pose.pose.position.x;
            zy = odom1.pose.pose.position.y;
            ztheta = Orientation2Elur(odom1);
            zvx = odom1.twist.twist.linear.x;
            zvy = odom1.twist.twist.linear.y;
            zw = odom1.twist.twist.angular.z;
            first++;
        }
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

        // to estimate system
        zx = zx + (cos(ztheta) * cmd[0] - sin(ztheta) * cmd[1]) * dt;
        zy = zy + (sin(ztheta) * cmd[0] + cos(ztheta) * cmd[1]) * dt;
        ztheta = ztheta + cmd[2] * dt;
        zvx = cmd[0];
        zvy = cmd[1];
        zw = cmd[2];

        /*Eigen::MatrixXd x_cur;
        x_cur.resize(3, 1);
        x_cur << x,
                 y,
                 theta;
        Eigen::MatrixXd A;
        A.resize(3, 3);
        A << 1.0, 0, 0,
             0, 1.0, 0,
             0, 0, 1.0;
        Eigen::MatrixXd B;
        B.resize(3, 3);
        B << cos(theta)*dt, -sin(theta)*dt, 0,
             sin(theta)*dt,  cos(theta)*dt, 0,
                   0,            0,    1.0*dt;
        DLQR *dlqr = new DLQR(A, B, 3, 3, dt, DISCRETE);

        dlqr->DLQRInit();
        dlqr->DLQRRun();

        auto k = dlqr->K;
        Eigen::MatrixXd x_now;
        x_now.resize(3,1);
        x_now<<x-zx,y-zy,theta-ztheta;
        Eigen::MatrixXd u = -k*x_now;*/


        Eigen::MatrixXd x_now, k;
        x_now.resize(3,1);
        k.resize(3, 3);
        double que = 10.0;
        k <<  cos(theta)*que, sin(theta)*que, 0,
             -sin(theta)*que, cos(theta)*que, 0,
                    0,            0,    1.0*que;
        x_now<<x-zx,y-zy,theta-ztheta;
        Eigen::MatrixXd u = -k*x_now;

        if (times % 20 == 0)
        {
            //rand_noisy = drand48() * drand48() * 2.0 - 1.0;
        }
        if (true)
        {
            rand_noisy = rand_noisy*drand48()+drand48()*0.5-0.25;
        }
        //rand_noisy=0;
        geometry_msgs::msg::Twist twist;
        if (robot_num == '1')
        {
            twist.linear.x = 1.0*cmd[0]+1.0*u(0);
            twist.linear.y = 1.0*cmd[1]+1.0*u(1)+rand_noisy;
            twist.angular.z = 1.0*cmd[2]+1.0*u(2);
        }
        else
        {
            twist.linear.x = 1.0*cmd[0]+0.0*u(0);
            twist.linear.y = 1.0*cmd[1]+0.0*u(1)+rand_noisy;
            twist.angular.z = 1.0*cmd[2]+0.0*u(2);
        }
        cmd_vel->publish(twist);
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N*6;
        int n_const = N*(4+node_number);

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);
#pragma omp parallel for
        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }

        state[0] = zx;
        state[N] = zy;
        state[2*N] = ztheta;
#pragma omp parallel for
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

            state_down[5*N+i] = -0.1;
            state_up[5*N+i] = 0.1;
        }

        Dvector const_up(n_const), const_down(n_const);
#pragma omp parallel for
        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        const_down[0] = zx;
        const_up[0] = zx;

        const_down[N] = zy;
        const_up[N] = zy;

        const_down[2*N] = ztheta;
        const_up[2*N] = ztheta;
#pragma omp parallel for
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
#pragma omp parallel for
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
        if (robot_num == '2')
            fg_eval.setXY(20, -20);

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
#endif
