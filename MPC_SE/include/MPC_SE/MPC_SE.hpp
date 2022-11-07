#pragma once

#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <functional>
#include <memory>
#include <string>
#include <cmath>
#include <list>
#include <sstream>
#include <fstream>
#include <vector>

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


#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#define NX          ACADO_NX	/* number of differential states */
#define NXA         ACADO_NXA	/* number of alg. states */
#define NU          ACADO_NU	/* number of control inputs */
#define N          	ACADO_N		/* number of control intervals */
#define NY			ACADO_NY	/* number of references, nodes 0..N - 1 */
#define NYN			ACADO_NYN
#define NUM_STEPS   100			/* number of simulation steps */
#define VERBOSE     1			/* show iterations: 1, silent: 0 */

//static size_t N = 6;
static double dt = 0.1;

static int node_number=2;

using namespace std::chrono_literals;

ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

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

        std::ifstream ifsx, ifsy, ifst, ifsv, ifsw;
        ifsx.open("x-t.txt");
        ifsy.open("y-t.txt");
        ifst.open("yaw-t.txt");
        ifsv.open("v-t.txt");
        ifsw.open("w-t.txt");
        std::string str;
        while (getline(ifsx, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a-0.15;
                i++;
            }
            xd.push_back(data);
        }
        while (getline(ifsy, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a-0.25;
                i++;
            }
            yd.push_back(data);
        }
        while (getline(ifst, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a;
                i++;
            }
            td.push_back(data);
        }
        while (getline(ifsv, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                //data[i] = a;
                i++;
            }
            vd.push_back(data);
        }
        while (getline(ifsw, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                //data[i] = a;
                i++;
            }
            wd.push_back(data);
        }

        for (int i = 0; i < xd.size(); i++)
        {
            geometry_msgs::msg::PoseStamped point;
            point.pose.position.x = xd[i][0]-0.15;
            point.pose.position.y = yd[i][0]-0.25;
            //desire_path.poses.push_back(point);
        }

        times = 0;

        memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
        memset(&acadoVariables, 0, sizeof( acadoVariables ));

        acado_initializeSolver();

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    std::string node_name;
    char robot_num;

    double x0, y0, theta0;

    std::vector<double *> xd, yd, td, vd, wd;

    int times;




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
        //
        // Prepare a consistent initial guess
        //
        int i;

        for (i = 0; i < N + 1; ++i)
        {
            acadoVariables.x[i * NX + 0] = x0;
            acadoVariables.x[i * NX + 1] = y0;
            acadoVariables.x[i * NX + 2] = theta0;
        }

        //
        // Prepare references
        //

        for (i = 0; i < N; ++i)
        {
            acadoVariables.y[i * NY + 0] = 0.0; // x
            acadoVariables.y[i * NY + 1] = 0.0; // y
            acadoVariables.y[i * NY + 2] = 0; // w
            acadoVariables.y[i * NY + 3] = 1.0;
            acadoVariables.y[i * NY + 4] = 0;
        }

        acadoVariables.yN[ 0 ] = 0.0; // x
        acadoVariables.yN[ 1 ] = 0.0; // y
        acadoVariables.yN[ 2 ] = 0; // w

        //
        // Current state feedback
        //
        for (i = 0; i < NX; ++i)
            acadoVariables.x0[ i ] = acadoVariables.x[ i ];

        acado_preparationStep();

        int status;

        status = acado_feedbackStep();

        if (status)
        {
            std::cout<<"error!!!"<<std::endl;
        }

        geometry_msgs::msg::Twist twist;
        twist.linear.x = acadoVariables.u[0];
        twist.angular.z = acadoVariables.u[1];
        cmd_vel->publish(twist);

    }
    /*std::vector<double> Solve()
    {

    }*/
};

#if 0
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

        std::ifstream ifsx, ifsy, ifst, ifsv, ifsw;
        ifsx.open("x-t.txt");
        ifsy.open("y-t.txt");
        ifst.open("yaw-t.txt");
        ifsv.open("v-t.txt");
        ifsw.open("w-t.txt");
        std::string str;
        while (getline(ifsx, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a-0.15;
                i++;
            }
            xd.push_back(data);
        }
        while (getline(ifsy, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a-0.25;
                i++;
            }
            yd.push_back(data);
        }
        while (getline(ifst, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                data[i] = a;
                i++;
            }
            td.push_back(data);
        }
        while (getline(ifsv, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                //data[i] = a;
                i++;
            }
            vd.push_back(data);
        }
        while (getline(ifsw, str))
        {
            double * data = new double[2]{0.0, 0.0};
            std::stringstream number(str);
            std::string out;
            int i = 0;
            while (number >> out)
            {
                double a = std::stod(out);
                //data[i] = a;
                i++;
            }
            wd.push_back(data);
        }

        for (int i = 0; i < xd.size(); i++)
        {
            geometry_msgs::msg::PoseStamped point;
            point.pose.position.x = xd[i][0]-0.15;
            point.pose.position.y = yd[i][0]-0.25;
            //desire_path.poses.push_back(point);
        }

        times = 0;

    }
private:
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_list[20];
    rclcpp::Publisher<geometry_msgs::msg::Twist>::SharedPtr cmd_vel;
    nav_msgs::msg::Odometry odom[20];
    rclcpp::TimerBase::SharedPtr timer;

    std::string node_name;
    char robot_num;

    double x0, y0, theta0;

    std::vector<double *> xd, yd, td, vd, wd;

    int times;

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

        DifferentialState xe, ye, thetae;
        Control v, w;
        DiscretizedDifferentialEquation f(dt);
        Function h;

        DMatrix Q;
        DVector r;

        xe.clearStaticCounters(); // clean all buffer...
        ye.clearStaticCounters();
        thetae.clearStaticCounters();
        v.clearStaticCounters();
        w.clearStaticCounters();
        f.clearBuffer();
        h.clearBuffer();

        f << next(xe) == cos(wd[times][0]*dt)*xe + sin(wd[times][0]*dt)*ye + cos(thetae-wd[times][0])*v*dt - vd[times][0]*dt; // BUG!!!
        f << next(ye) == -sin(wd[times][0]*dt)*xe + cos(wd[times][0]*dt)*ye + sin(thetae-wd[times][0])*v*dt;
        f << next(thetae) == thetae + w*dt - wd[times][0]*dt;

        h << xe;
        h << ye;
        h << thetae;

        Q.resize(3, 3);
        Q.setIdentity();

        r.resize(2);
        r.setAll(0.0);

        OCP ocp(0.0, 1.0, 10);

        ocp.minimizeLSQ(Q, h, r);
        //ocp.minimizeLagrangeTerm(xe*xe+ye*ye);
        ocp.subjectTo(f);
        ocp.subjectTo(AT_START, xe == cos(td[times][0])*x0+sin(td[times][0])*y0-xd[times][0]);
        ocp.subjectTo(AT_START, ye == -sin(td[times][0])*x0+cos(td[times][0])*y0-yd[times][0]);
        ocp.subjectTo(AT_START, thetae == theta0-td[times][0]);
        ocp.subjectTo(-1.0<=v<=1.0);
        ocp.subjectTo(-1.0<=w<=1.0);
        ocp.subjectTo(-1.0<=dot(v)<=1.0);

        RealTimeAlgorithm alg(ocp);
OCPexport mpc( ocp );
        Controller controller(alg);
        DVector asd(4);
        asd(0) = cos(td[times][0])*x0+sin(td[times][0])*y0-xd[times][0];
        asd(1) = -sin(td[times][0])*x0+cos(td[times][0])*y0-yd[times][0];
        asd(2) = theta0-td[times][0];
        controller.init(0.0, asd);
        controller.step(0.0, asd);

        /*OptimizationAlgorithm algorithm(ocp);
        algorithm.solve();

        VariablesGrid controls;

        algorithm.getControls(controls);

        DMatrix U = controls.getMatrix(0);

        geometry_msgs::msg::Twist twist;
        twist.linear.x = U(0);
        //twist.linear.y = U(1);
        twist.angular.z = U(1);

        cmd_vel->publish(twist);
        times++;*/
    }
    /*std::vector<double> Solve()
    {

    }*/
};

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
#endif
