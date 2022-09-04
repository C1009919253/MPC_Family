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

#include "Robust_MPC/gnuplot_i.hpp"
#include "Robust_MPC/DLQR.hpp"

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

using namespace std::chrono_literals;
using CppAD::AD;

class Data_Driven_Compute
{
public:
    Data_Driven_Compute(Eigen::MatrixXd H, Eigen::MatrixXd w)
    {
        H_ = H;
        w_ = w;
    }
    typedef CPPAD_TESTVECTOR( AD<double> ) ADvector;
    void operator ()(ADvector& fg, const ADvector& state)
    {
        fg[0] = 0;
        for (int i = 0; i < H_.rows(); i++)
        {
            AD<double> miws = 0;
            for (int j = 0; j < H_.rows(); j++)
            {
                miws += H_(i, j) * state[j];
            }
            fg[0] += CppAD::pow(miws, 2);
        }
    }
private:
    Eigen::MatrixXd H_, w_;
};

class Data_Driven_Model
{
public:
    Data_Driven_Model(int N, int L, int nu, Eigen::MatrixXd Hu, Eigen::MatrixXd w) // N:激励数量，L：预测视野,nu:输入阶数, Hu:输入汉科尔矩阵，w：预测输入序列
    {
        N_ = N;
        L_ = L;
        nu_ = nu;
        H_ = Hu;
        w_ = w;
    }
    Data_Driven_Model(int N, int L, int nu, int nx, Eigen::MatrixXd Hu, Eigen::MatrixXd Hx, Eigen::MatrixXd w)
    {
        N_ = N;
        L_ = L;
        nu_ = nu;
        nx_ = nx;
        H_ = Eigen::MatrixXd::Zero(Hu.rows()+Hx.rows(), Hu.cols());
        H_.block(0, 0, Hu.rows(), Hu.cols()) = Hu;
        H_.block(Hu.rows(), 0, Hx.rows(), Hx.cols()) = Hx;
        w_ = w;
    }
    std::vector<double> Solve()
    {
        typedef CPPAD_TESTVECTOR(double) Dvector;

        int n_state = N_ - L_ - nu_ + 1;
        int n_const = 0;

        Dvector state(n_state);
        Dvector state_up(n_state), state_down(n_state);
#pragma omp parallel for
        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = 0;
            state_up[i] = 0;
            state[i] = 0;
        }
#pragma omp parallel for
        for (int i = 0; i < n_state; i++)
        {
            state_down[i] = -1.0e19;
            state_up[i] = 1.0e19;
        }
        Dvector const_up(n_const), const_down(n_const);
#pragma omp parallel for
        for (int i = 0; i < n_const; i++)
        {
            const_down[i] = 0;
            const_up[i] = 0;
        }

        Data_Driven_Compute fg_eval(H_, w_);

        std::string options;
        options += "Integer print_level  0\n";
        options += "Sparse  true        forward\n";
        options += "Sparse  true        reverse\n";
        options += "Numeric max_cpu_time 0.5\n";

        CppAD::ipopt::solve_result<Dvector> solution;

        CppAD::ipopt::solve<Dvector, Data_Driven_Compute>(options, state, state_down, state_up, const_down, const_up, fg_eval, solution);

        std::vector<double> output;

        for (int i = 0; i < n_state; i++)
            output.push_back(solution.x[i]);

        return output;

    }
    int N_, L_, nu_, nx_;
    Eigen::MatrixXd H_, w_;
private:
};
