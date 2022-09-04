#include "Robust_MPC/Robust_MPC.hpp"

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<Robust_MPC>());
    rclcpp::shutdown();
    return 0;
}
