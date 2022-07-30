#include "GP_MPC/GP_MPC.hpp"

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<GP_MPC>());
    rclcpp::shutdown();
    return 0;
}
