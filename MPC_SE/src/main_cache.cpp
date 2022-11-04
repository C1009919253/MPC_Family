#include "MPC_SE/MPC_SE.hpp"

int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<MPC_SE>());
    rclcpp::shutdown();
    return 0;
}
