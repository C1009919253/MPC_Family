cmake_minimum_required(VERSION 3.5)
project(Robust_MPC)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(LINK_DIR /usr/local/lib
/usr/lib
/usr/lib/x86_64-linux-gnu/
/home/wsf/Git_Src/bayesopt/lib/
/home/wsf/Git_Src/bayesopt/
/home/wsf/Git_Src/bayesopt/include/
)

link_directories(${LINK_DIR})

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(tf2_msgs REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(blasfeo REQUIRED)
find_package(hpipm REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(ct_core REQUIRED)
find_package(ct_rbd REQUIRED)
find_package(ct_optcon REQUIRED)
find_package(ct_models REQUIRED)
find_package(manif REQUIRED)
#find_package(boost REQUIRED)
#find_package(bayesopt REQUIRED)

include_directories(include)

include(./FindGUROBI.cmake)
include_directories(${GUROBI_INCLUDE_DIRS}
/home/wsf/Git_Src/bayesopt/include/
/home/wsf/Git_Src/bayesopt/src/
)

add_executable(Robust_MPC
  src/main.cpp
  include/Robust_MPC/DLQR.cpp
)

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${manif_INCLUDE_DIRS})

target_link_libraries(Robust_MPC
  ${GUROBI_LIBRARIES}
)

ament_target_dependencies(Robust_MPC
  "rclcpp"
  "std_msgs"
  "geometry_msgs"
  "tf2_msgs"
  "tf2"
  "tf2_geometry_msgs"
  "nav_msgs"
  "Eigen3"
  "blasfeo"
  "hpipm"
  "ct_core"
  "ct_rbd"
  "ct_optcon"
  "ct_models"
  "manif"
)

install(TARGETS Robust_MPC
 DESTINATION lib/${PROJECT_NAME})
 
 install(DIRECTORY
launch
worlds
DESTINATION share/${PROJECT_NAME}/
)

include_directories(
  ${GUROBI_INCLUDE_DIRS}
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # uncomment the line when a copyright and license is not present in all source files
  #set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # uncomment the line when this package is not in a git repo
  #set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
