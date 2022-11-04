cmake_minimum_required(VERSION 3.5)
project(MPC_SE)

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

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(tf2_msgs REQUIRED)
find_package(tf2 REQUIRED)
find_package(tf2_geometry_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories( include/)

add_library(nmpc_lib
        solver/qpoases/SRC/Bounds.cpp
        solver/qpoases/SRC/Constraints.cpp
        solver/qpoases/SRC/CyclingManager.cpp
        solver/qpoases/SRC/Indexlist.cpp
        solver/qpoases/SRC/MessageHandling.cpp
        solver/qpoases/SRC/QProblem.cpp
        solver/qpoases/SRC/QProblemB.cpp
        solver/qpoases/SRC/SubjectTo.cpp
        solver/qpoases/SRC/Utils.cpp
        solver/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp
        solver/OCPexport/acado_qpoases_interface.cpp
        solver/OCPexport/acado_integrator.c
        solver/OCPexport/acado_solver.c
        solver/OCPexport/acado_auxiliary_functions.c
)

target_include_directories(nmpc_lib PUBLIC
        solver/
        solver/OCPexport
        solver/qpoases
        solver/qpoases/INCLUDE
        solver/qpoases/SRC
        )

add_executable(MPC_SE src/main.cpp)

target_link_libraries(MPC_SE nmpc_lib)

ament_target_dependencies(MPC_SE
  rclcpp
  std_msgs
  geometry_msgs
  tf2_msgs
  tf2
  tf2_geometry_msgs
  nav_msgs
  Eigen3
)

install(TARGETS MPC_SE
 DESTINATION lib/${PROJECT_NAME})

install(DIRECTORY
launch
worlds
DESTINATION share/${PROJECT_NAME}/
)

# add_dependencies(nonlinear_mpc_node ${${PROJECT_NAME}_EXPORTED_TARGETS})
