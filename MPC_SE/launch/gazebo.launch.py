# Copyright (c) 2021 Juan Miguel Jimeno
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http:#www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument, ExecuteProcess, IncludeLaunchDescription
from launch.substitutions import LaunchConfiguration, Command, PathJoinSubstitution
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch_ros.actions import Node
from launch_ros.substitutions import FindPackageShare
from ament_index_python.packages import get_package_share_directory
import math


def generate_launch_description():
    use_sim_time = True
    
    num = 2

    world_path = PathJoinSubstitution(
        [FindPackageShare("MPC_SE"), "worlds", "empty.world"]
    )

    description_launch_path = PathJoinSubstitution(
        [FindPackageShare('linorobot2_description'), 'launch', 'description.launch.py']
    )
    
    tello_description_path = get_package_share_directory('linorobot2_description')
    urdf_path = os.path.join(tello_description_path, 'urdf', 'robots', 'mecanum.urdf.xacro')
    
    entities = [
        ExecuteProcess(
            cmd=['gazebo', '--verbose', '-s', 'libgazebo_ros_factory.so',  '-s', 'libgazebo_ros_init.so', world_path],
            output='screen'
        ),        

        Node(
            package='linorobot2_gazebo',
            executable='command_timeout.py',
            name='command_timeout'
        ),

        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(description_launch_path),
            launch_arguments={
                'use_sim_time': str(use_sim_time),
                'publish_joints': 'false',
            }.items()
        )
    ]
    
    for i in range(num):
        entities.extend([
        Node(
            package='gazebo_ros',
            executable='spawn_entity.py',
            name='urdf_spawner',
            output='screen',
            arguments=["-topic", "robot_description", "-entity", "linorobot"+str(i+1), "-robot_namespace", "Vehicle"+str(i+1), "-x", str(0), "-y", str(-20*i)]# "-x", str(-10+10*math.cos(2*i*math.pi/num)), "-y", str(-10+10*math.sin(2*i*math.pi/num))]
        ),])

    return LaunchDescription(entities)

#sources: 
#https://navigation.ros.org/setup_guides/index.html#
#https://answers.ros.org/question/374976/ros2-launch-gazebolaunchpy-from-my-own-launch-file/
#https://github.com/ros2/rclcpp/issues/940
