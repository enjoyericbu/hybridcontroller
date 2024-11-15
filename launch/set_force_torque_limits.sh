#!/bin/bash

ros2 service call /service_server/set_force_torque_collision_behavior franka_msgs/srv/SetForceTorqueCollisionBehavior "{
  lower_torque_thresholds_nominal: [60.0, 60.0, 60.0, 6000.0, 12.0, 12.0, 12.0] * 100, 
  upper_torque_thresholds_nominal: [60.0, 60.0, 60.0, 60.0, 12.0, 12.0, 12.0] * 100,
  lower_force_thresholds_nominal: [60.0, 60.0, 60.0, 60.0, 60.0, 60.0] * 100,
  upper_force_thresholds_nominal: [60.0, 60.0, 60.0, 60.0, 60.0, 60.0] * 100
}" #0 added for every num


