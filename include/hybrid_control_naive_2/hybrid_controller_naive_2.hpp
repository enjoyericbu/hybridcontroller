// Copyright (c) 2021 Franka Emika GmbH
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <unistd.h>
#include <thread>
#include <chrono>      
#include <vector>


#include "hybrid_control_naive_2/user_input_server.hpp"

#include <rclcpp/rclcpp.hpp>
#include "rclcpp/subscription.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigen>

#include <controller_interface/controller_interface.hpp>

#include <franka/model.h>
#include <franka/robot.h>
#include <franka/robot_state.h>

#include "franka_hardware/franka_hardware_interface.hpp"
#include <franka_hardware/model.hpp>

#include "franka_msgs/msg/franka_robot_state.hpp"
#include "franka_msgs/msg/errors.hpp"
#include "messages_fr3/srv/set_pose.hpp"

#include "franka_semantic_components/franka_robot_model.hpp"
#include "franka_semantic_components/franka_robot_state.hpp"
#include "geometry_msgs/msg/point.hpp"


#define IDENTITY Eigen::MatrixXd::Identity(6, 6)
#define POSITION_CONTROL 0

using CallbackReturn = rclcpp_lifecycle::node_interfaces::LifecycleNodeInterface::CallbackReturn;
using Vector7d = Eigen::Matrix<double, 7, 1>;

namespace hybrid_control_naive_2 {

class HybridController : public controller_interface::ControllerInterface {
public:

HybridController();

  [[nodiscard]] controller_interface::InterfaceConfiguration command_interface_configuration()
      const override;

  [[nodiscard]] controller_interface::InterfaceConfiguration state_interface_configuration()
      const override;

  controller_interface::return_type update(const rclcpp::Time& time,
                                           const rclcpp::Duration& period) override;
  controller_interface::CallbackReturn on_init() override;

  controller_interface::CallbackReturn on_configure(
      const rclcpp_lifecycle::State& previous_state) override;

  controller_interface::CallbackReturn on_activate(
      const rclcpp_lifecycle::State& previous_state) override;

  controller_interface::CallbackReturn on_deactivate(
      const rclcpp_lifecycle::State& previous_state) override;

    void setPose(const std::shared_ptr<messages_fr3::srv::SetPose::Request> request, 
    std::shared_ptr<messages_fr3::srv::SetPose::Response> response);


    void publishPositionTarget();//for publisher
    void schmitt_trigger();

      

 private:
    //Nodes
    rclcpp::Subscription<franka_msgs::msg::FrankaRobotState>::SharedPtr franka_state_subscriber = nullptr;
    rclcpp::Service<messages_fr3::srv::SetPose>::SharedPtr pose_srv_;
    rclcpp::Subscription<geometry_msgs::msg::Pose>::SharedPtr desired_pose_sub;




    //Functions
    void update_environmental_stiffness();
    void calculating_n(); 
    void calculating_environmental_force_norm();
    void calculating_position_norm();
    void topic_callback_force();
    void topic_callback_position(); //hybrid control
    void reference_pose_callback(const geometry_msgs::msg::Pose::SharedPtr msg);

 


    void calculate_tau_friction();   //friction compensation

    void topic_callback(const std::shared_ptr<franka_msgs::msg::FrankaRobotState> msg);
    void updateJointStates();
    void update_stiffness_and_references();
    void arrayToMatrix(const std::array<double, 6>& inputArray, Eigen::Matrix<double, 6, 1>& resultMatrix);
    void arrayToMatrix(const std::array<double, 7>& inputArray, Eigen::Matrix<double, 7, 1>& resultMatrix);
    Eigen::Matrix<double, 7, 1> saturateTorqueRate(const Eigen::Matrix<double, 7, 1>& tau_d_calculated, const Eigen::Matrix<double, 7, 1>& tau_J_d);  
    std::array<double, 6> convertToStdArray(const geometry_msgs::msg::WrenchStamped& wrench);
    //State vectors and matrices
    std::array<double, 7> q_subscribed;
    std::array<double, 7> tau_J_d = {0,0,0,0,0,0,0};
    std::array<double, 6> O_F_ext_hat_K = {0,0,0,0,0,0};
    Eigen::Matrix<double, 7, 1> q_subscribed_M;
    Eigen::Matrix<double, 7, 1> tau_J_d_M = Eigen::MatrixXd::Zero(7, 1);
    Eigen::Matrix<double, 6, 1> O_F_ext_hat_K_M = Eigen::MatrixXd::Zero(6,1);
    Eigen::Matrix<double, 7, 1> q_;
    Eigen::Matrix<double, 7, 1> dq_;
    Eigen::MatrixXd jacobian_transpose_pinv; 
    Eigen::Matrix<double, 6, 7> jacobian; 


    //Robot parameters
    const int num_joints = 7;
    const std::string state_interface_name_{"robot_state"};
    const std::string robot_name_{"fr3"};
    const std::string k_robot_state_interface_name{"robot_state"};
    const std::string k_robot_model_interface_name{"robot_model"};
    franka_hardware::FrankaHardwareInterface interfaceClass;
    std::unique_ptr<franka_semantic_components::FrankaRobotModel> franka_robot_model_;
    const double delta_tau_max_{1.0};
    const double dt = 0.001;
                
    //Impedance control variables              
    Eigen::Matrix<double, 6, 6> Lambda = IDENTITY;                                           // operational space mass matrix
    Eigen::Matrix<double, 6, 6> Sm = IDENTITY;                                               // task space selection matrix for positions and rotation
    Eigen::Matrix<double, 6, 6> Sf = Eigen::MatrixXd::Zero(6, 6);                            // task space selection matrix for forces
    Eigen::Matrix<double, 6, 6> K =  (Eigen::MatrixXd(6,6) << 250,   0,   0,   0,   0,   0,
                                                                0, 250,   0,   0,   0,   0,
                                                                0,   0, 250,   0,   0,   0,  // impedance stiffness term
                                                                0,   0,   0, 30,   0,   0,
                                                                0,   0,   0,   0, 30,   0,
                                                                0,   0,   0,   0,   0,  12).finished() ;

    Eigen::Matrix<double, 6, 6> D =  Eigen::MatrixXd(6,6) ;// impedance damping term

    Eigen::Matrix<double, 6, 6> K_imp =  (Eigen::MatrixXd(6,6) << 800,   0,   0,   0,   0,   0,
                                                                0, 800,   0,   0,   0,   0,
                                                                0,   0, 800,   0,   0,   0,  // impedance stiffness term
                                                                0,   0,   0, 80,   0,   0,
                                                                0,   0,   0,   0, 80,   0,
                                                                0,   0,   0,   0,   0,  18).finished() ;

                                                      

    // Eigen::Matrix<double, 6, 6> D =  (Eigen::MatrixXd(6,6) <<  30,   0,   0,   0,   0,   0,
    //                                                             0,  30,   0,   0,   0,   0,
    //                                                             0,   0,  30,   0,   0,   0,  // impedance damping term
    //                                                             0,   0,   0,  18,   0,   0,
    //                                                             0,   0,   0,   0,  18,   0,
    //                                                             0,   0,   0,   0,   0,   9).finished();
    Eigen::Matrix<double, 6, 6> Theta = IDENTITY;
    Eigen::Matrix<double, 6, 6> T = (Eigen::MatrixXd(6,6) <<       10,   0,   0,   0,   0,   0,
                                                                   0,   10,   0,   0,   0,   0,
                                                                   0,   0,   1,   0,   0,   0,  // Inertia term
                                                                   0,   0,   0,   10,   0,   0,
                                                                   0,   0,   0,   0,   10,   0,
                                                                   0,   0,   0,   0,   0,   1).finished();  
                                                                   
    /*Eigen::Matrix<double, 6, 6> T = (Eigen::MatrixXd(6,6) <<       1,   0,   0,   0,   0,   0,
                                                                   0,   1,   0,   0,   0,   0,
                                                                   0,   0,   2.5,   0,   0,   0,  // Inertia term,before
                                                                   0,   0,   0,   1,   0,   0,
                                                                   0,   0,   0,   0,   1,   0,
                                                                   0,   0,   0,   0,   0,   2.5).finished();  */                                                                                                            // impedance inertia term

    Eigen::Matrix<double, 6, 6> cartesian_stiffness_target_;                                 // impedance damping term
    Eigen::Matrix<double, 6, 6> cartesian_damping_target_;                                   // impedance damping term
    Eigen::Matrix<double, 6, 6> cartesian_inertia_target_;                                   // impedance damping term
    Eigen::Vector3d position_d_target_ = {0.5, 0.0, 0.5};
    Eigen::Vector3d rotation_d_target_ = {M_PI, 0.0, 0.0};
    Eigen::Quaterniond orientation_d_target_;
    Eigen::Vector3d position_d_;
    Eigen::Quaterniond orientation_d_; 
    Eigen::Quaterniond x_d_orientation_quat; //orientation quaternion of virtual spring system
    Eigen::Matrix<double, 6, 1> F_output;   
    Eigen::Matrix<double, 6, 1> F_contact_des = Eigen::MatrixXd::Zero(6, 1);                 // desired contact force
    Eigen::Matrix<double, 6, 1> F_contact_target = Eigen::MatrixXd::Zero(6, 1);              // desired contact force used for filtering
    Eigen::Matrix<double, 6, 1> F_ext = Eigen::MatrixXd::Zero(6, 1);                         // external forces
    Eigen::Matrix<double, 6, 1> F_cmd = Eigen::MatrixXd::Zero(6, 1);                         // commanded contact force
    Eigen::Matrix<double, 7, 1> q_d_nullspace_;
    Eigen::Matrix<double, 6, 1> error;                                                      // pose error (6d)
    Eigen::Matrix<double, 6, 1> error_pd;                                                    //pose error for control in admittance control
    double nullspace_stiffness_{0.001};
    double nullspace_stiffness_target_{0.001};
    //admittance controller parameter
    Eigen::Matrix<double, 6, 1> w;
    int control_mode;
    Eigen::Matrix<double, 6, 6> Kd; //  Will be initialized as critically damped  in constructor  
    double elapsed_time = 1;
    double Kp_multiplier = 1;
    Eigen::Matrix<double, 6, 6> Kp = (Eigen::MatrixXd(6,6) << 2000,   0,   0,   0,   0,   0,
                                                                0, 2000,   0,   0,   0,   0,
                                                                0,   0, 2000,   0,   0,   0,  // Inner Position Loop Controller Gains
                                                                0,   0,   0,  50,   0,   0,
                                                                0,   0,   0,   0,  50,   0,
                                                                0,   0,   0,   0,   0,  15).finished();

                                                                    // position and velocity outer controlle reference
    Eigen::Matrix<double, 6, 1> x_ddot_d;
    Eigen::Matrix<double, 6, 1> x_dot_d;
    Eigen::Matrix<double, 6, 1> x_d;

    Eigen::Matrix<double, 6, 1> F_impedance_current;
    Eigen::Matrix<double, 6, 1> F_admittance_current;

    

    //Logging
    int outcounter = 0;
    const int update_frequency = 2; //frequency for update outputs

    //Integrator
    Eigen::Matrix<double, 6, 1> I_error = Eigen::MatrixXd::Zero(6, 1);                      // pose error (6d)
    Eigen::Matrix<double, 6, 1> I_F_error = Eigen::MatrixXd::Zero(6, 1);                    // force error integral
    Eigen::Matrix<double, 6, 1> integrator_weights = 
      (Eigen::MatrixXd(6,1) << 75.0, 75.0, 75.0, 75.0, 75.0, 4.0).finished();
    Eigen::Matrix<double, 6, 1> max_I = 
      (Eigen::MatrixXd(6,1) << 30.0, 30.0, 30.0, 50.0, 50.0, 2.0).finished();

   
  
    std::mutex position_and_orientation_d_target_mutex_;

    //Flags
    bool config_control = false;           // sets if we want to control the configuration of the robot in nullspace
    bool do_logging = false;               // set if we do log values

    //Filter-parameters
    double filter_params_{0.001};
    int mode_ = 4;



    //hybrid controller
    double X_norm_current = 0.0;
    double F_norm_current = 0.0;


    std::array<double, 2> F_env;
    std::array<double, 2> X_env;

    std::vector<double> k_e_history = {0.0,0.0,0.0,0.0,0.0};
    double k_e = 0;
    
    double max_k_e = 1.5;
    double delta_k_e = 0;

    //maybe to replace the current moving average with a more precise one
    //std::vector<double , 5> F_env = {0,0,0,0,0}; 
    //std::vector<double , 5> X_env = {0,0,0,0,0};
    //std::vector<double , 5> k_e = {0,0,0,0,0};

    rclcpp::Publisher<geometry_msgs::msg::Point>::SharedPtr position_publisher_;



    
    Eigen::Matrix<double, 7, 1> tau_admittance = Eigen::MatrixXd::Zero(7,1); // admittance torque
    Eigen::Matrix<double, 7, 1> tau_admittance_filtered = Eigen::MatrixXd::Zero(7,1); // admittance torque filtered
    Eigen::Matrix<double, 7, 1> tau_friction = Eigen::MatrixXd::Zero(7,1);
    Eigen::Matrix<double, 7, 1> tau_threshold = Eigen::MatrixXd::Zero(7,1);  //Creating and filtering a "fake" tau_admittance with own weights, optimized for friction compensation
    bool friction = false; // set if friciton compensation should be turned on
    Eigen::MatrixXd N; // nullspace projection matrix


    // friction compensation observer
    Eigen::Matrix<double, 7, 1> dz = Eigen::MatrixXd::Zero(7,1);
    Eigen::Matrix<double, 7, 1> z = Eigen::MatrixXd::Zero(7,1);
    Eigen::Matrix<double, 7, 1> f = Eigen::MatrixXd::Zero(7,1);
    Eigen::Matrix<double, 7, 1> g = Eigen::MatrixXd::Zero(7,1);
    Eigen::Matrix<double, 6,6> K_friction = IDENTITY; //impedance stiffness term for friction compensation
    Eigen::Matrix<double, 6,6> D_friction = IDENTITY; //impedance damping term for friction compensation
    const Eigen::Matrix<double, 7, 1> sigma_0 = (Eigen::VectorXd(7) << 76.95, 37.94, 71.07, 44.02, 21.32, 21.83, 53).finished();
    const Eigen::Matrix<double, 7, 1> sigma_1 = (Eigen::VectorXd(7) << 0.056, 0.06, 0.064, 0.073, 0.1, 0.0755, 0.000678).finished();
    
    //friction compensation model paramers (coulomb, viscous, stribeck)
    Eigen::Matrix<double, 7, 1> dq_s = (Eigen::VectorXd(7) << 0, 0, 0, 0.0001, 0, 0, 0.05).finished();
    Eigen::Matrix<double, 7, 1> static_friction = (Eigen::VectorXd(7) << 1.025412896, 1.259913793, 0.8380147058, 1.005214968, 1.2928, 0.41525, 0.5341655).finished();
    Eigen::Matrix<double, 7, 1> offset_friction = (Eigen::VectorXd(7) << -0.05, -0.70, -0.07, -0.13, -0.1025, 0.103, -0.02).finished();
    Eigen::Matrix<double, 7, 1> coulomb_friction = (Eigen::VectorXd(7) << 1.025412896, 1.259913793, 0.8380147058, 0.96, 1.2928, 0.41525, 0.5341655).finished();
    Eigen::Matrix<double, 7, 1> beta = (Eigen::VectorXd(7) << 1.18, 0, 0.55, 0.87, 0.935, 0.54, 0.45).finished();//component b of linear friction model (a + b*dq)

 

    double k_e_imp = 0.06;  
    double k_e_adm = 0.03;   

    double n = 0.0;
    double n_current = 0.0;

    int counter_ ;  //counter for the whole sys， to calculate env stiffness and n
    int counter = 0; // counter for the update loop
    std::vector<double> dynamic_n;
    std::string filename = "output.txt";


    int counter_smooth = 0; // for the smooth change bet. imp and adm, not used
    Eigen::Matrix<double, 6, 1> F_last;
    Eigen::Matrix<double, 6, 1> F_impedance_new;
    Eigen::Matrix<double, 6, 1> delta_F_impedance;
    Eigen::Matrix<double, 6, 1> desired_F;
    Eigen::Matrix<double, 6, 1> pos_switching_adm;
    Eigen::Matrix<double, 6, 1> pos_d_switching;
    Eigen::Matrix<double, 6, 1> pos_switching_imp;


    Eigen::Matrix<double, 6, 1> acceleration;
    Eigen::Matrix<double, 6, 1> w_last;
    Eigen::Matrix<double, 6, 1> error_dot;
    Eigen::Matrix<double, 6, 1> error_last;
    


    double max_F_impedance = 1.5;
    bool flag_biginterpolation = 0;
    bool smooth_imp_adm;  
    bool smooth_adm_imp;
    bool flag_imp = false;

    
    bool schmitt_adm = true;

  

};
}  // namespace cartesian_impedance_control
