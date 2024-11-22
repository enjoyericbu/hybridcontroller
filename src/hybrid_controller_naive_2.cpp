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

#include <hybrid_control_naive_2/hybrid_controller_naive_2.hpp>

#include <cassert>
#include <cmath>
#include <exception>
#include <string>

#include <Eigen/Eigen>

namespace {

template <class T, size_t N>
std::ostream& operator<<(std::ostream& ostream, const std::array<T, N>& array) {
  ostream << "[";
  std::copy(array.cbegin(), array.cend() - 1, std::ostream_iterator<T>(ostream, ","));
  std::copy(array.cend() - 1, array.cend(), std::ostream_iterator<T>(ostream));
  ostream << "]";
  return ostream;
}
}

namespace hybrid_control_naive_2 {
//to be added in hybrid_controller.hpptions are copied from the previous hybrid_control_naive package. To be aware that the
//switching process might have problems with continuity
void HybridController::update_environmental_stiffness(){ //calculation the environmental stiffness
    if (F_env.size() <= 1 || X_env.size() <= 1) {
        std::cout << "F_env or X_env has insufficient data to compute stiffness." << std::endl;
        return;
    }

    while (rclcpp::ok() && counter_ < F_env.size() - 1 && counter_ < X_env.size() - 1)
    {
        if (X_env[counter_] - X_env[counter_ - 1] == 0) {
            std::cout << "Division by zero encountered in environmental stiffness calculation." << std::endl;
            return;
        }
        F_env[counter_] = F_norm_current;
        X_env[counter_] = X_norm_current;
        //only a easy moving average
        k_e_temp[counter_] = ((F_env[counter_] - F_env[counter_ - 1]) / (X_env[counter_] - X_env[counter_-1])) * 0.5 + 0.5 * k_e[counter_]; 
        k_e[counter_] = (k_e_temp[counter_] + k_e_temp[counter_ - 1] + k_e_temp[counter_ -2]) / 3;
        counter_++;
    }
}


    void HybridController::calculating_n(){
        if (counter_ >= k_e.size()) {
            std::cout << "Index out of bounds: counter_ exceeds size of k_e." << std::endl;
            return;
        }

        else if(k_e[counter_] <= k_e_adm){
            n = 1; 
        }

        else if (k_e[counter_] >= k_e_imp){
            n = 0;
        }
        else {
            n = 1 * (k_e_imp - k_e[counter_])/(k_e_imp - k_e_adm);
        }
    }


void HybridController::calculating_environmental_force_norm() {
    F_norm_current = sqrt(pow(O_F_ext_hat_K[0], 2) + pow(O_F_ext_hat_K[1], 2) + pow(O_F_ext_hat_K[2], 2));
}// maybe t_f_ext_hat_k...idk


void HybridController::calculating_position_norm() {
    std::array<double, 16> pos_end = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
    X_norm_current = sqrt(pow(pos_end[3], 2) + pow(pos_end[7], 2) + pow(pos_end[11], 2));
} //O_T_EE: r11, r12, r13, x; ....y ....z ;0 0 0 1




/*void HybridController::topic_callback_position(const std::shared_ptr<franka_msgs::msg::FrankaRobotState> msg) {
  O_T_EE = convertToStdArray(msg->o_t_ee);
  //arrayToMatrix(O_F_ext_hat_K, O_F_ext_hat_K_M); 
}  */


  //from AdmittancecController, parameters to be added in hpp

  HybridController::HybridController()  { //if error, delete node

  elapsed_time = 0.0;
  Kp_multiplier = 1; // Initial multiplier for Kp
  Kp = Kp * Kp_multiplier;  // increasing Kp from 0.1 to 1 made robot far less compliant
  control_mode = POSITION_CONTROL; // sets control mode
  //input_control_mode = TARGET_POSITION; // sets position control mode
  D =  2* K.cwiseSqrt(); // set critical damping from the get go
  Kd = 2 * Kp.cwiseSqrt();

  }






void HybridController::update_stiffness_and_references(){
  //target by filtering
  /** at the moment we do not use dynamic reconfigure and control the robot via D, K and T **/
  //K = filter_params_ * cartesian_stiffness_target_ + (1.0 - filter_params_) * K;
  //D = filter_params_ * cartesian_damping_target_ + (1.0 - filter_params_) * D;
  nullspace_stiffness_ = filter_params_ * nullspace_stiffness_target_ + (1.0 - filter_params_) * nullspace_stiffness_;
  //std::lock_guard<std::mutex> position_d_target_mutex_lock(position_and_orientation_d_target_mutex_);
  position_d_ = filter_params_ * position_d_target_ + (1.0 - filter_params_) * position_d_;
  orientation_d_ = orientation_d_.slerp(filter_params_, orientation_d_target_);
  F_contact_des = 0.05 * F_contact_target + 0.95 * F_contact_des;
}


void HybridController::arrayToMatrix(const std::array<double,7>& inputArray, Eigen::Matrix<double,7,1>& resultMatrix)
{
 for(long unsigned int i = 0; i < 7; ++i){
     resultMatrix(i,0) = inputArray[i];
   }
}

void HybridController::arrayToMatrix(const std::array<double,6>& inputArray, Eigen::Matrix<double,6,1>& resultMatrix)
{
 for(long unsigned int i = 0; i < 6; ++i){
     resultMatrix(i,0) = inputArray[i];
   }
}

Eigen::Matrix<double, 7, 1> HybridController::saturateTorqueRate(
  const Eigen::Matrix<double, 7, 1>& tau_d_calculated,
  const Eigen::Matrix<double, 7, 1>& tau_J_d_M) {  
  Eigen::Matrix<double, 7, 1> tau_d_saturated{};
  for (size_t i = 0; i < 7; i++) {
  double difference = tau_d_calculated[i] - tau_J_d_M[i];
  tau_d_saturated[i] =
         tau_J_d_M[i] + std::max(std::min(difference, delta_tau_max_), -delta_tau_max_);
  }
  return tau_d_saturated;
}


inline void pseudoInverse(const Eigen::MatrixXd& M_, Eigen::MatrixXd& M_pinv_, bool damped = true) {
  double lambda_ = damped ? 0.2 : 0.0;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_, Eigen::ComputeFullU | Eigen::ComputeFullV);   
  Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType sing_vals_ = svd.singularValues();
  Eigen::MatrixXd S_ = M_;  // copying the dimensions of M_, its content is not needed.
  S_.setZero();

  for (int i = 0; i < sing_vals_.size(); i++)
     S_(i, i) = (sing_vals_(i)) / (sing_vals_(i) * sing_vals_(i) + lambda_ * lambda_);

  M_pinv_ = Eigen::MatrixXd(svd.matrixV() * S_.transpose() * svd.matrixU().transpose());
}


controller_interface::InterfaceConfiguration
HybridController::command_interface_configuration() const {
  controller_interface::InterfaceConfiguration config;
  config.type = controller_interface::interface_configuration_type::INDIVIDUAL;
  for (int i = 1; i <= num_joints; ++i) {
    config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/effort");
  }
  return config;
}


controller_interface::InterfaceConfiguration HybridController::state_interface_configuration()
  const {
  controller_interface::InterfaceConfiguration state_interfaces_config;
  state_interfaces_config.type = controller_interface::interface_configuration_type::INDIVIDUAL;

  for (int i = 1; i <= num_joints; ++i) {
    state_interfaces_config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/position");
    state_interfaces_config.names.push_back(robot_name_ + "_joint" + std::to_string(i) + "/velocity");
  }

  for (const auto& franka_robot_model_name : franka_robot_model_->get_state_interface_names()) {
    state_interfaces_config.names.push_back(franka_robot_model_name);
    std::cout << franka_robot_model_name << std::endl;
  }

  const std::string full_interface_name = robot_name_ + "/" + state_interface_name_;

  return state_interfaces_config;
}


CallbackReturn HybridController::on_init() {
   UserInputServer input_server_obj(&position_d_target_, &rotation_d_target_, &K, &D, &T);
   std::thread input_thread(&UserInputServer::main, input_server_obj, 0, nullptr);
   input_thread.detach();
   return CallbackReturn::SUCCESS;
}


CallbackReturn HybridController::on_configure(const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_ = std::make_unique<franka_semantic_components::FrankaRobotModel>(
  franka_semantic_components::FrankaRobotModel(robot_name_ + "/" + k_robot_model_interface_name,
                                               robot_name_ + "/" + k_robot_state_interface_name));
                                               
  try {
    rclcpp::QoS qos_profile(1); // Depth of the message queue
    qos_profile.reliability(RMW_QOS_POLICY_RELIABILITY_RELIABLE);
    franka_state_subscriber = get_node()->create_subscription<franka_msgs::msg::FrankaRobotState>(
    "franka_robot_state_broadcaster/robot_state", qos_profile, 
    std::bind(&HybridController::topic_callback, this, std::placeholders::_1));
    std::cout << "Succesfully subscribed to robot_state_broadcaster" << std::endl;
  }

  catch (const std::exception& e) {
    fprintf(stderr,  "Exception thrown during publisher creation at configure stage with message : %s \n",e.what());
    return CallbackReturn::ERROR;
    }


  RCLCPP_DEBUG(get_node()->get_logger(), "configured successfully");
  return CallbackReturn::SUCCESS;
}


CallbackReturn HybridController::on_activate(
  


  const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_->assign_loaned_state_interfaces(state_interfaces_);

    desired_pose_sub = get_node()->create_subscription<geometry_msgs::msg::Pose>(
        "admittance_controller/reference_pose", 
        10,  // Queue size
        std::bind(&HybridController::reference_pose_callback, this, std::placeholders::_1)
    );

  std::array<double, 16> initial_pose = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
  Eigen::Affine3d initial_transform(Eigen::Matrix4d::Map(initial_pose.data()));
  position_d_ = initial_transform.translation();
  orientation_d_ = Eigen::Quaterniond(initial_transform.rotation());
   //for admittance control
  x_d_orientation_quat.coeffs() << orientation_d_.coeffs();

  x_d.head(3) << position_d_;
  x_d.tail(3) << x_d_orientation_quat.toRotationMatrix().eulerAngles(0, 1, 2);
  std::cout << "position_d, orientation_d, on_activate is: " << position_d_.transpose() << " " << initial_transform.rotation().eulerAngles(0, 1, 2).transpose() <<  std::endl;    // Debugging
  std::cout << "position_d_target on activation is: " << position_d_target_.transpose() <<  std::endl;    // Debugging
  std::cout << "x_desired head on_activate is: " << x_d.head(3) <<  std::endl;    // Debugging
  std::cout << "x_desired tail on_activate is: " << x_d.tail(3) <<  std::endl;    // Debugging
  std::cout << "Completed Activation process" << std::endl;
  std::cout << "Completed Activation process" << std::endl;
  return CallbackReturn::SUCCESS;
}


controller_interface::CallbackReturn HybridController::on_deactivate(
  const rclcpp_lifecycle::State& /*previous_state*/) {
  franka_robot_model_->release_interfaces();
  return CallbackReturn::SUCCESS;
}

std::array<double, 6> HybridController::convertToStdArray(const geometry_msgs::msg::WrenchStamped& wrench) {
    std::array<double, 6> result;
    result[0] = wrench.wrench.force.x;
    result[1] = wrench.wrench.force.y;
    result[2] = wrench.wrench.force.z;
    result[3] = wrench.wrench.torque.x;
    result[4] = wrench.wrench.torque.y;
    result[5] = wrench.wrench.torque.z;
    return result;
}

void HybridController::topic_callback(const std::shared_ptr<franka_msgs::msg::FrankaRobotState> msg) {
  O_F_ext_hat_K = convertToStdArray(msg->o_f_ext_hat_k);
  arrayToMatrix(O_F_ext_hat_K, O_F_ext_hat_K_M); 
}


void HybridController::reference_pose_callback(const geometry_msgs::msg::Pose::SharedPtr msg) //desired position , to be changed
{
    // Handle the incoming pose message
    std::cout << "received reference posistion as " <<  msg->position.x << ", " << msg->position.y << ", " << msg->position.z << std::endl;
    position_d_target_ << msg->position.x, msg->position.y, msg->position.z;
    orientation_d_target_.coeffs() << msg->orientation.x, msg->orientation.y, msg->orientation.z, msg->orientation.w;
    



    // You can add more processing logic here
    // Update x_d to reflect the new reference poses
/*     x_d.head(3) = position_d_target_;  // New target position
    x_d.tail(3) << msg->orientation.x, msg->orientation.y, msg->orientation.z;  // New target orientation */
}






void HybridController::updateJointStates() {
  for (auto i = 0; i < num_joints; ++i) {
    const auto& position_interface = state_interfaces_.at(2 * i);
    const auto& velocity_interface = state_interfaces_.at(2 * i + 1);
    assert(position_interface.get_interface_name() == "position");
    assert(velocity_interface.get_interface_name() == "velocity");
    q_(i) = position_interface.get_value();
    dq_(i) = velocity_interface.get_value();
  }
}


controller_interface::return_type HybridController::update(const rclcpp::Time& /*time*/, const rclcpp::Duration& /*period*/) {  
  // if (outcounter == 0){
  // std::cout << "Enter 1 if you want to track a desired position or 2 if you want to use free floating with optionally shaped inertia" << std::endl;
  // std::cin >> mode_;
  // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  // std::cout << "Mode selected" << std::endl;
  // while (mode_ != 1 && mode_ != 2){
  //   std::cout << "Invalid mode, try again" << std::endl;
  //   std::cin >> mode_;
  // }
  // }



  std::array<double, 49> mass = franka_robot_model_->getMassMatrix(); //M(q)
  std::array<double, 7> coriolis_array = franka_robot_model_->getCoriolisForceVector();  //h(q,q_dot)
  std::array<double, 42> jacobian_array =  franka_robot_model_->getZeroJacobian(franka::Frame::kEndEffector); //J
  std::array<double, 16> pose = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector); //state_of_
  
  Eigen::Map<Eigen::Matrix<double, 7, 1>> coriolis(coriolis_array.data());
  Eigen::Map<Eigen::Matrix<double, 6, 7>> jacobian(jacobian_array.data());
  Eigen::Map<Eigen::Matrix<double, 7, 7>> M(mass.data());
  
  Eigen::Affine3d transform(Eigen::Matrix4d::Map(pose.data()));
  Eigen::Vector3d position(transform.translation());
  Eigen::Quaterniond orientation(transform.rotation());
  orientation_d_target_ = Eigen::AngleAxisd(rotation_d_target_[0], Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(rotation_d_target_[1], Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(rotation_d_target_[2], Eigen::Vector3d::UnitZ());
  updateJointStates(); 

  
  //std::cout << "Error outer loop is: " << error.transpose() <<  std::endl;
  // Now, align error.tail(3) to use virtual_error like position error
  // You want to use the virtual_error for rotational error (tail):
  // Set current state  
  Theta = Lambda;
  D = 2.05* K.cwiseSqrt() * Lambda.diagonal().cwiseSqrt().asDiagonal(); // Admittance control law to compute desired trajectory
  w = jacobian * dq_; // cartesian velocity

  


  error.head(3) << position - position_d_;   
 


  Lambda = (jacobian * M.inverse() * jacobian.transpose()).inverse();
  // Theta = T*Lambda;
  // F_impedance = -1*(Lambda * Theta.inverse() - IDENTITY) * F_ext;
  //Inertia of the robot

   error_pd.head(3) << position -x_d.head(3); //error for pd controller in admittance control
   //copied froom previous error

  //calculating the velocity of end effector for the admittance control:

//Force Updates
F_ext = 0.001 * F_ext + 0.999 * O_F_ext_hat_K_M; // noFiltering
F_ext.head(3) = -F_ext.head(3);


 // Convert x_d.tail(3) (Euler angles) to a quaternion
 Eigen::Matrix<double, 6, 1> virtual_error = Eigen::MatrixXd::Zero(6, 1);
  x_d_orientation_quat = Eigen::AngleAxisd(x_d.tail(3)(0), Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(x_d.tail(3)(1), Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(x_d.tail(3)(2), Eigen::Vector3d::UnitZ());


  if (orientation_d_.coeffs().dot(x_d_orientation_quat.coeffs()) < 0.0) {
      x_d_orientation_quat.coeffs() << -x_d_orientation_quat.coeffs();
  }

  Eigen::Quaterniond virtual_error_quat = orientation_d_.inverse() * x_d_orientation_quat;
  // Extract the virtual orientation error as a 3D vector (imaginary part for small-angle approximation)
  // Transform the error into the desired orientation frame
  virtual_error.tail(3) << virtual_error_quat.x(), virtual_error_quat.y(), virtual_error_quat.z();
  virtual_error.tail(3) << -x_d_orientation_quat.toRotationMatrix() * virtual_error.tail(3);
  //virtual_error.tail(3) = -x_d_orientation_quat.toRotationMatrix() * virtual_error_quat.vec(); // vec for rotational part
  virtual_error.head(3) = x_d.head(3) - position_d_; //linear error


  x_ddot_d = Lambda.inverse() * (-F_ext - D * x_dot_d - K * virtual_error);
  // Integrate once to get velocities
  x_dot_d += x_ddot_d * dt;
  // Calculate the angle and axis for the quaternion rotation
  Eigen::Vector3d angular_displacement = x_dot_d.tail(3) * dt; // Angular displacement vector
  // Create the rotation increment quaternion
  Eigen::Quaterniond rotation_increment = Eigen::Quaterniond(Eigen::AngleAxisd(angular_displacement.norm(), angular_displacement.normalized()));
  // Update the orientation quaternion case 1:
    //Theta = Lambda;
    //F_impedance = -1 * (D * (jacobian * dq_) + K * error /*+ I_error*/); 

      
        //also copied from the hybird_control_naive, to be added in hpp
        //auto loop_duration_imp = std::chrono::milliseconds(static_cast<int>(10 * (1.0 - n))) ; // the program will run for 10ms
        //auto loop_duration_adm = std::chrono::milliseconds(10) - loop_duration_imp; //unsure how to use the n from hybridcontrol.cpp
        // it's highly possible that the rob cannot follow the desired motion due to frequency.
        //std::cout <<loop_duration_imp<< std::endl;

        

  x_d_orientation_quat = (rotation_increment * x_d_orientation_quat).normalized(); // Normalize to avoid drift
  //Convert final x_d_orientation_quat back to Euler angles if needed
  x_d.tail(3) = x_d_orientation_quat.toRotationMatrix().eulerAngles(0, 1, 2);
  x_d.head(3) += x_dot_d.head(3) * dt; 
  
  //error_pd.head(3) << position - x_d.head(3);//error for pd controller in admittance control
  if (x_d_orientation_quat.coeffs().dot(orientation.coeffs()) < 0.0) {
  orientation.coeffs() << -orientation.coeffs();
  }
  Eigen::Quaterniond error_quaternion(orientation.inverse() * x_d_orientation_quat);
  error_pd.tail(3) << orientation.toRotationMatrix().eulerAngles(0, 1, 2) - x_d.tail(3);//, error_quaternion.y(), error_quaternion.z();

  //error_pd.tail(3) << transform.rotation() * error.tail(3);


  switch (mode_)
  {
  case 1:
    //Theta = Lambda;
    //F_impedance = -1 * (D * (jacobian * dq_) + K * error /*+ I_error*/); 

      
        //also copied from the hybird_control_naive, to be added in hpp
        //auto loop_duration_imp = std::chrono::milliseconds(static_cast<int>(10 * (1.0 - n))) ; // the program will run for 10ms
        //auto loop_duration_adm = std::chrono::milliseconds(10) - loop_duration_imp; //unsure how to use the n from hybridcontrol.cpp
        // it's highly possible that the rob cannot follow the desired motion due to frequency.
        //std::cout <<loop_duration_imp<< std::endl;



        // for historical reason the f output for impedance and admittance control are called F_impedance

        if(counter < (100 - n * 100)) { //impedance control 
          Theta = Lambda;
          F_impedance = -1 * (D * (jacobian * dq_) + K * error ); 
          error_pd.head(3) = position - (position + 
              (Kp.inverse() * ((Kd * (jacobian * dq_)) - F_ext)).head(3));  //x-k_d
          //F_impedance_last = F_impedance;
          ++counter;
        }

        /*else if(100 - counter >=5){
          flag_biginterpolation = 1;
        }

        else if(flag_biginterpolation == 1 && counter_smooth <= 5){
          F_impedance_new = - Kp * error_pd - Kd * (x_dot_d - w);
          double alpha = static_cast<double>(counter_smooth) / 5.0; 
          F_impedance = (1.0 - alpha) * F_impedance_last + alpha * F_impedance_new;
          ++counter_smooth;
          ++counter;
        } */
        
        
        else if(counter < 100){ //admittance control


        
          //- Kp * error_pd wrong , -Kd * w okay
        
          F_impedance = - Kp * error_pd  - Kd * w; //position controller
          
          //F_impedance = - Kp * error_pd - Kd * (x_dot_d - w); //position - x_d.head(3)
          //F_impedance = - Kd * (w - x_dot_d);  //velocity control
          //w = jacobian * dq_;
          ++counter;
                 
        }  //unsure if the robot can get the F_inpedance from the while loop
        else{
          //F_impedance = - Kp * error - Kd * w;  
          counter = 0;
          flag_biginterpolation = 0;
        }
      



        
    break;
  case 2:
    Theta = T*Lambda;
    F_impedance = -1*(Lambda * Theta.inverse() - IDENTITY) * F_ext;
    break;
  
  default:
    break;
  }

  //position_d_target_(1) += 0000.1 * cos( 2 * 3.1415926 * 0.0001 * counter);
  //position_d_target_(2) += 0000.1 * sin( 2 * 3.1415926 * 0.0001 * counter); // change the reference position

  //in admittance, F_ext.head(3) = -F_ext.head(3), why?
 // F_ext = 0.9 * F_ext + 0.1 * O_F_ext_hat_K_M; //Filtering 
  /*I_F_error += dt * Sf* (F_contact_des - F_ext);
  F_cmd = Sf*(0.4 * (F_contact_des - F_ext) + 0.9 * I_F_error + 0.9 * F_contact_des);*/
 //F_cmd in admittance not needed?
  Eigen::VectorXd tau_task(7), tau_nullspace(7), tau_d(7), tau_impedance(7);
  pseudoInverse(jacobian.transpose(), jacobian_transpose_pinv);

  tau_nullspace << (Eigen::MatrixXd::Identity(7, 7) -
                    jacobian.transpose() * jacobian_transpose_pinv) *
                    (nullspace_stiffness_ * config_control * (q_d_nullspace_ - q_) - //if config_control = true we control the whole robot configuration
                    (2.0 * sqrt(nullspace_stiffness_)) * dq_);  // if config control ) false we don't care about the joint position

  //calculate_tau_friction();
  tau_impedance = jacobian.transpose() * Sm * (F_impedance /*+ F_repulsion + F_potential*/)  /* + jacobian.transpose() * Sf * F_cmd*/;
  auto tau_d_placeholder = tau_impedance + tau_nullspace + coriolis  /*tau_friction*/; //add nullspace and coriolis components to desired torque
  //friction not needed?
  tau_d << tau_d_placeholder;
  tau_d << saturateTorqueRate(tau_d, tau_J_d_M);  // Saturate torque rate to avoid discontinuities
  tau_J_d_M = tau_d;

  for (size_t i = 0; i < 7; ++i) {
    command_interfaces_[i].set_value(tau_d(i));
  }
  
  if (outcounter % 1000/update_frequency == 0){
    std::cout << "F_ext_robot [N]" << std::endl;
    std::cout << O_F_ext_hat_K << std::endl;
    std::cout << O_F_ext_hat_K_M << std::endl;
    std::cout << "Lambda  Thetha.inv(): " << std::endl;
    std::cout << Lambda*Theta.inverse() << std::endl;
    std::cout << "tau_d" << std::endl;
    std::cout << tau_d << std::endl;
    std::cout << "--------" << std::endl;
    std::cout << tau_nullspace << std::endl;counter >=n && 
    std::cout << "--------" << std::endl;
    std::cout << tau_impedance << std::endl;
    std::cout << "--------" << std::endl;
    std::cout << coriolis << std::endl;
    std::cout << "Inertia scaling [m]: " << std::endl;
    std::cout << T << std::endl;
    std::cout << "current n for hybrid control: " << std::endl;
    std::cout << n << std::endl;
    //std::cout << "loop_duration_imp [ms]: " << std::endl;
    //std::cout << loop_duration_imp << std::endl;
  }
  outcounter++;
  update_stiffness_and_references();
  return controller_interface::return_type::OK;
}
}

// namespace cartesian_impedance_control
#include "pluginlib/class_list_macros.hpp"
// NOLINTNEXTLINE
PLUGINLIB_EXPORT_CLASS(hybrid_control_naive_2::HybridController,
                       controller_interface::ControllerInterface)