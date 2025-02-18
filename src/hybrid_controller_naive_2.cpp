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

#include <iostream>
#include <fstream> 


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
void HybridController::topic_callback(const std::shared_ptr<franka_msgs::msg::FrankaRobotState> msg) {
  O_F_ext_hat_K = convertToStdArray(msg->o_f_ext_hat_k);
  arrayToMatrix(O_F_ext_hat_K, O_F_ext_hat_K_M); 
}

void HybridController::calculating_environmental_force_norm() {
    //F_norm_current = sqrt(pow(O_F_ext_hat_K[0], 2) + pow(O_F_ext_hat_K[1], 2) + pow(O_F_ext_hat_K[2], 2));
    F_norm_current = F_ext.head(3).norm();
}// maybe t_f_ext_hat_k...idk


void HybridController::calculating_position_norm() {
    //std::array<double, 16> pos_end = franka_robot_model_->getPoseMatrix(franka::Frame::kEndEffector);
    X_norm_current = error.head(3).norm();
  } 



void HybridController::update_environmental_stiffness(){ //calculation the environmental stiffness
        F_env[0] = F_env[1];
        X_env[0] = X_env[1];
        //k_e_temp[0] = k_e_temp[1];

        F_env[1] = F_norm_current;
        X_env[1] = X_norm_current;
        //only a easy moving average
        //k_e_temp[1] = (F_env[1]  / X_env[1] ) * 0.5 + 0.5 * k_e; 


        //k_e = (F_env[1] - F_env[0])/(X_env[1] - X_env[0] + 0.001);
        k_e = F_env[1]/ (X_env[1]+ 0.001);

        delta_k_e = k_e - k_e_last;
        delta_k_e = std::max(-max_k_e, std::min(delta_k_e, max_k_e));
        k_e = k_e_last + delta_k_e;
        k_e_last = k_e; //avoid great change in k_e

      
        if(k_e < 0){
          k_e = -k_e;
        }
        //k_e = F_norm_current / X_norm_current;
}





void HybridController::calculating_n(){
        if(k_e <= k_e_adm){
            n = 1.0; 
        }

        else if (k_e >= k_e_imp){
            n = 0.0;
        }
        else {
            n = 1.0 * (k_e_imp - k_e)/(k_e_imp - k_e_adm);
        }
} //n = 1 => total admittance; n = 0 => total impedance

void write_to_file(const std::string& filename, double value) { //output n and F_ext to external file
    std::ofstream output_file(filename, std::ios::app);
    if (!output_file) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }
    output_file << value << std::endl;
    output_file.close();
}

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
  //x_d_orientation_quat.normalize();


  x_d.head(3) << position_d_;
  x_d.tail(3) << x_d_orientation_quat.toRotationMatrix().eulerAngles(0, 1, 2); //might be wrong
  
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



void HybridController::reference_pose_callback(const geometry_msgs::msg::Pose::SharedPtr msg) //desired position , to be changed
{
    // Handle the incoming pose message
    std::cout << "received reference posistion as " <<  msg->position.x << ", " << msg->position.y << ", " << msg->position.z << std::endl;
    position_d_target_ << msg->position.x, msg->position.y, msg->position.z;
    orientation_d_target_.coeffs() << msg->orientation.x, msg->orientation.y, msg->orientation.z, msg->orientation.w;
    //orientation_d_target_.normalize();




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
  calculating_environmental_force_norm();
  calculating_position_norm();
  update_environmental_stiffness();
  calculating_n();

  Theta = Lambda;
  D = 2.05* K.cwiseSqrt() * Lambda.diagonal().cwiseSqrt().asDiagonal(); // Admittance control law to compute desired trajectory
  w = jacobian * dq_; // cartesian velocity


  //calculating acceleration
  acceleration = (w - w_last)/dt;
  w_last = w;

  //calculating error_dot
  error_dot = (error-error_last)/dt;
  error_last = error;

  



  error.head(3) << position - position_d_;
  if (orientation_d_.coeffs().dot(orientation.coeffs()) < 0.0) {
    orientation.coeffs() << -orientation.coeffs();
  }
  Eigen::Quaterniond error_quaternion(orientation.inverse() * orientation_d_);
  error.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
  error.tail(3) << -transform.rotation() * error.tail(3);
  
  Lambda = (jacobian * M.inverse() * jacobian.transpose()).inverse();
  //Force Updates
  F_ext = 0.9 * F_ext + 0.1 * O_F_ext_hat_K_M; // noFiltering
  F_ext.head(3) = -F_ext.head(3);


 // Convert x_d.tail(3) (Euler angles) to a quaternion
 Eigen::Matrix<double, 6, 1> virtual_error = Eigen::MatrixXd::Zero(6, 1);
  x_d_orientation_quat = Eigen::AngleAxisd(x_d.tail(3)(0), Eigen::Vector3d::UnitX())
                        * Eigen::AngleAxisd(x_d.tail(3)(1), Eigen::Vector3d::UnitY())
                        * Eigen::AngleAxisd(x_d.tail(3)(2), Eigen::Vector3d::UnitZ());


  if (orientation_d_.coeffs().dot(x_d_orientation_quat.coeffs()) < 0.0) {
      x_d_orientation_quat.coeffs() << -x_d_orientation_quat.coeffs();
  }
Eigen::Quaterniond virtual_error_quat = (x_d_orientation_quat.inverse() * orientation_d_);

// Normalize and check small norm
if (virtual_error_quat.norm() > 1e-6) {
    virtual_error_quat.normalize();
    // Use angle-axis representation to extract error vector
    Eigen::AngleAxisd angle_axis_virtual_error(virtual_error_quat);
    Eigen::Vector3d rotation_virtual_error_vector = angle_axis_virtual_error.angle() * angle_axis_virtual_error.axis();
    // Transform error into target frame
    virtual_error.tail(3) = -transform.rotation() * rotation_virtual_error_vector;
} else {
    virtual_error.tail(3).setZero();
}


  //virtual_error.tail(3) = -x_d_orientation_quat.toRotationMatrix() * virtual_error_quat.vec(); // vec for rotational part
  virtual_error.head(3) = x_d.head(3) - position_d_; //linear erro

  //F_ext?
  x_ddot_d = Lambda.inverse() * ( /*0 * */ F_ext - D * x_dot_d - K * virtual_error);
  // Integrate once to get velocities
  x_dot_d += x_ddot_d * dt;
  // Calculate the angle and axis for the quaternion rotation
  Eigen::Vector3d angular_displacement = x_dot_d.tail(3) * dt; // Angular displacement vector
  // Create the rotation increment quaternion
  Eigen::Quaterniond rotation_increment = Eigen::Quaterniond(Eigen::AngleAxisd(angular_displacement.norm(), angular_displacement.normalized()));

  x_d_orientation_quat = (rotation_increment * x_d_orientation_quat).normalized(); // Normalize to avoid drift
  //Convert final x_d_orientation_quat back to Euler angles if needed
  x_d.tail(3) = x_d_orientation_quat.toRotationMatrix().eulerAngles(0, 1, 2);
  x_d.head(3) += x_dot_d.head(3) * dt; 
  //x_d.head(3) << 0,0,0;
  error_pd.head(3) << position - x_d.head(3);
  //error_pd.head(3) <<0,0,0;
  if (x_d_orientation_quat.coeffs().dot(orientation.coeffs()) < 0.0) {
  orientation.coeffs() << -orientation.coeffs();
  }
    
  //error here
  
  
  //x - x_d
  Eigen::Quaterniond error_quaternion_current (orientation.inverse() * x_d_orientation_quat);
  error_pd.tail(3) << error_quaternion_current.x(), error_quaternion_current.y(), error_quaternion_current.z();
  error_pd.tail(3) << -transform.rotation() * error_pd.tail(3); 
  

  switch (mode_)
  {
  case 1:

    if(counter_ == 0){
      n_current = n;
    } //make sure that within every loop the n is constant

    

    if(counter < (20 - n_current * 20)) { // Impedance control

      /*if (smooth_adm_imp == 1){
        pos_d_switching.head(3) << position_d_;
        pos_d_switching.tail(3) << orientation_d_.toRotationMatrix().eulerAngles(0, 1, 2);

        pos_switching_imp = K.inverse() * (Kp * error_pd + Kd *w - D*w) + pos_d_switching;
        
        error.head(3) << pos_switching_imp.head(3) - position_d_;
        orientation = Eigen::AngleAxisd(pos_switching_imp[4], Eigen::Vector3d::UnitZ()) *
                          Eigen::AngleAxisd(pos_switching_imp[5], Eigen::Vector3d::UnitY()) *
                          Eigen::AngleAxisd(pos_switching_imp[6], Eigen::Vector3d::UnitX());
        Eigen::Quaterniond error_quaternion(orientation.inverse() * orientation_d_);
        error.tail(3) << error_quaternion.x(), error_quaternion.y(), error_quaternion.z();
        error.tail(3) << -transform.rotation() * error.tail(3);


        //error = Lambda.inverse() * (Kp * error_pd + Kd * w - D * w);

        smooth_adm_imp = 0;
      } */
      Eigen::Matrix<double, 6, 1> F_impedance_current = -1 * (D * (jacobian * dq_) + K * error);
      //F_impedance = (1 - n_current) * F_impedance_current + n_current * F_impedance_last; 
      desired_F = F_impedance_current; //desired_F for controlling the settling time, no impact on output
      //double blend_weight = 0.5 * (1 - cos(n_current * M_PI));
      double blend_weight = 0.3;
      F_output = blend_weight * F_impedance_current + (1 - blend_weight) * F_last;
      //delta_F_impedance = F_impedance_current - F_impedance_last;
      //delta_F_impedance = delta_F_impedance.cwiseMax(-max_F_impedance).cwiseMin(max_F_impedance);
      //F_impedance = F_impedance_last + delta_F_impedance;
      F_last = F_output;
      ++counter;
      ++counter_;
      flag_imp = true;
      smooth_imp_adm = 1;

    } else if(counter < 20) { // Admittance control

       if(smooth_imp_adm == 1){

        pos_switching_adm.head(3) << position;
        pos_switching_adm.tail(3) << orientation.toRotationMatrix().eulerAngles(0, 1, 2);
        x_d = pos_switching_adm - Kp.inverse() * ((D * (jacobian * dq_) + K * error) - Kd * w);
        //x - x_d
        /*error_pd.head(3) << position - x_d.head(3);
        Eigen::Quaterniond x_d_quat = Eigen::AngleAxisd(x_d[4], Eigen::Vector3d::UnitZ()) *
                          Eigen::AngleAxisd(x_d[5], Eigen::Vector3d::UnitY()) *
                          Eigen::AngleAxisd(x_d[6], Eigen::Vector3d::UnitX());
        Eigen::Quaterniond error_quaternion_x_d(orientation.inverse()*x_d_quat );
        error_pd.tail(3) << error_quaternion_x_d.x(), error_quaternion_x_d.y(), error_quaternion_x_d.z();
        error_pd.tail(3) << -transform.rotation() * error_pd.tail(3); */

        error_pd = Kp.inverse() * (D * w + K * error - Kd * w);
        x_dot_d =  Kp.inverse() * (D * acceleration + K * error_dot - Kd * w) + w;

        /*
        pos_switching_adm.head(3) << position;
        pos_switching_adm.tail(3) << orientation.toRotationMatrix().eulerAngles(0, 1, 2);
        Eigen::Matrix<double, 6, 1> A = F_ext - K * error - D * w;
        x_d = pos_switching_adm + Kp.inverse() * (Kd * w- F_ext + Lambda.inverse() * A);
        error_pd.head(3) << position - x_d.head(3);
        Eigen::Quaterniond x_d_quat = Eigen::AngleAxisd(x_d[4], Eigen::Vector3d::UnitZ()) *
                     Eigen::AngleAxisd(x_d[5], Eigen::Vector3d::UnitY()) *
                      Eigen::AngleAxisd(x_d[6], Eigen::Vector3d::UnitX());
        Eigen::Quaterniond error_quaternion_x_d(orientation.inverse() * x_d_quat );
        error_pd.tail(3) << error_quaternion_x_d.x(), error_quaternion_x_d.y(), error_quaternion_x_d.z();
        error_pd.tail(3) << -transform.rotation() * error_pd.tail(3);
        virtual_error.head(3) = x_d.head(3) - position_d_;
        */
        //x_dot_d.setZero();//naive Initialization 

        smooth_imp_adm = 0;
       }  
      
      Eigen::Matrix<double, 6, 1> F_admittance_current = - Kp * error_pd - Kd * w;
      //F_impedance = n_current * F_impedance_current + (1 - n_current) * F_impedance_last;
      desired_F = F_admittance_current;
      double blend_weight = 0.3;
      F_output = blend_weight * F_admittance_current + (1 - blend_weight) * F_last;

    
      //delta_F_impedance = F_admittance_current - F_impedance_last;
      //delta_F_impedance = delta_F_impedance.cwiseMax(-max_F_impedance).cwiseMin(max_F_impedance);
      //F_impedance = F_impedance_last + delta_F_impedance;
      F_last = F_output;
      ++counter;
      ++counter_;
      flag_imp = false;
      smooth_adm_imp = 1;
    } else {
      counter = 0;
      counter_ = 0;
    }
    //F_output.tail(3) = Eigen::Vector3d(0, 0, 0);//test


    
  break;

  case 2:
    Theta = T*Lambda;
    F_output = -1*(Lambda * Theta.inverse() - IDENTITY) * F_ext;
    break;
  
  default:
    break;
  }

   
  write_to_file("ke_data.txt", k_e);
  write_to_file("n.txt", n);
  write_to_file("output_x", F_output[1]);
  write_to_file("desired_F_x", desired_F[1]);
  write_to_file("Flag_imp", flag_imp);
  


    //F_impedance = - Kp * error_pd  - Kd * w;
    
  //in admittance, F_ext.head(3) = -F_ext.head(3), why?
 // F_ext = 0.9 * F_ext + 0.1 * O_F_ext_hat_K_M; //Filtering 
  /*I_F_error += dt * Sf* (F_contact_des - F_ext);
  F_cmd = Sf*(0.4 * (F_contact_des - F_ext)orientation_d_transpose_pinv); */
  Eigen::VectorXd tau_task(7), tau_nullspace(7), tau_d(7), tau_impedance(7);
  pseudoInverse(jacobian.transpose(), jacobian_transpose_pinv);

  tau_nullspace << (Eigen::MatrixXd::Identity(7, 7) -
                    jacobian.transpose() * jacobian_transpose_pinv) *
                    (nullspace_stiffness_ * config_control * (q_d_nullspace_ - q_) - //if config_control = true we control the whole robot configuration
                    (2.0 * sqrt(nullspace_stiffness_)) * dq_);  // if config control ) false we don't care about the joint position

  //calculate_tau_friction();
  tau_impedance = jacobian.transpose() * Sm * (F_output /*+ F_repulsion + F_potential*/)  /* + jacobian.transpose() * Sf * F_cmd*/;
  auto tau_d_placeholder = tau_impedance + tau_nullspace + coriolis  /*tau_friction*/; //add nullspace and coriolis components to desired torque
  tau_d << tau_d_placeholder;
  tau_d << saturateTorqueRate(tau_d, tau_J_d_M);  // Saturate torque rate to avoid discontinuities
  tau_J_d_M = tau_d;

  for (size_t i = 0; i < 7; ++i) {
    command_interfaces_[i].set_value(tau_d(i));
  }
  
  if (outcounter % 1000/update_frequency == 0){

    //std::cout << "loop_duration_imp [ms]: " << std::endl;
    //std::cout << loop_duration_imp << std::endl;
    //std::cout << "Target orientation: " << orientation_d_target_.toRotationMatrix().eulerAngles(0, 1, 2).transpose() << std::endl; //correct
    //std::cout << "x_d rotation is: " << x_d_orientation_quat.toRotationMatrix().eulerAngles(0, 1, 2).transpose() <<  std::endl; //wrong
    //std::cout << "Orientation virtual_error is: " << virtual_error.tail(3).transpose() <<  std::endl; //correct
    //std::cout << "Orientation Error_pd is: " << error.tail(3).transpose() <<  std::endl; //correct
    //std::cout << "Orientation_d_ is: " << orientation_d_.toRotationMatrix().eulerAngles(0, 1, 2).transpose() <<  std::endl;
    //std::cout << "Orientation is: " << orientation.toRotationMatrix().eulerAngles(0, 1, 2).transpose() <<  std::endl;
    std::cout << "Current n is: " << n <<  std::endl;
    //std::cout << "Current F_norm_current is: " << F_norm_current <<  std::endl;
    //std::cout << "Current X_norm_current is: " << X_norm_current <<  std::endl;
    std::cout << "Current k_e is: " << k_e <<  std::endl;
    std::cout << "Current tau_impedance is: " << tau_impedance <<  std::endl;
    std::cout << "Current counter is: " << counter <<  std::endl;
    std::cout << "If impedance " << flag_imp <<  std::endl;



    
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
