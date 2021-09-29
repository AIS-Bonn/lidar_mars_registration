/**
BSD 3-Clause License

This file is part of the Basalt project.
https://gitlab.com/VladyslavUsenko/basalt-headers.git

Copyright (c) 2019, Vladyslav Usenko and Nikolaus Demmel.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@file
@brief Type definitions for IMU preintegration
*/

#pragma once

#include <basalt/utils/sophus_utils.hpp>
#include <memory>
#include <basalt/utils/assert.h>

namespace basalt {

constexpr size_t POSE_SIZE = 6;  ///< Dimentionality of the pose state
constexpr size_t POSE_VEL_SIZE =
    9;  ///< Dimentionality of the pose-velocity state
constexpr size_t POSE_VEL_BIAS_SIZE =
    15;  ///< Dimentionality of the pose-velocity-bias state

/// @brief State that consists of SE(3) pose at a certain time.
struct PoseState {
  using VecN = Eigen::Matrix<double, POSE_SIZE, 1>;

  /// @brief Default constructor with Identity pose and zero timestamp.
  PoseState() { t_ns = 0; }

  /// @brief Constructor with timestamp and pose.
  ///
  /// @param t_ns timestamp of the state in nanoseconds
  /// @param T_w_i transformation from the body frame to the world frame
  PoseState(int64_t t_ns, const Sophus::SE3d& T_w_i)
      : t_ns(t_ns), T_w_i(T_w_i) {}

  /// @brief Apply increment to the pose
  ///
  /// For details see \ref incPose
  /// @param[in] inc 6x1 increment vector
  void applyInc(const VecN& inc) { incPose(inc, T_w_i); }

  /// @brief Apply increment to the pose
  ///
  /// The incremernt vector consists of translational and rotational parts \f$
  /// [\upsilon, \omega]^T \f$. Given the current pose \f$ R \in
  /// SO(3), p \in \mathbb{R}^3\f$ the updated pose is: \f{align}{ R' &=
  /// \exp(\omega) R
  /// \\ p' &= p + \upsilon
  /// \f}
  ///  The increment is consistent with \ref
  /// Se3Spline::applyInc.
  ///
  /// @param[in] inc 6x1 increment vector
  /// @param[in,out] T the pose to update
  inline static void incPose(const Sophus::Vector6d& inc, Sophus::SE3d& T) {
    T.translation() += inc.head<3>();
    T.so3() = Sophus::SO3d::exp(inc.tail<3>()) * T.so3();
  }

  int64_t t_ns;        ///< timestamp of the state in nanoseconds
  Sophus::SE3d T_w_i;  ///< pose of the state
};

/// @brief State that consists of SE(3) pose and linear velocity at a certain
/// time.
struct PoseVelState : public PoseState {
  using VecN = Eigen::Matrix<double, POSE_VEL_SIZE, 1>;

  /// @brief Default constructor with Identity pose and zero other values.
  PoseVelState() { vel_w_i.setZero(); };

  /// @brief Constructor with timestamp, pose and linear velocity.
  ///
  /// @param t_ns timestamp of the state in nanoseconds
  /// @param T_w_i transformation from the body frame to the world frame
  /// @param vel_w_i linear velocity in world coordinate frame
  PoseVelState(int64_t t_ns, const Sophus::SE3d& T_w_i,
               const Eigen::Vector3d& vel_w_i)
      : PoseState(t_ns, T_w_i), vel_w_i(vel_w_i) {}

  /// @brief Apply increment to the state.
  ///
  /// For pose see \ref incPose. For velocity simple addition.
  /// @param[in] inc 9x1 increment vector [trans, rot, vel]
  void applyInc(const VecN& inc) {
    PoseState::applyInc(inc.head<6>());
    vel_w_i += inc.tail<3>();
  }

  /// @brief Compute difference to other state.
  ///
  /// ```
  ///      PoseVelState::VecN inc;
  ///      PoseVelState p0, p1;
  ///      // Initialize p0 and inc
  ///      p1 = p0;
  ///      p1.applyInc(inc);
  ///      p0.diff(p1) == inc; // Should be true.
  /// ```
  /// @param other state to compute difference.
  VecN diff(const PoseVelState& other) const {
    VecN res;
    res.segment<3>(0) = other.T_w_i.translation() - T_w_i.translation();
    res.segment<3>(3) = (other.T_w_i.so3() * T_w_i.so3().inverse()).log();
    res.tail<3>() = other.vel_w_i - vel_w_i;
    return res;
  }

  Eigen::Vector3d vel_w_i;  ///< Linear velocity of the state
};

/// @brief State that consists of SE(3) pose, linear velocity, gyroscope and
/// accelerometer biases at a certain time.
struct PoseVelBiasState : public PoseVelState {
  using Ptr = std::shared_ptr<PoseVelBiasState>;
  using VecN = Eigen::Matrix<double, POSE_VEL_BIAS_SIZE, 1>;

  /// @brief Default constructor with Identity pose and zero other values.
  PoseVelBiasState() {
    bias_gyro.setZero();
    bias_accel.setZero();
  };

  /// @brief Constructor with timestamp, pose, linear velocity, gyroscope and
  /// accelerometer biases.
  ///
  /// @param t_ns timestamp of the state in nanoseconds
  /// @param T_w_i transformation from the body frame to the world frame
  /// @param vel_w_i linear velocity in world coordinate frame
  /// @param bias_gyro gyroscope bias
  /// @param bias_accel accelerometer bias
  PoseVelBiasState(int64_t t_ns, const Sophus::SE3d& T_w_i,
                   const Eigen::Vector3d& vel_w_i,
                   const Eigen::Vector3d& bias_gyro,
                   const Eigen::Vector3d& bias_accel)
      : PoseVelState(t_ns, T_w_i, vel_w_i),
        bias_gyro(bias_gyro),
        bias_accel(bias_accel) {}

  /// @brief Apply increment to the state.
  ///
  /// For pose see \ref incPose. For velocity and biases simple addition.
  /// @param[in] inc 15x1 increment vector [trans, rot, vel, bias_gyro,
  /// bias_accel]
  void applyInc(const VecN& inc) {
    PoseVelState::applyInc(inc.head<9>());
    bias_gyro += inc.segment<3>(9);
    bias_accel += inc.segment<3>(12);
  }

  /// @brief Compute difference to other state.
  ///
  /// ```
  ///      PoseVelBiasState::VecN inc;
  ///      PoseVelBiasState p0, p1;
  ///      // Initialize p0 and inc
  ///      p1 = p0;
  ///      p1.applyInc(inc);
  ///      p0.diff(p1) == inc; // Should be true.
  /// ```
  /// @param other state to compute difference.
  VecN diff(const PoseVelBiasState& other) const {
    VecN res;
    res.segment<9>(0) = PoseVelState::diff(other);
    res.segment<3>(9) = other.bias_gyro - bias_gyro;
    res.segment<3>(12) = other.bias_accel - bias_accel;
    return res;
  }

  Eigen::Vector3d bias_gyro;   ///< Gyroscope bias
  Eigen::Vector3d bias_accel;  ///< Accelerometer bias
};

/// @brief Timestamped gyroscope and accelerometer measurements.
struct ImuData {
  using Ptr = std::shared_ptr<ImuData>;

  int64_t t_ns;           ///< timestamp in nanoseconds
  Eigen::Vector3d accel;  ///< Accelerometer measurement
  Eigen::Vector3d gyro;   ///< Gyroscope measurement

  /// @brief Default constructor with zero measurements.
  ImuData() {
    accel.setZero();
    gyro.setZero();
  }

};

namespace constants {
static const Eigen::Vector3d g(0, 0, -9.81);
static const Eigen::Vector3d g_dir(0, 0, -1);
}  // namespace constants

struct PoseVelBiasStateWithLin {
  using VecN = PoseVelBiasState::VecN;

  PoseVelBiasStateWithLin() {
    linearized = false;
    delta.setZero();
  }

  PoseVelBiasStateWithLin(int64_t t_ns, const Sophus::SE3d& T_w_i,
                          const Eigen::Vector3d& vel_w_i,
                          const Eigen::Vector3d& bias_gyro,
                          const Eigen::Vector3d& bias_accel, bool linearized)
      : linearized(linearized),
        state_linearized(t_ns, T_w_i, vel_w_i, bias_gyro, bias_accel) {
    delta.setZero();
    state_current = state_linearized;
  }

  PoseVelBiasStateWithLin(const PoseVelBiasState& other)
      : linearized(false), state_linearized(other) {
    delta.setZero();
    state_current = other;
  }

  void setLinFalse() {
    linearized = false;
    delta.setZero();
  }

  void setLinTrue() {
    linearized = true;
    BASALT_ASSERT(delta.isApproxToConstant(0));
    state_current = state_linearized;
  }

  void applyInc(const VecN& inc) {
    if (!linearized) {
      state_linearized.applyInc(inc);
    } else {
      delta += inc;
      state_current = state_linearized;
      state_current.applyInc(delta);
    }
  }

  inline const PoseVelBiasState& getState() const {
    if (!linearized) {
      return state_linearized;
    } else {
      return state_current;
    }
  }

  inline const PoseVelBiasState& getStateLin() const {
    return state_linearized;
  }

  inline bool isLinearized() const { return linearized; }
  inline const VecN& getDelta() const { return delta; }
  inline int64_t getT_ns() const { return state_linearized.t_ns; }

  friend struct PoseStateWithLin;

  inline void backup() {
    backup_delta = delta;
    backup_state_linearized = state_linearized;
    backup_state_current = state_current;
  }

  inline void restore() {
    delta = backup_delta;
    state_linearized = backup_state_linearized;
    state_current = backup_state_current;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 private:
  bool linearized;
  VecN delta;
  PoseVelBiasState state_linearized, state_current;

  VecN backup_delta;
  PoseVelBiasState backup_state_linearized, backup_state_current;
};

struct PoseStateWithLin {
  using VecN = PoseState::VecN;

  PoseStateWithLin() {
    linearized = false;
    delta.setZero();
  }

  PoseStateWithLin(int64_t t_ns, const Sophus::SE3d& T_w_i,
                   bool linearized = false)
      : linearized(linearized), pose_linearized(t_ns, T_w_i) {
    delta.setZero();
    T_w_i_current = T_w_i;
  }

  PoseStateWithLin(const PoseVelBiasStateWithLin& other)
      : linearized(other.linearized),
        delta(other.delta.head<6>()),
        pose_linearized(other.state_linearized.t_ns,
                        other.state_linearized.T_w_i) {
    T_w_i_current = pose_linearized.T_w_i;
    PoseState::incPose(delta, T_w_i_current);
  }

  void setLinTrue() {
    linearized = true;
    BASALT_ASSERT(delta.isApproxToConstant(0));
    T_w_i_current = pose_linearized.T_w_i;
  }

  inline const Sophus::SE3d& getPose() const {
    if (!linearized) {
      return pose_linearized.T_w_i;
    } else {
      return T_w_i_current;
    }
  }

  inline const Sophus::SE3d& getPoseLin() const {
    return pose_linearized.T_w_i;
  }

  inline void applyInc(const VecN& inc) {
    if (!linearized) {
      PoseState::incPose(inc, pose_linearized.T_w_i);
    } else {
      delta += inc;
      T_w_i_current = pose_linearized.T_w_i;
      PoseState::incPose(delta, T_w_i_current);
    }
  }

  inline bool isLinearized() const { return linearized; }
  inline const VecN& getDelta() const { return delta; }
  inline int64_t getT_ns() const { return pose_linearized.t_ns; }

  inline void backup() {
    backup_delta = delta;
    backup_pose_linearized = pose_linearized;
    backup_T_w_i_current = T_w_i_current;
  }

  inline void restore() {
    delta = backup_delta;
    pose_linearized = backup_pose_linearized;
    T_w_i_current = backup_T_w_i_current;
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 private:
  bool linearized;
  VecN delta;
  PoseState pose_linearized;
  Sophus::SE3d T_w_i_current;

  VecN backup_delta;
  PoseState backup_pose_linearized;
  Sophus::SE3d backup_T_w_i_current;

};

}  // namespace basalt

