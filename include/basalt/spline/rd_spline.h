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
@brief Uniform B-spline for euclidean vectors
*/

#pragma once

#include <basalt/spline/spline_common.h>
#include <basalt/utils/assert.h>
#include <basalt/utils/sophus_utils.hpp>

#include <Eigen/Dense>

#include <array>

namespace basalt {

/// @brief Uniform B-spline for euclidean vectors with dimention DIM of order
/// N
///
/// For example, in the particular case scalar values and order N=5, for a time
/// \f$t \in [t_i, t_{i+1})\f$ the value of \f$p(t)\f$ depends only on 5 control
/// points at \f$[t_i, t_{i+1}, t_{i+2}, t_{i+3}, t_{i+4}]\f$. To
/// simplify calculations we transform time to uniform representation \f$s(t) =
/// (t - t_0)/\Delta t \f$, such that control points transform into \f$ s_i \in
/// [0,..,N] \f$. We define function \f$ u(t) = s(t)-s_i \f$ to be a time since
/// the start of the segment. Following the matrix representation of De Boor -
/// Cox formula, the value of the function can be
/// evaluated as follows: \f{align}{
///    p(u(t)) &=
///    \begin{pmatrix} p_{i}\\ p_{i+1}\\ p_{i+2}\\ p_{i+3}\\ p_{i+4}
///    \end{pmatrix}^T M_5 \begin{pmatrix} 1 \\ u \\ u^2 \\ u^3 \\ u^4
///    \end{pmatrix},
/// \f}
/// where \f$ p_{i} \f$ are knots and  \f$ M_5 \f$ is a blending matrix computed
/// using \ref computeBlendingMatrix \f{align}{
///    M_5 = \frac{1}{4!}
///    \begin{pmatrix} 1 & -4 & 6 & -4 & 1 \\ 11 & -12  & -6 & 12  & -4 \\11 &
///    12 &  -6 &  -12  &  6 \\ 1  &  4  &  6  &  4  & -4 \\ 0  &  0  &  0  &  0
///    &  1 \end{pmatrix}.
/// \f}
/// Given this formula, we can evaluate derivatives with respect to time
/// (velocity, acceleration) in the following way:
/// \f{align}{
///    p'(u(t)) &= \frac{1}{\Delta t}
///    \begin{pmatrix} p_{i}\\ p_{i+1}\\ p_{i+2}\\ p_{i+3}\\ p_{i+4}
///    \end{pmatrix}^T
///    M_5
///    \begin{pmatrix} 0 \\ 1 \\ 2u \\ 3u^2 \\ 4u^3 \end{pmatrix},
/// \f}
/// \f{align}{
///    p''(u(t)) &= \frac{1}{\Delta t^2}
///    \begin{pmatrix} p_{i}\\ p_{i+1}\\ p_{i+2}\\ p_{i+3}\\ p_{i+4}
///    \end{pmatrix}^T
///    M_5
///    \begin{pmatrix} 0 \\ 0 \\ 2 \\ 6u \\ 12u^2 \end{pmatrix}.
/// \f}
/// Higher time derivatives are evaluated similarly. This class supports
/// vector values for knots \f$ p_{i} \f$. The corresponding derivative vector
/// on the right is computed using \ref baseCoeffsWithTime.
template <int _DIM, int _N, typename _Scalar = double>
class RdSpline {
 public:
  static constexpr int N = _N;        ///< Order of the spline.
  static constexpr int DEG = _N - 1;  ///< Degree of the spline.

  static constexpr int DIM = _DIM;  ///< Dimension of euclidean vector space.

  static constexpr _Scalar ns_to_s = 1e-9;  ///< Nanosecond to second conversion
  static constexpr _Scalar s_to_ns = 1e9;   ///< Second to nanosecond conversion

  using MatN = Eigen::Matrix<_Scalar, _N, _N>;
  using VecN = Eigen::Matrix<_Scalar, _N, 1>;

  using VecD = Eigen::Matrix<_Scalar, _DIM, 1>;
  using MatD = Eigen::Matrix<_Scalar, _DIM, _DIM>;

  /// @brief Struct to store the Jacobian of the spline
  ///
  /// Since B-spline of order N has local support (only N knots infuence the
  /// value) the Jacobian is zero for all knots except maximum N for value and
  /// all derivatives.
  struct JacobianStruct {
    size_t
        start_idx;  ///< Start index of the non-zero elements of the Jacobian.
    std::array<_Scalar, N> d_val_d_knot;  ///< Value of nonzero Jacobians.
  };

  /// @brief Default constructor
  RdSpline() : dt_ns(0), start_t_ns(0) {}

  /// @brief Constructor with knot interval and start time
  ///
  /// @param[in] time_interval_ns knot time interval in nanoseconds
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  RdSpline(int64_t time_interval_ns, int64_t start_time_ns = 0)
      : dt_ns(time_interval_ns), start_t_ns(start_time_ns) {
    pow_inv_dt[0] = 1.0;
    pow_inv_dt[1] = s_to_ns / dt_ns;

    for (size_t i = 2; i < N; i++) {
      pow_inv_dt[i] = pow_inv_dt[i - 1] * pow_inv_dt[1];
    }
  }

  /// @brief Cast to different scalar type
  template <typename Scalar2>
  inline RdSpline<_DIM, _N, Scalar2> cast() const {
    RdSpline<_DIM, _N, Scalar2> res;

    res.dt_ns = dt_ns;
    res.start_t_ns = start_t_ns;

    for (int i = 0; i < _N; i++) res.pow_inv_dt[i] = pow_inv_dt[i];

    for (const auto k : knots)
      res.knots.emplace_back(k.template cast<Scalar2>());

    return res;
  }

  /// @brief Set start time for spline
  ///
  /// @param[in] start_time_ns start time of the spline in nanoseconds
  inline void setStartTimeNs(int64_t start_time_ns) {
    start_t_ns = start_time_ns;
  }

  /// @brief Maximum time represented by spline
  ///
  /// @return maximum time represented by spline in nanoseconds
  int64_t maxTimeNs() const {
    return start_t_ns + (knots.size() - N + 1) * dt_ns - 1;
  }

  /// @brief Minimum time represented by spline
  ///
  /// @return minimum time represented by spline in nanoseconds
  int64_t minTimeNs() const { return start_t_ns; }

  /// @brief Gererate random trajectory
  ///
  /// @param[in] n number of knots to generate
  /// @param[in] static_init if true the first N knots will be the same
  /// resulting in static initial condition
  void genRandomTrajectory(int n, bool static_init = false) {
    if (static_init) {
      VecD rnd = VecD::Random() * 5;

      for (int i = 0; i < N; i++) knots.push_back(rnd);
      for (int i = 0; i < n - N; i++) knots.push_back(VecD::Random() * 5);
    } else {
      for (int i = 0; i < n; i++) knots.push_back(VecD::Random() * 5);
    }
  }

  /// @brief Add knot to the end of the spline
  ///
  /// @param[in] knot knot to add
  inline void knots_push_back(const VecD& knot) { knots.push_back(knot); }

  /// @brief Remove knot from the back of the spline
  inline void knots_pop_back() { knots.pop_back(); }

  /// @brief Return the first knot of the spline
  ///
  /// @return first knot of the spline
  inline const VecD& knots_front() const { return knots.front(); }

  /// @brief Remove first knot of the spline and increase the start time
  inline void knots_pop_front() {
    start_t_ns += dt_ns;
    knots.pop_front();
  }

  /// @brief Resize containter with knots
  ///
  /// @param[in] n number of knots
  inline void resize(size_t n) { knots.resize(n); }

  /// @brief Return reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return reference to the knot
  inline VecD& getKnot(int i) { return knots[i]; }

  /// @brief Return const reference to the knot with index i
  ///
  /// @param i index of the knot
  /// @return const reference to the knot
  inline const VecD& getKnot(int i) const { return knots[i]; }

  /// @brief Return const reference to deque with knots
  ///
  /// @return const reference to deque with knots
  const Eigen::aligned_deque<VecD>& getKnots() const { return knots; }

  /// @brief Return time interval in nanoseconds
  ///
  /// @return time interval in nanoseconds
  int64_t getTimeIntervalNs() const { return dt_ns; }

  /// @brief Evaluate value or derivative of the spline
  ///
  /// @param Derivative derivative to evaluate (0 for value)
  /// @param[in] time_ns time for evaluating of the spline in nanoseconds
  /// @param[out] J if not nullptr, return the Jacobian of the value with
  /// respect to knots
  /// @return value of the spline or derivative. Euclidean vector of dimention
  /// DIM.
  template <int Derivative = 0>
  VecD evaluate(int64_t time_ns, JacobianStruct* J = nullptr) const {
    int64_t st_ns = (time_ns - start_t_ns);

    BASALT_ASSERT_STREAM(st_ns >= 0, "st_ns " << st_ns << " time_ns " << time_ns
                                              << " start_t_ns " << start_t_ns);

    int64_t s = st_ns / dt_ns;
    double u = double(st_ns % dt_ns) / double(dt_ns);

    BASALT_ASSERT_STREAM(s >= 0, "s " << s);
    BASALT_ASSERT_STREAM(size_t(s + N) <= knots.size(), "s " << s << " N " << N
                                                             << " knots.size() "
                                                             << knots.size());

    VecN p;
    baseCoeffsWithTime<Derivative>(p, u);

    VecN coeff = pow_inv_dt[Derivative] * (blending_matrix_ * p);

    // std::cerr << "p " << p.transpose() << std::endl;
    // std::cerr << "coeff " << coeff.transpose() << std::endl;

    VecD res;
    res.setZero();

    for (int i = 0; i < N; i++) {
      res += coeff[i] * knots[s + i];

      if (J) J->d_val_d_knot[i] = coeff[i];
    }

    if (J) J->start_idx = s;

    return res;
  }

  /// @brief Alias for first derivative of spline. See \ref evaluate.
  inline VecD velocity(int64_t time_ns, JacobianStruct* J = nullptr) const {
    return evaluate<1>(time_ns, J);
  }

  /// @brief Alias for second derivative of spline. See \ref evaluate.
  inline VecD acceleration(int64_t time_ns, JacobianStruct* J = nullptr) const {
    return evaluate<2>(time_ns, J);
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

 protected:
  /// @brief Vector of derivatives of time polynomial.
  ///
  /// Computes a derivative of \f$ \begin{bmatrix}1 & t & t^2 & \dots &
  /// t^{N-1}\end{bmatrix} \f$ with repect to time. For example, the first
  /// derivative would be \f$ \begin{bmatrix}0 & 1 & 2 t & \dots & (N-1)
  /// t^{N-2}\end{bmatrix} \f$.
  /// @param Derivative derivative to evaluate
  /// @param[out] res_const vector to store the result
  /// @param[in] t
  template <int Derivative, class Derived>
  static void baseCoeffsWithTime(const Eigen::MatrixBase<Derived>& res_const,
                                 _Scalar t) {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, N);
    Eigen::MatrixBase<Derived>& res =
        const_cast<Eigen::MatrixBase<Derived>&>(res_const);

    res.setZero();

    if (Derivative < N) {
      res[Derivative] = base_coefficients_(Derivative, Derivative);

      _Scalar _t = t;
      for (int j = Derivative + 1; j < N; j++) {
        res[j] = base_coefficients_(Derivative, j) * _t;
        _t = _t * t;
      }
    }
  }

  template <int, int, typename>
  friend class RdSpline;

  static const MatN
      blending_matrix_;  ///< Blending matrix. See \ref computeBlendingMatrix.

  static const MatN base_coefficients_;  ///< Base coefficients matrix.
                                         ///< See \ref computeBaseCoefficients.

  Eigen::aligned_deque<VecD> knots;    ///< Knots
  int64_t dt_ns;                       ///< Knot interval in nanoseconds
  int64_t start_t_ns;                  ///< Start time in nanoseconds
  std::array<_Scalar, _N> pow_inv_dt;  ///< Array with inverse powers of dt
};

template <int _DIM, int _N, typename _Scalar>
const typename RdSpline<_DIM, _N, _Scalar>::MatN
    RdSpline<_DIM, _N, _Scalar>::base_coefficients_ =
        computeBaseCoefficients<_N, _Scalar>();

template <int _DIM, int _N, typename _Scalar>
const typename RdSpline<_DIM, _N, _Scalar>::MatN
    RdSpline<_DIM, _N, _Scalar>::blending_matrix_ =
        computeBlendingMatrix<_N, _Scalar, false>();

}  // namespace basalt
