/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2021, Computer Science Institute VI, University of Bonn.

This file contains adapted code from
https://gitlab.com/libeigen/eigen

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
*/
#pragma once

#include <Eigen/Dense>
#include "sophus/se3.hpp"

namespace Eigen
{
    typedef Matrix<int64_t, Dynamic, 1 > VectorXt;
    typedef Matrix<double, 6, 6 > Matrix6d;
    typedef Matrix<double, 3, 6 > Matrix36d;
    typedef Matrix<double, 3, 4 > Matrix34d;
    typedef Matrix<double, 6, 2 > Matrix62d;
    typedef Matrix<double, 6, 1 > Vector6d;
    typedef Array<double, 6, 2 > Array62d;
    typedef Matrix<float, 6, 6 > Matrix6f;
    typedef Matrix<float, 3, 6 > Matrix36f;
    typedef Matrix<float, 3, 4 > Matrix34f;
    typedef Matrix<float, 6, 2 > Matrix62f;
    typedef Matrix<float, 6, 1 > Vector6f;
    typedef Array<float, 6, 2 > Array62f;
    typedef Matrix<uint16_t,Eigen::Dynamic,1> VectorXu;
}

// https://math.stackexchange.com/questions/233378/inverse-of-a-3-x-3-covariance-matrix-or-any-positive-definite-pd-matrix
template<typename Scalar>
inline Eigen::Matrix<Scalar,3,3> symDirectInverse( const Eigen::Matrix<Scalar,3,3> & m)
{
    Eigen::Matrix<Scalar,3,3> a = Eigen::Matrix<Scalar,3,3>::Zero();
    a(0,0) = m(2,2) * m(1,1) - m(1,2) * m(1,2);
    a(0,1) = m(0,2) * m(1,2) - m(2,2) * m(0,1);

    a(0,2) = m(0,1) * m(1,2) - m(0,2) * m(1,1);

    a(1,1) = m(2,2) * m(0,0) - m(0,2) * m(0,2);
    a(1,2) = m(0,1) * m(0,2) - m(0,0) * m(1,2);

    a(2,2) = m(0,0) * m(1,1) - m(0,1) * m(0,1);

    const double inv_det = 1./(m(0,0) * a(0,0) + m(0,1) * a(0,1) + m(0,2) * a(0,2));

    return (a*inv_det).template selfadjointView<Eigen::Upper>();
}

template<typename Scalar>
inline Eigen::Matrix<Scalar,3,3> sym2Inverse( const Eigen::Matrix<Scalar,3,3> & m)
{
    const Eigen::Array<Scalar,6,2> b = (Eigen::Array<Scalar,6,2>() <<
                               m(2,2), - m(1,2),
                               m(0,2), - m(2,2),
                               m(0,1), - m(0,2),
                               m(2,2), - m(0,2),
                               m(0,1), - m(0,0),
                               m(0,0), - m(0,1)).finished();
    const Eigen::Array<Scalar,6,2> c = (Eigen::Array<Scalar,6,2>() <<
                               m(1,1), m(1,2),
                               m(1,2), m(0,1),
                               m(1,2), m(1,1),
                               m(0,0), m(0,2),
                               m(0,2), m(1,2),
                               m(1,1), m(0,1)).finished();
    const Eigen::Matrix<Scalar,6,1> a_ = (b * c).rowwise().sum();
    const Eigen::Matrix<Scalar,6,1> a2 = a_ / (a_.template head<3>().dot(m.col(0)));
    return (Eigen::Matrix<Scalar,3,3>() << a2(0),a2(1),a2(2),
                                 a2(1),a2(3),a2(4),
                                 a2(2),a2(4),a2(5)).finished();
}

template<typename Scalar>
inline Eigen::Matrix<Scalar,3,3> symInverse( const Eigen::Matrix<Scalar,3,3> & m )
{
    const Eigen::Matrix<Scalar,6,1> a_ = ((Eigen::Array<Scalar,6,2>() <<
                                 m(2,2), - m(1,2),
                                 m(0,2), - m(2,2),
                                 m(0,1), - m(0,2),
                                 m(2,2), - m(0,2),
                                 m(0,1), - m(0,0),
                                 m(0,0), - m(0,1)).finished() * (Eigen::Array<Scalar,6,2>() <<
                                                                 m(1,1), m(1,2),
                                                                 m(1,2), m(0,1),
                                                                 m(1,2), m(1,1),
                                                                 m(0,0), m(0,2),
                                                                 m(0,2), m(1,2),
                                                                 m(1,1), m(0,1)).finished()).rowwise().sum();
    const Eigen::Matrix<Scalar,6,1> a2 = a_ / (a_.template head<3>().dot(m.col(0)));
    return (Eigen::Matrix<Scalar,3,3>() << a2(0),a2(1),a2(2),
                                           a2(1),a2(3),a2(4),
                                           a2(2),a2(4),a2(5)).finished();
}

template<typename Scalar>
inline Eigen::Matrix<Scalar,3,1> symInverseMultiply( const Eigen::Matrix<Scalar,3,3> & m, const Eigen::Matrix<Scalar,3,1> & b)
{
    const Eigen::Matrix<Scalar,6,1> a_ = ((Eigen::Array<Scalar,6,2>() <<
                                 m(2,2), - m(1,2),
                                 m(0,2), - m(2,2),
                                 m(0,1), - m(0,2),
                                 m(2,2), - m(0,2),
                                 m(0,1), - m(0,0),
                                 m(0,0), - m(0,1)).finished() * (Eigen::Array<Scalar,6,2>() <<
                                                                 m(1,1), m(1,2),
                                                                 m(1,2), m(0,1),
                                                                 m(1,2), m(1,1),
                                                                 m(0,0), m(0,2),
                                                                 m(0,2), m(1,2),
                                                                 m(1,1), m(0,1)).finished()).rowwise().sum();
    const Eigen::Matrix<Scalar,6,1> a = a_ / (a_.template head<3>().dot(m.col(0))); // a(0,0), a(0,1), a(0,2), a(1,1), a(1,2), a(2,2)

    const Eigen::Matrix<Scalar,3,1> a2 ( a(1), a(3), a(4) ); // a(1,0), a(1,1), a(1,2)
    const Eigen::Matrix<Scalar,3,1> a3 ( a(2), a(4), a(5) ); // a(2,0), a(2,1), a(2,2)

    return Eigen::Matrix<Scalar,3,1>( b.dot(a.template head<3>()), b.dot(a2), b.dot(a3) );
}

template<typename Scalar>
inline Eigen::Matrix<Scalar,3,1> symInverseMultiplyDiagAdd( const Eigen::Matrix<Scalar,3,3> & m, const Scalar & s, const Eigen::Matrix<Scalar,3,1> & b)
{
    const Eigen::Matrix<Scalar,6,1> a_ = ((Eigen::Array<Scalar,6,2>() <<
                                 m(2,2)+s, - m(1,2),
                                 m(0,2)  , - m(2,2)-s,
                                 m(0,1)  , - m(0,2),
                                 m(2,2)+s, - m(0,2),
                                 m(0,1)  , - m(0,0)-s,
                                 m(0,0)+s, - m(0,1)).finished() * (Eigen::Array<Scalar,6,2>() <<
                                                                 m(1,1)+s, m(1,2),
                                                                 m(1,2)  , m(0,1),
                                                                 m(1,2)  , m(1,1)+s,
                                                                 m(0,0)+s, m(0,2),
                                                                 m(0,2)  , m(1,2),
                                                                 m(1,1)+s, m(0,1)).finished()).rowwise().sum();
    //const Eigen::Matrix<Scalar,6,1> a = a_ / (a_.template head<3>().dot(m.col(0)));
    const Eigen::Matrix<Scalar,6,1> a = a_ / (a_.template head<3>().dot(Eigen::Matrix<Scalar,3,1>(m(0,0)+s,m(1,0),m(2,0))));

    const Eigen::Matrix<Scalar,3,1> a2 ( a(1), a(3), a(4) );
    const Eigen::Matrix<Scalar,3,1> a3 ( a(2), a(4), a(5) );

    return Eigen::Matrix<Scalar,3,1>( b.dot(a.template head<3>()), b.dot(a2), b.dot(a3) );
}

template<typename Scalar>
inline Scalar symDetAdd( const Eigen::Matrix<Scalar,3,3> & m, const Scalar & s )
{
    return (m(0,0)+s)*(m(1,1)+s)*(m(2,2)+s) + 2*m(1,0)*m(2,0)*m(2,1)  - (m(1,1)+s)*m(2,0)*m(2,0) -(m(2,2)+s)*m(1,0)*m(1,0) - (m(0,0)+s)*m(2,1)*m(2,1);
}
template<typename Scalar>
inline Scalar symDet( const Eigen::Matrix<Scalar,3,3> & m )
{
    return m(0,0)*m(1,1)*m(2,2) + 2*m(1,0)*m(2,0)*m(2,1)  - m(1,1)*m(2,0)*m(2,0) -m(2,2)*m(1,0)*m(1,0) - m(0,0)*m(2,1)*m(2,1);
}
