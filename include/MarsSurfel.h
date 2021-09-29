/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2021, Computer Science Institute VI, University of Bonn.
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
#include <memory>
#include "EigenFwdDec.h"

class MarsSurfelBase
{
public:
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar,3,1> Vec3;
    typedef Eigen::Matrix<Scalar,3,3> Mat3;
    typedef Eigen::Matrix<double,3,1> Vec3d;
    typedef std::shared_ptr<MarsSurfelBase> Ptr;
    const unsigned int MAX_NUM_SURFELS = 6;
    const Scalar MIN_SURFEL_POINTS = 10.0;
    const Scalar MAX_SURFEL_POINTS = 10000.0;

    MarsSurfelBase();
    MarsSurfelBase(const MarsSurfelBase & surfel);

    ~MarsSurfelBase(){}

    static Ptr create( Ptr other = nullptr );

    void clear();
    MarsSurfelBase& operator+=(const MarsSurfelBase& rhs);

    void add(const Vec3& point);
    void add(const Vec3& rhs_sum, const Mat3& rhs_sum_squares, const Scalar & rhs_num_points);

    // applies rotation and translation separately
    void transform( const Sophus::SE3d & pose )
    {
        up_to_date_ = false;
        valid_ = false;
        first_view_dir_ = pose.so3().template cast<Scalar>() * first_view_dir_;
        sum_squares_ = (pose.so3().template cast<Scalar>().matrix() * sum_squares_ * pose.so3().inverse().template cast<Scalar>().matrix());
        sum_ += num_points_ * pose.translation().cast<Scalar>();
    }

    void translate( const Vec3 & t )
    {
        up_to_date_ = false;
        valid_ = false;
        sum_ += num_points_ * t;
    }

    void evaluate();

    void unevaluate();
    inline Scalar getNumPoints ( ) const
    {
        return num_points_;
    }

    Scalar num_points_;
    Vec3 sum_;
    Mat3 sum_squares_;

    Vec3 mean_;
    Vec3 normal_;
    Vec3 first_view_dir_;

    Scalar c1_ = 0;
    Vec3 eigen_values_;
    Mat3 eigen_vectors_;
    Mat3 cov_;
    bool up_to_date_ = false;
    bool evaluated_ = false;
    bool valid_ = false;
};

typedef MarsSurfelBase Surfel;
