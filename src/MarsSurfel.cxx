/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2016-2021, Computer Science Institute VI, University of Bonn.

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

#include "MarsSurfel.h"
#include <Eigen/Eigenvalues>

MarsSurfelBase::MarsSurfelBase()
{
    clear();
}

MarsSurfelBase::MarsSurfelBase(const MarsSurfelBase& surfel)
    : first_view_dir_(surfel.first_view_dir_)
    , num_points_(surfel.num_points_)
    , mean_(surfel.mean_)
    , sum_(surfel.sum_)
    , normal_(surfel.normal_)
    , cov_(surfel.cov_)
    , sum_squares_(surfel.sum_squares_)
    , up_to_date_(surfel.up_to_date_)
    , valid_(surfel.valid_)
    , evaluated_(surfel.evaluated_)
    , c1_ ( surfel.c1_ )
    , eigen_values_ ( surfel.eigen_values_ )
{
}
typename MarsSurfelBase::Ptr MarsSurfelBase::create( MarsSurfelBase::Ptr other )
{
    return ( other != nullptr ) ? std::make_shared<MarsSurfelBase>(*other) : std::make_shared<MarsSurfelBase>();
}

void MarsSurfelBase::clear()
{
    num_points_ = 0.0;
    mean_.setZero();
    cov_.setZero();
    sum_.setZero();
    sum_squares_.setZero();
    first_view_dir_.setZero();
    normal_.setZero();
    eigen_values_.setZero();
    c1_ = 0;

    valid_ = false;
    up_to_date_ = false;
    evaluated_ = false;
}

MarsSurfelBase& MarsSurfelBase::operator+=(const MarsSurfelBase& rhs)
{
    if (rhs.num_points_ > 0 && num_points_ < MAX_SURFEL_POINTS)
    {
        // numerically stable one-pass update scheme
        if (num_points_ < std::numeric_limits<Scalar>::epsilon())
        {
            sum_squares_ = rhs.sum_squares_;
            sum_ = rhs.sum_;
            num_points_ = rhs.num_points_;

            first_view_dir_ = rhs.first_view_dir_;
        }
        else
        {
            const Vec3 deltaS = rhs.num_points_ * sum_ - num_points_ * rhs.sum_;
            sum_squares_ += rhs.sum_squares_ +
                    1.0 / (num_points_ * rhs.num_points_ * (rhs.num_points_ + num_points_)) * deltaS * deltaS.transpose();
            sum_ += rhs.sum_;
            num_points_ += rhs.num_points_;
        }
        up_to_date_ = false;
    }
    return *this;
}

void MarsSurfelBase::add(const Vec3& rhs_sum, const Mat3& rhs_sum_squares, const Scalar & rhs_num_points)
{
    if (rhs_num_points > 0 && num_points_ < MAX_SURFEL_POINTS)
    {
        // numerically stable one-pass update scheme
        if (num_points_ < std::numeric_limits<Scalar>::epsilon())
        {
            sum_squares_ = rhs_sum_squares;
            sum_ = rhs_sum;
            num_points_ = rhs_num_points;
        }
        else
        {
            const Vec3 deltaS = rhs_num_points * sum_ - num_points_ * rhs_sum;
            sum_squares_ += rhs_sum_squares +
                    1.0 / (num_points_ * rhs_num_points * (rhs_num_points + num_points_)) * deltaS * deltaS.transpose();
            sum_ += rhs_sum;
            num_points_ += rhs_num_points;
        }
        up_to_date_ = false;
    }
}

void MarsSurfelBase::add(const Vec3& point)
{
    // numerically stable one-pass update scheme
    if (num_points_ < std::numeric_limits<Scalar>::epsilon())
    {
        sum_ += point;
        num_points_ += 1.0;
        up_to_date_ = false;
        valid_ = false;
    }
    else if (num_points_ < MAX_SURFEL_POINTS)
    {
        const Vec3 deltaS = (sum_ - num_points_ * point);
        sum_squares_ += 1.0 / (num_points_ * (num_points_ + 1.0)) * deltaS * deltaS.transpose();
        sum_ += point;
        num_points_ += 1.0;
        up_to_date_ = false;
    }
    else
    {
        //ROS_WARN_STREAM("numPoints > MAX_SURFEL_POINTS: " << num_points_ << " sum: " << sum_.transpose() << " sum_squares=[" << sum_squares_.row(0) << "; " << sum_squares_.row(1) << "; " << sum_squares_.row(2) << "] pt=[" << point.transpose() << "] " );
        //ROS_WARN_STREAM_THROTTLE(1.0, "num_points_ > MAX_SURFEL_POINTS");
    }
}

void MarsSurfelBase::evaluate()
{
    valid_ = false;
    if (num_points_ > (1.-std::numeric_limits<Scalar>::epsilon()))
    {
        mean_ = sum_ / num_points_;
    }

    if (num_points_ >= MIN_SURFEL_POINTS)
    {
        cov_ = (sum_squares_ / (num_points_ - 1.0)).template cast<Scalar>();

        // enforce symmetry..
        for ( int r = 0; r < cov_.rows(); ++r )
            for ( int c = r; c < cov_.cols(); ++c)
                if ( r != c )
                    cov_(r,c) = cov_(c,r);
    }

    up_to_date_ = true;
    evaluated_ = true;

    // not enough points or surfel degenerate
    if ( num_points_ >= MIN_SURFEL_POINTS )
    {
        // eigen vectors are stored in the columns
        Eigen::SelfAdjointEigenSolver<Mat3> eig(cov_);
        eigen_vectors_ = eig.eigenvectors();
        eigen_values_ = eig.eigenvalues().normalized();

        if ( (eigen_values_.array() > 1e-6).template cast<int>().sum() >= 2 )
            valid_ = true;

        normal_ = eigen_vectors_.col(0).normalized().template cast<Scalar>();
        if (normal_.dot(first_view_dir_) > 0.0)
            normal_ *= -1.0;
    }
}

void MarsSurfelBase::unevaluate()
{
    evaluated_ = false;
    valid_ = false;
}
