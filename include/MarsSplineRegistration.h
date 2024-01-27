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

#pragma once

#include <memory>

#include "EigenFwdDec.h"
#include "MarsFwdDec.h"
#include "MarsPointTypes.h"

struct SceneAssociationsVec;
typedef std::shared_ptr<SceneAssociationsVec> SceneAssociationsVecPtr;

class MarsSplineRegistration
{
public:
    typedef std::shared_ptr<MarsSplineRegistration> Ptr;
    static Ptr create( )
    {
        return std::make_shared<MarsSplineRegistration>();
    }

    MarsSplineRegistration();

    ~MarsSplineRegistration() {}

    bool estimate ( MarsMap* model, const std::vector<MarsMap*> & scene, const Eigen::VectorXt & times, MarsSplineTypePtr & spline )
    {
        return estimate ( model, scene, times, spline, m_max_iterations );
    }
    bool estimate ( MarsMap* model, const std::vector<MarsMap*> & scene, const Eigen::VectorXt & times, MarsSplineTypePtr & spline, const int & maxIterations );

    bool getCovariance( const int & scene_num, const Sophus::SE3d & pose_ws, Eigen::Matrix6d & kf_covar ) const;

    void setMaxIterations ( const int & max_iters )
    {
        m_max_iterations = max_iters;
    }

    void optimizeSplineKnotsFromPoses ( MarsSplineTypePtr & spline, const Eigen::VectorXt & times, const std::vector<Sophus::SE3d> & poses ) const;

    void setUseAdaptive ( const bool & adaptive )
    {
        m_use_adaptive = adaptive;
    }

    void setScaConstraint ( const bool & sca_use_constraint, const double & sca_factor )
    {
        m_sca_use_constraint = sca_use_constraint;
        m_sca_use_constraint_factor = sca_factor;
    }

protected:
    // exposed parameters  
    int m_max_iterations = 3;
    bool m_use_adaptive = false;
    bool m_sca_use_constraint = false;
    double m_sca_use_constraint_factor = 0.001;
    int m_model_num_points;
    std::vector<int> m_scene_num_points;
    SceneAssociationsVecPtr m_prev_scene_assocs;
};
