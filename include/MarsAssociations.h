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
#include "EigenFwdDec.h"
#include "MarsSurfelInfo.h"

struct MarsGMMParameters
{
    static constexpr double sigma_size_factor = 0.75;
    static constexpr double prior_prob = 0.25;
    static constexpr double soft_assoc_c1 = 0.9;
    static constexpr double soft_assoc_c2 = 10.0;
    static constexpr double soft_assoc_c3 = 1.0;
};

typedef float ScalarType;

struct ModelAssociationData
{
    typedef ScalarType Scalar;
    Eigen::Matrix<Scalar,3,3> m_model_cov;
    Eigen::Matrix<Scalar,3,1> m_model_mean;
    Eigen::Matrix<Scalar,3,1> m_model_normal;
    Scalar m_view_dir_weight = 0;
    Scalar m_normal_weight = 0;
    size_t m_model_num_points = 0;
};

struct SceneAssociationData
{
    typedef ScalarType Scalar;
    SurfelInfoConstPtr m_cell_scene = nullptr;
    size_t m_model_num_points = 0;
    size_t m_scene_num_points = 0;
    Scalar m_sigma2  = std::numeric_limits<Scalar>::signaling_NaN();
    Eigen::Matrix<Scalar,3,3> m_scene_cov = Eigen::Matrix<Scalar,3,3>::Constant(std::numeric_limits<Scalar>::signaling_NaN());
    Eigen::Matrix<Scalar,3,1> m_scene_mean = Eigen::Matrix<Scalar,3,1>::Constant(std::numeric_limits<Scalar>::signaling_NaN());
    Eigen::Vector2i m_indices;
};

struct ModelAssociationInfoOnly
{
    typedef ScalarType Scalar;
    ModelAssociationInfoOnly( ) {}
    ~ModelAssociationInfoOnly(){}
    Eigen::Matrix<Scalar,3,1> m_diff;
    Eigen::Matrix<Scalar,3,3> m_W;
};
