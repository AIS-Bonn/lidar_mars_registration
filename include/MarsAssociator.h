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

#include "MarsSurfelInfo.h"
#include "MarsFwdDec.h"
#include "sophus/se3.hpp"

struct SceneAssociationData;

template<typename ModelAssociationT>
class MarsAssociator
{
public:
    MarsAssociator( MarsMap * model, const SurfelInfoConstPtrVector* scene_cells,
                    std::vector<SceneAssociationData> * scene_associations,
                    std::vector<ModelAssociationT> * model_associations,
                    int * m_num_valid_scene,
                    int * m_num_valid_model,
                    const Sophus::SE3d & transform, const int & neighbors = 1 );


    ~MarsAssociator(){ }

    void for_each();

    Sophus::SE3d m_transform;

    int m_neighbors = 0;
    int * m_num_valid_scene = nullptr;
    int * m_num_valid_model = nullptr;

    MarsMap * m_model = nullptr;
    const SurfelInfoConstPtrVector* m_scene_cells = nullptr;
    std::vector<SceneAssociationData>* m_scene_associations = nullptr;
    std::vector<ModelAssociationT>* m_model_associations = nullptr;

    template <typename Scalar>
    inline Scalar faster_acos( const Scalar & x_ ) const;
    template <typename Scalar>
    inline Eigen::Array<Scalar,2,1> faster_vec_acos( const Eigen::Array<Scalar,2,1> & x_ ) const;

    inline void operator()( const SurfelInfoConstPtr & sceneCell ) const;
};



