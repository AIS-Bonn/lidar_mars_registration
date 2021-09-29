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
#include <vector>
#include <deque>
#include <Eigen/Dense>
#include "MarsFwdDec.h"
#include "EigenFwdDec.h"

struct MarsMapWindow
{
    typedef std::shared_ptr<MarsMapWindow> Ptr;
    static Ptr create( const size_t & cur_queue_size = 2, const bool & store_map_per_cloud = false, const bool & use_adaptive = true );
    std::deque<MarsMapTypePtr> maps;
    std::deque<MarsMapTypePtr> oldMaps;
    std::deque<MarsMapPointCloud::Ptr> clouds;
    std::deque<Sophus::SE3d> clouds_pose;

    size_t queue_size = 2;
    uint64_t last_ts = 0;
    bool m_empty = true;

    bool store_map_per_cloud = false;
    bool m_use_adaptive = true;
    int m_min_num_cells_for_moving_map = 1;

    bool local_map_moved = false;
    Sophus::SE3d local_map_pose ;
    Sophus::SE3d local_map_last_moving_pose ;

    MarsMapTypePtr asOneMap = nullptr;
    MarsMapTypePtr asOneTmpMap = nullptr;

    void translateMap( const Eigen::Vector3d & t );

    MarsMapWindow( const size_t & cur_queue_size = 2, const bool & store_map_per_cloud = false, const bool & use_adaptive = true );

    Sophus::SE3d getClosestKeyFramePose( const Sophus::SE3d & scene_pose ) const;

    MarsMapTypePtr getAsOneCenteredMap();
    std::vector<MarsMapType *> getMapPtrVec();
    std::vector<MarsMapType *> getTransformedMapPtrVec( const Sophus::SE3d & pose  );

    std::vector<MarsMapPointCloud*> getCloudPtrVec();
    std::vector<MarsMapPointCloud::Ptr> getCloudSharedPtrVec();
    std::vector<Sophus::SE3d> getCloudPosesVec() const;
    Eigen::VectorXd getTimeFactors() const;
    Eigen::VectorXt getTimesSince ( const uint64_t & last_ts ) const;
    double getLastTime() const;
    std::vector<uint64_t> getTimes( ) const;
    void addCloud ( MarsMapPointCloud::Ptr newCloud, const Sophus::SE3d & map_pose );

    void transform ( const Sophus::SE3d & tf);

    bool moveWindowToNext();
    bool isEmpty() const;
    void showMap(const std::string & prefix, const Sophus::SE3d & tf , const Sophus::SE3d & linear_transform ) const;
    bool hasFullWindow() const;
    bool hasOutOfWindow() const;
    MarsMapPointCloud::Ptr getCurrentOriginCloud() const;
    MarsMapPointCloud::Ptr getOutOfWindowCloud() const;
    Sophus::SE3d getOutOfWindowCloudPose() const;
    MarsMapPointCloud::Ptr m_out_of_window_cloud = nullptr;
    Sophus::SE3d m_out_of_window_cloud_pose;

    Sophus::SE3d m_next_pose ;
};


