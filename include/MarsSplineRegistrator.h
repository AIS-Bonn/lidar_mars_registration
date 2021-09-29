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
#include "EigenFwdDec.h"
#include "MarsFwdDec.h"
#include "MarsSurfelInfo.h"
#include <memory>
#include <deque>
#include <mutex>

class MarsSplineRegistration;
using MarsSplineRegistrationPtr = std::shared_ptr<MarsSplineRegistration>;

namespace easy_pbr{
    class LabelMngr;
}
enum class SurfelType : int;

class MarsSplineRegistrator: public std::enable_shared_from_this<MarsSplineRegistrator>
{
public:

    void init_params(const std::string & config_file);

    //https://stackoverflow.com/questions/29881107/creating-objects-only-as-shared-pointers-through-a-base-class-create-method
    template <class ...Args>
    static std::shared_ptr<MarsSplineRegistrator> create( Args&& ...args ){
        return std::shared_ptr<MarsSplineRegistrator>( new MarsSplineRegistrator(std::forward<Args>(args)...) );
    }

    void setCurPose ( const Sophus::SE3d & cur_pose, const bool & cur_pose_valid = false )
    {
        m_cur_pose = cur_pose;
        m_cur_pose_valid = cur_pose_valid;
    }
    void setCurPos ( const Sophus::SE3d & cur_pos, const bool & cur_pos_valid = false  )
    {
        m_cur_pos = cur_pos;
        m_cur_pos_valid = cur_pos_valid;
    }
    Sophus::SE3d m_cur_pose; // gps pose
    Sophus::SE3d m_cur_pos; // vo pos
    bool m_cur_pose_valid = false;
    bool m_cur_pos_valid = false;

    Sophus::SE3d register_cloud ( MarsMapPointCloudPtr scene_cloud );
    MarsMapPointCloud::Ptr getLocalMapSurfelCenters();
    MarsMapPointCloud::Ptr getSceneMapSurfelCenters();
    Eigen::Matrix6d getOdomCovariance() const{ return m_current_pose_cov; }
    Eigen::Vector6d getOdomVelocity() const;
    std::vector<Sophus::SE3d> getSplinePoses() const;
    Sophus::SE3d getInterpolatedSplinePose() const;
    Sophus::SE3d getMapPose() const;
    void setFirstPose ( const Sophus::SE3d & first_pose_world_sensor );
    void setFirstOrientation ( const Eigen::Quaterniond & first_ori_world_sensor );

    bool was_keyframe ( ) const { return m_was_keyframe; }
    int64_t getLastKeyFrameStamp() const { return m_last_key_frame_stamp; }

    ~MarsSplineRegistrator();
private:
    MarsSplineRegistrator(const std::string & config_file);

    bool m_was_keyframe = false;
    int64_t m_last_key_frame_stamp = -1;
    bool m_firstMarsMesh = true;
    uint16_t m_scan_id = 0;
    int m_scene_map_window = 1;
    int m_local_map_window = 10;
    int m_init_spline_window = 10;
    double m_min_range = 1.;
    double m_scan_frequency = 10;
    double m_min_initial_scans = 1.;
    int m_max_iterations = 2;
    bool m_use_closest_kf = false;
    bool m_use_adaptive = true;

    bool m_show_mars_map = false;
    bool m_show_scene_map = false;
    bool m_show_large_errors = false;
    bool m_show_scene_surfels = false;
    bool m_show_mars_surfels = false;
    bool m_show_initial_surfels = false;
    bool m_show_subdiv_surfels = false;
    bool m_show_knots = false;
    bool m_show_scene_knots = false;
    bool m_show_opt_scene = false;

    SurfelType m_surfel_type;

    double m_second_ev_threshold = 1;
    double m_plane_scale_factor = 0.1;

    bool m_keyframe_use_rotation = false;
    double m_min_keyframe_rotation = 10.*M_PI/180.;
    double m_min_keyframe_distance = 1.;
    double m_interpolationTimeOffset = 0.;
    Sophus::SE3d m_interpolated_scene_pose;

    Eigen::Quaterniond m_first_ori_world_sensor = Eigen::Quaterniond::Identity();
    Sophus::SE3d m_local_map_pose;
    Sophus::SE3d m_first_pose_world_sensor;
    std::vector<Sophus::SE3d> m_scene_poses;
    std::vector<Eigen::Vector6d> m_scene_velocity;

    std::shared_ptr<easy_pbr::LabelMngr> m_label_manager = nullptr;

    std::deque<std::pair<int64_t,Sophus::SE3d>> m_old_poses;

    std::deque<SurfelInfoVector> m_scene_surfels;
    Eigen::Matrix6d m_current_pose_cov;

    MarsMapWindowPtr m_localMarsMap = nullptr;
    MarsMapWindowPtr m_sceneMarsMap = nullptr;
    MarsSplineRegistrationPtr m_marsRegistrator = nullptr;
    MarsMapPointCloudPtr m_mars_map_mesh = nullptr;
    MarsMapPointCloudPtr m_mars_local_map_mesh = nullptr;
    MarsMapPointCloudPtr m_mars_scene_map_mesh = nullptr;
    MarsSplineTypePtr m_trajectory_spline = nullptr;

    std::mutex m_localMapMutex;
    std::mutex m_sceneMapMutex;
};


