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
#include <stdarg.h>

#ifdef USE_EASY_PBR
#include "easy_pbr/Mesh.h"
namespace easy_pbr{
    class Mesh;
    typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
typedef easy_pbr::Mesh Cloud;
typedef easy_pbr::MeshSharedPtr CloudPtr;
#else
#include "MarsFwdDec.h"
typedef MarsMapPointCloud Cloud;
typedef MarsMapPointCloudPtr CloudPtr;
#endif
#include "sophus/se3.hpp"
#include <ros/time.h>
#include <geometry_msgs/TransformStamped.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Vector3Stamped.h>

//forward declarations
namespace bag_reader{
    class BagReader;
}
class BagAdaptor;

class BagAdaptor : public std::enable_shared_from_this<BagAdaptor>
{
public:
    template <class ...Args>
    static std::shared_ptr<BagAdaptor> create( Args&& ...args ){
        return std::shared_ptr<BagAdaptor>( new BagAdaptor(std::forward<Args>(args)...) );
    }
    
    CloudPtr get_next_cloud( const bool & shouldSkipOne = false );
    CloudPtr get_next_mav() const;
    std::vector<sensor_msgs::Imu> get_all_imu() const;
    std::vector<sensor_msgs::Imu> get_imu_msgs( const ros::Time & first_ts, const ros::Time & last_ts ) const;
    std::vector<nav_msgs::Odometry> get_all_gps() const;
    std::vector<nav_msgs::Odometry> get_gps_msgs( const ros::Time & first_ts, const ros::Time & last_ts ) const;

    sensor_msgs::PointCloud2::ConstPtr get_cloud_by_stamp ( const ros::Time & bag_stamp );
#ifndef USE_EASY_PBR
    sensor_msgs::Imu::ConstPtr get_next_imu();
    nav_msgs::Odometry::ConstPtr get_next_gps();
#else
    std::shared_ptr<sensor_msgs::Imu> get_next_imu();
    std::shared_ptr<nav_msgs::Odometry> get_next_gps();
#endif

    bool has_imu_until_next_cloud() const;
    bool has_gps_until_next_cloud() const;
    bool has_next_cloud() const;
    void restart();
    int get_current_index() const;
    int get_total_num() const;

    Sophus::SO3d get_lidar_imu_orientation() const { return m_pose_lidar_imu.so3(); }

    Sophus::SE3d get_pose_lidar_imu() const { return m_pose_lidar_imu; }

private:
    BagAdaptor(const std::string config_file); // we put the constructor as private so as to dissalow creating this object on the stack because we want to only used shared ptr for it
    void init_params(const std::string config_file);
    void initializeBagReader();

    CloudPtr get_cloud_by_stamp ( const ros::Time & bag_stamp, const geometry_msgs::TransformStamped & poseMsg );
    CloudPtr cloudMsgToCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, const geometry_msgs::TransformStamped & poseMsg );

    int m_imu_index = 0; int64_t m_last_imu_ts = 0;
    int m_gps_index = 0; int64_t m_last_gps_ts = 0;
    int m_cloud_index = 0; int64_t m_last_cloud_ts = 0;
    std::vector<CloudPtr> m_clouds;
    std::vector<ros::Time> m_bag_stamps;
    std::vector<geometry_msgs::TransformStamped> m_poses;
    std::vector<sensor_msgs::Imu> m_imu_msgs;
    std::vector<nav_msgs::Odometry> m_gps_msgs;
    Sophus::SE3d m_pose_lidar_imu;

    CloudPtr m_mav_mesh;

    bool m_load_imu = false;
    bool m_load_gps = false;
    bool m_compensate_orientation = false;
    double m_compensate_offset = 0.;
    bool m_store_times = false;
    bool m_inc_mode = false;
    bool m_store_compensated = false;
    std::string m_bag_file = "";
    std::string m_pose_file = "";
    std::string m_compensated_bag_file = "";
    bool m_use_pose = false;
    bool m_use_pose_file = false;
    bool m_use_pose_bag_file = false;
    bool m_fix_ouster_tf = false;
    bool m_fix_gps_to_baselink_tf = false;
    std::string m_imu_topic = "";
    std::string m_gps_topic = "";
    std::string m_world_frame = "";
    std::string m_baselink_frame = "";
    std::string m_sensor_frame = "";
    std::string m_imu_frame = "";
    std::string m_cloud_topic = "";
    int m_num_clouds = 0;
    int m_num_zeros = 0;
    int m_cloud_skip_factor = 0;
    std::shared_ptr<bag_reader::BagReader> m_bag_reader;
};
