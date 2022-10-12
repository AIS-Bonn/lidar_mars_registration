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
#include <Eigen/Dense>
#include "MarsFwdDec.h"
#include <deque>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/TransformStamped.h>

#include <thread>
#include <mutex>
#include <condition_variable>

#ifdef USE_EASY_PBR
namespace easy_pbr{
    class Mesh;
    typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
class BagAdaptor;
#endif

typedef std::map<int64_t,sensor_msgs::Imu> ImuMsgs;
typedef std::shared_ptr<ImuMsgs> ImuMsgsPtr;
typedef std::map<int64_t,nav_msgs::Odometry> GpsMsgs;
typedef std::shared_ptr<GpsMsgs> GpsMsgsPtr;

namespace ros
{
    class Publisher;
    class Subscriber;
    class ServiceServer;
}
typedef std::shared_ptr<ros::Publisher> PublisherPtr;
typedef std::shared_ptr<ros::Subscriber> SubscriberPtr;
typedef std::shared_ptr<ros::ServiceServer> ServiceServerPtr;

namespace tf2_ros
{
    class Buffer;
    class TransformListener;
    class TransformBroadcaster;
}
typedef std::shared_ptr<tf2_ros::Buffer> TransformBufferPtr;
typedef std::shared_ptr<tf2_ros::TransformListener> TransformListenerPtr;
typedef std::shared_ptr<tf2_ros::TransformBroadcaster> TransformBroadcasterPtr;

class MarsSplineRegistrator;
typedef std::shared_ptr<MarsSplineRegistrator> MarsSplineRegistratorPtr;

class MarsSplineRegistrationAdaptor: public std::enable_shared_from_this<MarsSplineRegistrationAdaptor>
{
public:

#ifdef USE_EASY_PBR
    typedef easy_pbr::Mesh Cloud;
    typedef easy_pbr::MeshSharedPtr CloudPtr;
#else
    typedef MarsMapPointCloud Cloud;
    typedef Cloud::Ptr CloudPtr;
#endif
    typedef std::shared_ptr<MarsSplineRegistrationAdaptor> Ptr;

    void init_params(const std::string & config_file);

    //https://stackoverflow.com/questions/29881107/creating-objects-only-as-shared-pointers-through-a-base-class-create-method
    template <class ...Args>
    static Ptr create( Args&& ...args ){
        return std::shared_ptr<MarsSplineRegistrationAdaptor>( new MarsSplineRegistrationAdaptor(std::forward<Args>(args)...) );
    }

    void register_cloud_ros ( const sensor_msgs::PointCloud2::ConstPtr& inputCloud, const ImuMsgsPtr & imu_msgs = nullptr, const GpsMsgsPtr & gps_msgs = nullptr );
    void cloud_msgs_ros ( const sensor_msgs::PointCloud2::ConstPtr& inputCloud );
    void imu_msgs_ros ( const sensor_msgs::Imu::ConstPtr & imu_msg );
    void gps_msgs_ros ( const nav_msgs::Odometry::ConstPtr & gps_msg );

    void register_cloud ( CloudPtr mesh, const ImuMsgsPtr & imu_msgs = nullptr, const GpsMsgsPtr & gps_msgs = nullptr );
    void imu_msgs ( const std::shared_ptr<sensor_msgs::Imu> & imu_msg );
    void gps_msgs ( const std::shared_ptr<nav_msgs::Odometry> & gps_msg );

    void cloud_msg ( CloudPtr mesh );

#ifdef USE_EASY_PBR
    void register_module(const std::shared_ptr<BagAdaptor> bagAdaptor);
#endif

    bool was_keyframe ( ) const;
    void trigger_reset_first_pose();

    ~MarsSplineRegistrationAdaptor();
private:
    MarsSplineRegistrationAdaptor(const std::string & config_file);

    template<typename InputPointCloud, typename MarsPointCloud>
    void convertToMapCloud( const std::shared_ptr<InputPointCloud>& mesh,
                            const typename MarsPointCloud::Ptr& outputCloud,
                            const uint16_t & scan_id );

    void showPoses ( const Sophus::SE3d & cur_pose, const int64_t & cur_time );
    void showPoses ( const Sophus::SE3d & cur_pose, const Sophus::SE3d & reg_pose, const int64_t & cur_time );



    ImuMsgsPtr m_imu_msgs = nullptr;
    GpsMsgsPtr m_gps_msgs = nullptr;

    bool m_center_transform = true;
    bool m_use_tf_for_field_baselink = true;
    bool m_use_tf_for_baselink_sensor = true;
    uint16_t m_scan_id = 0;
    double m_min_range = 1.;
    double m_validityTimeOffsetTF = 0.1;
    int m_num_scan_lines = 128;
    int m_max_scan_line = 128;
    int m_num_scan_columns = 1024;
    int m_num_scans_per_second = 20;
    std::string m_corrected_sensor_frame = "";
    std::string m_corrected_baselink_frame = "";
    std::string m_baselink_gps_frame = "";
    std::string m_baselink_frame = "";
    std::string m_sensor_frame = "";
    std::string m_world_frame = "";
    std::string m_map_frame = "";

    bool m_send_map_baselink_tf = false;

    bool m_compensate_orientation = false;
    bool m_use_gyro_directly = false;
    bool m_organized_scans = true;
    double m_fov_up = 45;
    double m_fov_down = -45;
    double m_inv_vert_sensor_res = 1./0.7;
    double m_inv_hori_sensor_res = m_num_scan_columns / 360.;
    int64_t m_scan_time_offset_ns = 0;

    double m_gps_invalid_threshold = 6;
    double m_downsampled_scan_res = 0.1;

    bool m_publish_fine_occ = true;
    bool m_scan_line_from_organized = false;
    bool m_use_box_self_filter = false;
    Eigen::Vector3f m_self_filter_lower_bound = -Eigen::Vector3f::Ones();
    Eigen::Vector3f m_self_filter_upper_bound = Eigen::Vector3f::Ones();
    bool m_use_angular_self_filter = false;
    int m_max_angular_scan_line = 128;
    int m_min_angular_scan_line = 0;
    int m_angular_self_filter_lower_bound = m_num_scan_columns;
    int m_angular_self_filter_upper_bound = m_num_scan_columns;

    bool m_field_map_init = false;
    Sophus::SE3d m_pose_field_map_init;
    Sophus::SE3d m_gps_world_base_init;
    Sophus::SE3d m_cur_pose_baselink_sensor;

#ifdef USE_EASY_PBR
    std::shared_ptr<BagAdaptor> m_bag_adaptor;
#endif

    MarsSplineRegistratorPtr m_marsRegistrator = nullptr;
    PublisherPtr m_map_publisher = nullptr;
    PublisherPtr m_scene_publisher = nullptr;
    PublisherPtr m_scan_publisher = nullptr;
    PublisherPtr m_scan_downsampled_publisher = nullptr;
    PublisherPtr m_odom_publisher = nullptr;
    SubscriberPtr m_imu_subscriber = nullptr;
    SubscriberPtr m_scan_subscriber = nullptr;
    SubscriberPtr m_gps_subscriber = nullptr;
    ServiceServerPtr m_reset_first_pose_trigger_service = nullptr;
    TransformBufferPtr m_tf_buffer = nullptr;
    TransformListenerPtr m_tf_listener = nullptr;
    TransformBroadcasterPtr m_tf_broadcaster = nullptr;


    std::mutex m_visCloudDequeMutex;
    std::condition_variable m_visCloudDequeCv;
    std::deque<MarsMapPointCloudPtr> m_visCloudDeque;
    std::shared_ptr<std::thread> m_vis_thread;
    std::mutex m_visDequeMutex;

    bool m_run_vis_run = true;
    void runVisMap ( );
    void triggerVis ( MarsMapPointCloudPtr sceneCloud );
};

