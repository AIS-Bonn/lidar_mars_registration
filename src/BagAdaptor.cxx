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
#include "BagAdaptor.h"
#include "BagReader.h"
#include "CloudTransform.h"
#include "sophusToTf.h"

#ifdef USE_EASY_PBR
#include "easy_pbr/Mesh.h"
#include "eigen_utils.h"
using namespace easy_pbr;
using namespace radu::utils;
#endif

//loguru
//#define LOGURU_REPLACE_GLOG 1
#include <loguru.hpp>

//configuru
#define CONFIGURU_WITH_EIGEN 1
#define CONFIGURU_IMPLICIT_CONVERSIONS 1
#include <configuru.hpp>
using namespace configuru;

//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <Tracy.hpp>

#define USE_TBB
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

BagAdaptor::BagAdaptor(const std::string config_file)
{
    ros::Time::init();
    init_params(config_file);
    initializeBagReader( );
}

void BagAdaptor::init_params(const std::string config_file){
    //read all the parameters
    std::string config_file_abs;
    if (fs::path(config_file).is_relative()){
        config_file_abs=(fs::path(PROJECT_SOURCE_DIR) / config_file).string();
    }else{
        config_file_abs=config_file;
    }
    Config cfg = configuru::parse_file(config_file_abs, CFG);
    Config dyn_config=cfg["bag_adaptor"];
    m_bag_file = (std::string) dyn_config["bag_file"];
    m_pose_file = (std::string) dyn_config["pose_file"];
    m_use_pose = dyn_config["use_pose"];
    m_use_pose_file = dyn_config["use_pose_file"];
    m_world_frame = (std::string) dyn_config["world_frame"];
    m_baselink_frame = (std::string) dyn_config["baselink_frame"];
    m_sensor_frame = (std::string) dyn_config["sensor_frame"];
    m_imu_frame = (std::string) dyn_config["imu_frame"];
    m_cloud_topic = (std::string) dyn_config["cloud_topic"];
    m_imu_topic = (std::string) dyn_config["imu_topic"];
    m_gps_topic = (std::string) dyn_config["gps_topic"];
    m_num_clouds = dyn_config["num_clouds"];
    m_cloud_skip_factor = dyn_config["cloud_skip_factor"];
    m_fix_ouster_tf = dyn_config["fix_ouster_tf"];
    m_fix_gps_to_baselink_tf = dyn_config["fix_gps_to_baselink_tf"];
    m_inc_mode = dyn_config["inc_mode"];
    m_store_times = dyn_config["store_times"];
    m_load_imu = dyn_config["load_imu"];
    m_load_gps = dyn_config["load_gps"];
    m_compensate_orientation = dyn_config["compensate_orientation"];
    m_compensate_offset = dyn_config["compensate_offset"];
    m_store_compensated = dyn_config["store_compensated"];
    m_compensated_bag_file = (std::string) dyn_config["compensated_bag_file"];

#ifdef USE_EASY_PBR
    m_mav_mesh = Cloud::create();
    if ( dyn_config["use_mav_mesh"] )
        m_mav_mesh->load_from_file( dyn_config["mav_mesh_file"] );
#endif
}

void BagAdaptor::initializeBagReader( )
{
    m_clouds.clear();
    m_bag_stamps.clear();
    m_poses.clear();

    m_bag_reader = std::make_shared<bag_reader::BagReader>( m_bag_file );
    if ( m_use_pose_file ) m_bag_reader->populateTfBufferFromFile( m_pose_file, m_world_frame, m_sensor_frame );
    if ( m_fix_ouster_tf ) m_bag_reader->addStaticTransform ( "os1_sensor", "os_sensor", Sophus::SE3d() );
    if ( m_fix_gps_to_baselink_tf ) m_bag_reader->addStaticTransform ( "base_link_gps", "base_link", Sophus::SE3d() );

    LOG(INFO) << "Reading clouds... fix: " << m_fix_ouster_tf << " " << m_fix_gps_to_baselink_tf << " min: " << m_bag_reader->getMinTimeStamp() << " max: " << m_bag_reader->getMaxTimeStamp();
    std::vector<geometry_msgs::TransformStamped> poses;

    std::vector<ros::Time> ts;
    if ( m_inc_mode )
    {
        std::vector<ros::Time> bag_stamps;
        m_bag_reader->readCloudBagStamps( bag_stamps, m_cloud_topic, m_num_clouds, m_cloud_skip_factor );
        LOG(INFO) << "bag_stamp size: " << bag_stamps.size();
        std::vector<ros::Time> cloud_stamps;
        cloud_stamps = bag_stamps;
        LOG(INFO) << "cloud_stamp size: " << cloud_stamps.size();
        if ( cloud_stamps.empty() ) LOG(FATAL) << "no clouds.";

        int stamp_idx = 0;
        ros::Time prevStamps = *bag_stamps.begin();
        for ( const auto & stamps : bag_stamps )
        {
            LOG(INFO) << "st["<< stamp_idx<<"]: " << stamps.toNSec() << " p: " << (stamps-prevStamps).toSec();
            ++stamp_idx;
            prevStamps = stamps;
        }

        constexpr bool show_accel = false;
        if ( m_compensate_orientation || show_accel || m_load_imu )
        {
            LOG(INFO) << "loading imu.";
            std::vector<double> times;
            std::vector<sensor_msgs::Imu::ConstPtr> imus;
            m_bag_reader->readMessages<sensor_msgs::Imu>(imus, times, m_imu_topic);
            LOG(INFO) << "copying imu.";
            for ( auto imu : imus ) m_imu_msgs.emplace_back( *imu );
            LOG(INFO) << "loaded imu.";
        }
        if ( m_load_gps )
        {
            LOG(INFO) << "loading gps.";
            std::vector<double> times;
            std::vector<nav_msgs::Odometry::ConstPtr> gps;
            m_bag_reader->readMessages<nav_msgs::Odometry>(gps, times, m_gps_topic);
            LOG(INFO) << "copying gps.";
            for ( auto gps : gps ) m_gps_msgs.emplace_back( *gps );
            LOG(INFO) << "loaded gps.";
        }
        if ( show_accel && !m_imu_msgs.empty() )
        {
            std::vector<double> lin_accel_norms_v(m_imu_msgs.size(),0);
            std::vector<double> ang_vels_norms_v(m_imu_msgs.size(),0);
            #pragma omp parallel for
            for ( size_t idx = 0; idx < m_imu_msgs.size(); ++idx )
            {
                const sensor_msgs::Imu & msg = m_imu_msgs[idx];
                const Eigen::Quaterniond R = Eigen::Quaterniond(0,0,0,1).inverse() * Eigen::Quaterniond(msg.orientation.w, msg.orientation.x, msg.orientation.y, msg.orientation.z );
                const Eigen::Vector3d lin_accel_plus_g ( msg.linear_acceleration.x, msg.linear_acceleration.y, msg.linear_acceleration.z );
                const Eigen::Vector3d lin_accel = lin_accel_plus_g - R * Eigen::Vector3d(0,0,9.806);
                const Eigen::Vector3d ang_vel ( msg.angular_velocity.x, msg.angular_velocity.y, msg.angular_velocity.z );
                lin_accel_norms_v[idx] = lin_accel.norm();
                ang_vels_norms_v[idx] = ang_vel.norm();
            }
            const bool divByTwo = m_imu_msgs.size() % 2 == 0;
            const size_t med_idx = m_imu_msgs.size() / 2;
            auto m1 = lin_accel_norms_v.begin() + (divByTwo ? med_idx+1 : med_idx);
            auto m2 = ang_vels_norms_v.begin() + (divByTwo ? med_idx+1 : med_idx);
            std::nth_element(lin_accel_norms_v.begin(), m1, lin_accel_norms_v.end());
            std::nth_element(ang_vels_norms_v.begin(), m2, ang_vels_norms_v.end());

            const double lin_accel_med = divByTwo ? (lin_accel_norms_v[med_idx]+lin_accel_norms_v[med_idx+1])/2 : lin_accel_norms_v[med_idx];
            const double ang_vels_med = divByTwo ? (ang_vels_norms_v[med_idx]+ang_vels_norms_v[med_idx+1])/2 : ang_vels_norms_v[med_idx];
            Eigen::Map<Eigen::VectorXd> lin_accel_norms(lin_accel_norms_v.data(),m_imu_msgs.size(),1);
            Eigen::Map<Eigen::VectorXd> ang_vels_norms(ang_vels_norms_v.data(),m_imu_msgs.size(),1);

            LOG(INFO) << "\nLinAccel: " << m_imu_msgs[0].linear_acceleration << "\nOri: " << m_imu_msgs[0].orientation;
            LOG(INFO) << "MeanAccel: " << lin_accel_norms.mean() << " MedianAccel: " << lin_accel_med << " stdAccel: " << (std::sqrt((lin_accel_norms.array() - lin_accel_norms.mean()).square().sum()/(lin_accel_norms.size()-1))) << " Max: " << lin_accel_norms.maxCoeff();
            LOG(INFO) << "MeanAngVel: " << ang_vels_norms.mean() << " MedianAngVel: " << ang_vels_med << " stdAngVel: " << (std::sqrt((ang_vels_norms.array() - ang_vels_norms.mean()).square().sum()/(ang_vels_norms.size()-1))) << " Max: " << ang_vels_norms.maxCoeff();
        }
        LOG(INFO) << "num imu_msgs: " << m_imu_msgs.size();

        // store times:
        if ( m_store_times )
        {
            std::ofstream file ( "./cloud_stamps.txt" );
            std::vector<sensor_msgs::PointCloud2::ConstPtr> clouds;
            for ( size_t cloud_index = 0; cloud_index < cloud_stamps.size(); ++cloud_index )
            {
                clouds.clear();
                m_bag_reader->readClouds ( clouds, m_cloud_topic, m_num_clouds, m_cloud_skip_factor, cloud_stamps[cloud_index]-ros::Duration(0.001), cloud_stamps[cloud_index]+ros::Duration(0.001) );
                if ( clouds.empty() ) continue;
                sensor_msgs::PointCloud2::ConstPtr cloudMsg = clouds.front();
                file << cloudMsg->header.stamp.toNSec() << " 0 0 0" << " 0 0 0 1" << std::endl;
            }
            file.close();
            LOG(FATAL) <<"done";
        }
        const geometry_msgs::TransformStamped tf = m_bag_reader->getStaticTransform ( m_baselink_frame, m_sensor_frame );
        LOG(INFO) << "tf: "<< tf;
        const Eigen::Quaterniond tfq(tf.transform.rotation.w, tf.transform.rotation.x, tf.transform.rotation.y, tf.transform.rotation.z);
        if ( tfq.norm() > 0.01 )
        {
            m_q_lidar_imu = tfq.normalized();
            m_pose_lidar_imu.so3().setQuaternion(m_q_lidar_imu);
            m_pose_lidar_imu.translation() << tf.transform.translation.x, tf.transform.translation.y, tf.transform.translation.z;
        }
        else
            m_q_lidar_imu.setIdentity();

        if ( m_bag_reader->hasNonStatic() )
        {
            m_bag_reader->readPoses( cloud_stamps, m_poses, m_bag_reader->hasNonStatic() ? m_world_frame : m_baselink_frame, m_sensor_frame, true );
            for ( size_t idx = 0; idx < m_poses.size(); ++idx )
                LOG(INFO) << "st: " << cloud_stamps[idx].toNSec() << " pose["<<idx<<"]: " << m_poses[idx].transform.translation.x << " " << m_poses[idx].transform.translation.y << " " << m_poses[idx].transform.translation.z << " f: " << (m_bag_reader->hasNonStatic() ? m_world_frame : m_baselink_frame) << " " << m_sensor_frame;
        }
        else
        {
            m_poses.resize(cloud_stamps.size(), tf );
        }

        LOG(INFO) << "cloud_stamp size: " << cloud_stamps.size() << " min: " << m_bag_reader->getMinTimeStamp() << " max: " << m_bag_reader->getMaxTimeStamp();
        m_bag_stamps = cloud_stamps;
        ts = cloud_stamps;

        if ( m_store_compensated && m_compensate_orientation )
        {
            int64_t last_time = (ts[0] - (ts[1]-ts[0])).toNSec();
            std::vector<sensor_msgs::PointCloud2::ConstPtr> clouds;
            std::vector<ros::Time> times;
            times.reserve(m_bag_stamps.size());
            std::vector<sensor_msgs::PointCloud2::ConstPtr> compensated_clouds;
            compensated_clouds.reserve(m_bag_stamps.size());
            for ( size_t cloud_index = 0; cloud_index < m_bag_stamps.size(); ++cloud_index )
            {
                clouds.clear();
                m_bag_reader->readClouds ( clouds, m_cloud_topic, m_num_clouds, m_cloud_skip_factor, m_bag_stamps[cloud_index]-ros::Duration(0.001), m_bag_stamps[cloud_index]+ros::Duration(0.001) );
                if ( clouds.empty() ) continue;
                sensor_msgs::PointCloud2::ConstPtr cloudMsg = clouds.front();
                Eigen::VectorXi point_times(1,1);
                CloudPtr cloud = Cloud::create();
#ifdef USE_EASY_PBR
                copyCloud ( cloudMsg, cloud, &point_times );
#else
                copyCloud<Cloud> ( cloudMsg, cloud, &point_times );
#endif
                point_times.array() += m_compensate_offset;
                const int64_t cur_time = cloudMsg->header.stamp.toNSec();
                compensateOrientation( cloud, point_times, cur_time, last_time, m_imu_msgs, m_q_lidar_imu );
                last_time = cur_time;

                if ( m_store_compensated )
                {
                    sensor_msgs::PointCloud2::Ptr new_cloud_msg ( new sensor_msgs::PointCloud2() );
                    new_cloud_msg->header = cloudMsg->header;
                    new_cloud_msg->data = cloudMsg->data;
                    new_cloud_msg->width = cloudMsg->width;
                    new_cloud_msg->height = cloudMsg->height;
                    new_cloud_msg->fields = cloudMsg->fields;
                    new_cloud_msg->is_bigendian = cloudMsg->is_bigendian;
                    new_cloud_msg->point_step = cloudMsg->point_step;
                    new_cloud_msg->row_step = cloudMsg->row_step;
                    new_cloud_msg->is_dense = cloudMsg->is_dense;
                    uint32_t offset_x = 0;
                    uint32_t offset_y = 0;
                    uint32_t offset_z = 0;
                    for( size_t i = 0; i < cloudMsg->fields.size(); ++i )
                    {
                        if ( cloudMsg->fields[i].name=="x" ) offset_x = cloudMsg->fields[i].offset;
                        if ( cloudMsg->fields[i].name=="y" ) offset_y = cloudMsg->fields[i].offset;
                        if ( cloudMsg->fields[i].name=="z" ) offset_z = cloudMsg->fields[i].offset;
                    }
#ifdef USE_EASY_PBR
                    for ( int i = 0; i < cloud->V.rows(); ++i )
                    {
                        const uint32_t point_start = cloudMsg->point_step * i;
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_x]) = cloud->V(i,0);
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_y]) = cloud->V(i,1);
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_z]) = cloud->V(i,2);
                    }
#else
                    for ( int i = 0; i < cloud->size(); ++i )
                    {
                        const uint32_t point_start = cloudMsg->point_step * i;
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_x]) = cloud->m_points.col(i)(0);
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_y]) = cloud->m_points.col(i)(1);
                        *reinterpret_cast<float*>(&new_cloud_msg->data[point_start + offset_z]) = cloud->m_points.col(i)(2);
                    }
#endif
                    compensated_clouds.emplace_back(new_cloud_msg);
                    times.emplace_back(m_bag_stamps[cloud_index]);
                }
            }
            if ( m_store_compensated ) m_bag_reader->writeMessages<sensor_msgs::PointCloud2>( compensated_clouds, times, m_compensated_bag_file, m_cloud_topic+"_compensated", true);
        }
    }
    else
    {
        std::vector<sensor_msgs::PointCloud2::ConstPtr> clouds;
        m_bag_reader->readCloudsWithPoses( clouds, poses, m_cloud_topic, m_bag_reader->hasNonStatic() ? m_world_frame : m_baselink_frame, m_sensor_frame, m_num_clouds, m_cloud_skip_factor );
        m_clouds.reserve(clouds.size());

        LOG(INFO) << "converting " << int(clouds.size()) << " clouds...";
        for ( size_t idx = 0; idx < clouds.size(); ++idx )
        {
            LOG(INFO) << "idx: " << idx;
            CloudPtr cloud = cloudMsgToCloud( clouds[idx], poses[idx] );
            if ( cloud == nullptr ) continue;
            m_clouds.emplace_back(cloud);
#ifdef USE_EASY_PBR
            LOG(INFO) << "cloud->t: " << cloud->t;
            ts.emplace_back(ros::Time().fromNSec(cloud->t));
#else
            ts.emplace_back(cloud->m_stamp);
#endif
        }
    }
    m_cloud_index = 0;
    m_imu_index = 0;

    LOG(INFO) << "BagAdaptor ready to play. Clouds: " << m_clouds.size() << " stamps: " << m_bag_stamps.size() << " poses: " << m_poses.size();
}

CloudPtr BagAdaptor::get_next_mav() const
{
    CloudPtr mav_mesh = Cloud::create();
#ifdef USE_EASY_PBR
    mav_mesh = std::make_shared<Cloud>(m_mav_mesh->clone());
    Eigen::Affine3d pose = Eigen::Affine3d::Identity();
    if ( m_cloud_index == 0 && !m_clouds.empty() )
        pose = m_clouds[0]->cur_pose();
    else
    {
        if ( m_cloud_index > 0 && m_cloud_index-1 < int(m_clouds.size()) )
             pose = m_clouds[m_cloud_index]->cur_pose();
    }
    mav_mesh->transform_vertices_cpu(pose);
    mav_mesh->worldROS2worldGL();
#endif
    return mav_mesh;
}

bool BagAdaptor::has_next_cloud() const
{
    return m_cloud_index >= 0 && m_inc_mode ? m_cloud_index < int(m_bag_stamps.size()) : m_cloud_index < int(m_clouds.size());
}

CloudPtr BagAdaptor::cloudMsgToCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, const geometry_msgs::TransformStamped & poseMsg )
{
    if ( ! cloudMsg ) return nullptr;

    Sophus::SE3d pose = transformToSophus( poseMsg.transform );
    if ( !m_use_pose ) pose = Sophus::SE3d();

    LOG(INFO) << "pose_stamp: " << poseMsg.header.stamp;
    Eigen::VectorXi times(1,1);

    static int64_t last_time = 0;

    CloudPtr cloud = Cloud::create();
    {
        ZoneScopedN("CopyMsgToCloud");
#ifdef USE_EASY_PBR
        copyCloud ( cloudMsg, cloud, ( m_compensate_orientation ? & times : nullptr) );
#else
        copyCloud<Cloud> ( cloudMsg, cloud, ( m_compensate_orientation ? & times : nullptr) );
#endif
        if ( m_compensate_orientation )
        {
            ZoneScopedN("CopyMsgToCloud::compensateOrientation");
            times.array() += m_compensate_offset;
            compensateOrientation( cloud, times, cloudMsg->header.stamp.toNSec(), last_time, m_imu_msgs, m_q_lidar_imu );
        }
    }
#ifdef USE_EASY_PBR
    LOG(INFO) << "stamp: " << cloudMsg->header.stamp;
    cloud->t = cloudMsg->header.stamp.toNSec();
    cloud->m_is_dirty = true;
    last_time = cloud->t;
    cloud->transform_vertices_cpu(Eigen::Affine3d(pose.matrix()));
    cloud->worldROS2worldGL();// from ROS to GL view frame.
#else
    cloud->m_stamp = cloudMsg->header.stamp.toNSec();
    last_time = cloud->m_stamp;
    cloud->m_pose = pose;
#endif
    return cloud;
}

CloudPtr BagAdaptor::get_cloud_by_stamp ( const ros::Time & bag_stamps, const geometry_msgs::TransformStamped & poseMsg )
{
    if ( ! m_bag_reader ) return nullptr;
    std::vector<sensor_msgs::PointCloud2::ConstPtr> clouds;
    m_bag_reader->readClouds ( clouds, m_cloud_topic, m_num_clouds, m_cloud_skip_factor, bag_stamps-ros::Duration(0.001), bag_stamps+ros::Duration(0.001) );
    LOG(INFO) << "clouds: " << clouds.size() << " stamp: " << bag_stamps.toNSec();
    if ( clouds.empty() ) return nullptr;
    sensor_msgs::PointCloud2::ConstPtr closestCloudMsg = nullptr;
    int64_t closestTimeNs = std::numeric_limits<int64_t>::max();
    for ( sensor_msgs::PointCloud2::ConstPtr cloudMsg : clouds )
    {
        const int64_t curDiff = std::abs<int64_t>((cloudMsg->header.stamp-bag_stamps).toNSec());
        if ( closestTimeNs > curDiff && cloudMsg->header.stamp.toNSec() > m_last_cloud_ts )
        {
            closestTimeNs = curDiff;
            closestCloudMsg = cloudMsg;
        }
    }
    if ( closestCloudMsg != nullptr )
    {
        if ( closestCloudMsg->header.stamp.toNSec() > m_last_cloud_ts )
        {
            m_last_cloud_ts = closestCloudMsg->header.stamp.toNSec();
        }
        return cloudMsgToCloud( closestCloudMsg, poseMsg );
    }
    return nullptr;
}

sensor_msgs::PointCloud2::ConstPtr BagAdaptor::get_cloud_by_stamp ( const ros::Time & bag_stamp )
{
    if ( ! m_bag_reader ) { LOG(INFO) << "no bag_reader loaded."; return nullptr; }
    std::vector<sensor_msgs::PointCloud2::ConstPtr> clouds;
    m_bag_reader->readClouds ( clouds, m_cloud_topic, m_num_clouds, m_cloud_skip_factor, bag_stamp-ros::Duration(0.1), bag_stamp+ros::Duration(0.1) );
    LOG(INFO) << "clouds found: " << clouds.size();
    // sort by cloud stamps
    for ( sensor_msgs::PointCloud2::ConstPtr cloudMsg : clouds )
        if ( cloudMsg->header.stamp.toNSec() > m_last_cloud_ts )
        {
            m_last_cloud_ts = cloudMsg->header.stamp.toNSec();
            return cloudMsg;
        }
    return nullptr;
}

CloudPtr BagAdaptor::get_next_cloud( const bool & shouldSkipOne )
{
    LOG(INFO) << "has next cloud: " << has_next_cloud();
    if ( ! has_next_cloud () ) return nullptr;
    CloudPtr cloud = nullptr;
    if ( ! shouldSkipOne )
    {
        if ( m_inc_mode )
            cloud = get_cloud_by_stamp ( m_bag_stamps[m_cloud_index], m_poses[m_cloud_index] );
        else
            cloud = m_clouds[m_cloud_index];
    }
    ++m_cloud_index;
    return cloud;
}


void BagAdaptor::restart()
{
    m_cloud_index = 0;
    m_imu_index = 0;
}

int BagAdaptor::get_current_index() const
{
    return m_cloud_index;
}

int BagAdaptor::get_total_num() const
{
    return m_inc_mode ? m_bag_stamps.size() : m_clouds.size();
}

std::vector<sensor_msgs::Imu> BagAdaptor::get_all_imu() const
{
    return m_imu_msgs;
}

std::vector<sensor_msgs::Imu> BagAdaptor::get_imu_msgs( const ros::Time & first_ts, const ros::Time & last_ts ) const
{
    if ( ! m_bag_reader ) return std::vector<sensor_msgs::Imu>();
    std::vector<double> times;
    std::vector<sensor_msgs::Imu::ConstPtr> imus;
    m_bag_reader->readMessages<sensor_msgs::Imu>(imus, times, m_imu_topic, std::numeric_limits<int>::max(), m_cloud_skip_factor, first_ts-ros::Duration(0.001), last_ts+ros::Duration(0.001));
    std::vector<sensor_msgs::Imu> imu_msgs;
    for ( const auto & data : imus ) if ( data ) imu_msgs.emplace_back(*data);
    //if ( !imu_msgs.empty() ) m_last_imu_ts = imu_msgs.back().header.stamp.toNSec();
    return imu_msgs;
}

#ifndef USE_EASY_PBR
sensor_msgs::Imu::ConstPtr
#else
std::shared_ptr<sensor_msgs::Imu>
#endif
BagAdaptor::get_next_imu()
{
    if ( ! has_imu_until_next_cloud() ) return nullptr;
#ifndef USE_EASY_PBR
    sensor_msgs::Imu::Ptr msg (new sensor_msgs::Imu);
    *msg = (m_imu_msgs[m_imu_index]);
#else
    std::shared_ptr<sensor_msgs::Imu> msg = std::make_shared<sensor_msgs::Imu>(m_imu_msgs[m_imu_index]);
#endif
    if ( msg ) m_last_imu_ts = msg->header.stamp.toNSec();
    ++m_imu_index;
    return msg;
}

bool BagAdaptor::has_imu_until_next_cloud() const
{
    //LOG(INFO) << "imu_msgs : "<< (!m_imu_msgs.empty()) << " clouds: " << !m_clouds.empty()<< " c<#c: " (m_cloud_index < int(m_clouds.size())) << " i<#i: " (m_imu_index < int(m_imu_msgs.size());
    return !m_imu_msgs.empty() && ((m_inc_mode && !m_bag_stamps.empty() && m_cloud_index < int(m_bag_stamps.size())) || (!m_inc_mode && !m_clouds.empty() && m_cloud_index < int(m_clouds.size())))
                 && m_imu_index < int(m_imu_msgs.size()) &&
            ( ( m_inc_mode && m_imu_msgs[m_imu_index].header.stamp < m_bag_stamps[m_cloud_index])
              ||  (! m_inc_mode &&
           #ifdef USE_EASY_PBR
               m_imu_msgs[m_imu_index].header.stamp.toNSec() < m_clouds[m_cloud_index]->t
           #else
               m_imu_msgs[m_imu_index].header.stamp.toNSec() < m_clouds[m_cloud_index]->m_stamp
           #endif
               ) );
}

std::vector<nav_msgs::Odometry> BagAdaptor::get_all_gps() const
{
    return m_gps_msgs;
}

std::vector<nav_msgs::Odometry> BagAdaptor::get_gps_msgs( const ros::Time & first_ts, const ros::Time & last_ts ) const
{
    if ( ! m_bag_reader ) return std::vector<nav_msgs::Odometry>();
    std::vector<double> times;
    std::vector<nav_msgs::Odometry::ConstPtr> gps;
    m_bag_reader->readMessages<nav_msgs::Odometry>(gps, times, m_gps_topic, std::numeric_limits<int>::max(), m_cloud_skip_factor, first_ts-ros::Duration(0.001), last_ts+ros::Duration(0.001));
    std::vector<nav_msgs::Odometry> gps_msgs;
    for ( const auto & data : gps ) if ( data ) gps_msgs.emplace_back(*data);
    return gps_msgs;
}

#ifndef USE_EASY_PBR
nav_msgs::Odometry::ConstPtr
#else
std::shared_ptr<nav_msgs::Odometry>
#endif
BagAdaptor::get_next_gps()
{
    if ( ! has_gps_until_next_cloud() ) return nullptr;
#ifndef USE_EASY_PBR
    nav_msgs::Odometry::Ptr msg (new nav_msgs::Odometry);
    *msg = (m_gps_msgs[m_gps_index]);
#else
    std::shared_ptr<nav_msgs::Odometry> msg = std::make_shared<nav_msgs::Odometry>(m_gps_msgs[m_gps_index]);
#endif
    if ( msg ) m_last_gps_ts = msg->header.stamp.toNSec();
    ++m_gps_index;
    return msg;
}

bool BagAdaptor::has_gps_until_next_cloud() const
{
    //LOG(INFO) << "gps_msgs : "<< (!m_gps_msgs.empty()) << " clouds: " << !m_clouds.empty()<< " c<#c: " (m_cloud_index < int(m_clouds.size())) << " g<#g: " (m_gps_index < int(m_gps_msgs.size());
    return !m_gps_msgs.empty() && ((m_inc_mode && !m_bag_stamps.empty() && m_cloud_index < int(m_bag_stamps.size())) || (!m_inc_mode && !m_clouds.empty() && m_cloud_index < int(m_clouds.size())))
                 && m_gps_index < int(m_gps_msgs.size()) &&
            ( ( m_inc_mode && m_gps_msgs[m_gps_index].header.stamp < m_bag_stamps[m_cloud_index])
              ||  (! m_inc_mode &&
           #ifdef USE_EASY_PBR
               m_gps_msgs[m_gps_index].header.stamp.toNSec() < m_clouds[m_cloud_index]->t
           #else
               m_gps_msgs[m_gps_index].header.stamp.toNSec() < m_clouds[m_cloud_index]->m_stamp
           #endif
               ) );
}
