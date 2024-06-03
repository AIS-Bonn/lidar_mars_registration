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
#include "MarsSplineRegistrationAdaptor.h"
#include "CloudTransform.h"
#include "sophusToTf.h"
#ifdef USE_EASY_PBR
#include "BagAdaptor.h"
#include "VisUtils.h"
#endif
#include <Eigen/Geometry>

#define LOGURU_REPLACE_GLOG 1
#include "loguru.hpp"
//configuru
#define CONFIGURU_WITH_EIGEN 1
#define CONFIGURU_IMPLICIT_CONVERSIONS 1
#include <configuru.hpp>
using namespace configuru;

//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


#define TRACY_ONLY_LOCALHOST
#include <Tracy.hpp>

#ifdef USE_EASY_PBR
#include "easy_pbr/Mesh.h"
#include "easy_pbr/Scene.h"
using namespace easy_pbr;
using namespace radu::utils;
#endif

#include <deque>
#include "MarsSplineRegistrator.h"
#include "basalt/spline/se3_spline.h"

#include <ros/ros.h>
#include <ros/service.h>
#include <std_msgs/Empty.h>
#include <std_srvs/Trigger.h>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include <tf2_ros/transform_broadcaster.h>
#include <nav_msgs/Odometry.h>


MarsSplineRegistrationAdaptor::MarsSplineRegistrationAdaptor(const std::string & config_file)
{
    init_params(config_file);
}

MarsSplineRegistrationAdaptor::~MarsSplineRegistrationAdaptor()
{
    m_run_vis_run = false;
    if ( m_vis_thread )
        m_vis_thread->join();
}

#ifdef USE_EASY_PBR
void MarsSplineRegistrationAdaptor::register_module(const std::shared_ptr<BagAdaptor> bagAdaptor){
    m_bag_adaptor = bagAdaptor;
}
#endif

struct FirstPoseStatus
{
    typedef std::shared_ptr<FirstPoseStatus> Ptr;
    bool m_field_map_init = false; // will update it with the next scan.

    FirstPoseStatus() {}

    void set_initialized ( )
    {
        m_field_map_init = true;
    }
    bool is_initialized ( ) const
    {
        return m_field_map_init;
    }

    bool trigger_reset_cb( std_srvs::Trigger::Request & req, std_srvs::Trigger::Response & res )
    {
        m_field_map_init = false;
        res.success = true;
        return true;
    }
};
FirstPoseStatus::Ptr g_first_pose_status;

void MarsSplineRegistrationAdaptor::init_params(const std::string & config_file){
    //read all the parameters
    std::string config_file_abs;
    if (fs::path(config_file).is_relative()){
        config_file_abs=(fs::path(PROJECT_SOURCE_DIR) / config_file).string();
    }else{
        config_file_abs=config_file;
    }
    Config cfg = configuru::parse_file(config_file_abs, CFG);
    Config dyn_config=cfg["registration_adaptor"];
    m_min_range = dyn_config["min_range"];
    m_center_transform = dyn_config["center_transform"];

    m_max_scan_line = dyn_config["max_scan_line"];
    m_num_scan_lines = dyn_config["num_scan_lines"];
    m_num_scan_columns = dyn_config["num_scan_columns"];
    m_num_scans_per_second = dyn_config["num_scans_per_second"];

    m_scan_time_offset_ns = TimeConversion::to_ns(double(dyn_config["scan_time_offset"]));
    m_stamp_at_front = dyn_config["stamp_at_front"];
    m_compensate_orientation = dyn_config["compensate_orientation"];
    m_use_gyro_directly = dyn_config["use_gyro_directly"];
    m_organized_scans = dyn_config["organized_scans"];
    m_fov_up = dyn_config["fov_up"];
    m_fov_down = dyn_config["fov_down"];
    m_inv_vert_sensor_res = 1./ ((std::abs(m_fov_up) + std::abs(m_fov_down)) / m_num_scan_lines);
    m_inv_hori_sensor_res = m_num_scan_columns / 360.;

    m_use_box_self_filter = dyn_config["use_box_self_filter"];
    m_self_filter_lower_bound = dyn_config["box_self_filter_lower_bound"];
    m_self_filter_upper_bound = dyn_config["box_self_filter_upper_bound"];

    m_use_angular_self_filter = dyn_config["use_angular_self_filter"];
    m_max_angular_scan_line = dyn_config["max_angular_scan_line"];
    m_min_angular_scan_line = dyn_config["min_angular_scan_line"];
    m_angular_self_filter_lower_bound = dyn_config["angular_self_filter_lower_bound"];
    m_angular_self_filter_upper_bound = dyn_config["angular_self_filter_upper_bound"];

#ifdef USE_EASY_PBR
    if ( m_use_box_self_filter )
    {
        // show box_self_filter:
        VisMesh::Ptr v = VisMesh::create();
        const float lx = m_self_filter_lower_bound.x(), ly = m_self_filter_lower_bound.y(), lz = m_self_filter_lower_bound.z();
        const float px = m_self_filter_upper_bound.x(), py = m_self_filter_upper_bound.y(), pz = m_self_filter_upper_bound.z();
        v->addEdge(Eigen::Vector3f(lx,ly,lz),Eigen::Vector3f(px,ly,lz),VisMesh::red);
        v->addEdge(Eigen::Vector3f(lx,ly,lz),Eigen::Vector3f(lx,py,lz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(lx,ly,lz),Eigen::Vector3f(lx,ly,pz),VisMesh::blue);

        v->addEdge(Eigen::Vector3f(px,py,pz),Eigen::Vector3f(lx,py,pz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(px,py,pz),Eigen::Vector3f(px,ly,pz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(px,py,pz),Eigen::Vector3f(px,py,lz),VisMesh::green);

        v->addEdge(Eigen::Vector3f(lx,py,pz),Eigen::Vector3f(lx,ly,pz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(px,ly,pz),Eigen::Vector3f(lx,ly,pz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(lx,py,lz),Eigen::Vector3f(lx,py,pz),VisMesh::green);

        v->addEdge(Eigen::Vector3f(px,ly,pz),Eigen::Vector3f(px,ly,lz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(lx,py,lz),Eigen::Vector3f(px,py,lz),VisMesh::green);
        v->addEdge(Eigen::Vector3f(px,ly,lz),Eigen::Vector3f(px,py,lz),VisMesh::green);

        v->showEdgeMesh("box_self_filter", true);
    }
#endif


    m_marsRegistrator = MarsSplineRegistrator::create( config_file );

    g_first_pose_status = std::make_shared<FirstPoseStatus>();

    m_map_frame = (std::string)dyn_config["map_frame"];
    const bool use_ros = dyn_config["use_ros"];
    if ( use_ros )
    {
        ros::NodeHandle nh;
        ros::NodeHandle pnh("~");
        m_validityTimeOffsetTF = dyn_config["tf_validity_time_offset"];
        m_use_tf_for_field_baselink = dyn_config["use_tf_for_field_baselink"];
        m_use_tf_for_baselink_sensor = dyn_config["use_tf_for_baselink_sensor"];
        const std::string imu_topic = (std::string)dyn_config["imu_topic"];
        m_imu_subscriber = std::make_shared<ros::Subscriber>(nh.subscribe(imu_topic, 50, &MarsSplineRegistrationAdaptor::imu_msgs_ros, this, ros::TransportHints().tcpNoDelay()));
        const std::string cloud_topic = (std::string)dyn_config["cloud_topic"];
        m_scan_subscriber = std::make_shared<ros::Subscriber>(nh.subscribe(cloud_topic, 1, &MarsSplineRegistrationAdaptor::cloud_msgs_ros, this, ros::TransportHints().tcpNoDelay()));
        LOG(INFO) << "Subscribed for clouds to topic: " << cloud_topic;
        const std::string gps_topic = (std::string)dyn_config["gps_topic"];
        m_gps_subscriber = std::make_shared<ros::Subscriber>(nh.subscribe(gps_topic, 10, &MarsSplineRegistrationAdaptor::gps_msgs_ros, this, ros::TransportHints().tcpNoDelay()));
        const std::string map_topic = (std::string)dyn_config["map_topic"];
        m_map_publisher = std::make_shared<ros::Publisher>( pnh.advertise<sensor_msgs::PointCloud2>(map_topic, 1));
        const std::string scene_topic = (std::string)dyn_config["scene_topic"];
        m_scene_publisher = std::make_shared<ros::Publisher>( pnh.advertise<sensor_msgs::PointCloud2>(scene_topic, 1));
        const std::string scan_topic = (std::string)dyn_config["scan_topic"];
        m_scan_publisher = std::make_shared<ros::Publisher>( pnh.advertise<sensor_msgs::PointCloud2>(scan_topic, 1));
        m_downsampled_scan_res = dyn_config["scan_downsample_resolution"];
        const std::string scan_downsampled_topic = (std::string)dyn_config["scan_downsampled_topic"];
        m_scan_downsampled_publisher = std::make_shared<ros::Publisher>( pnh.advertise<sensor_msgs::PointCloud2>(scan_downsampled_topic, 1));
        const std::string odom_topic = (std::string)dyn_config["odom_topic"];
        m_odom_publisher = std::make_shared<ros::Publisher>( pnh.advertise<nav_msgs::Odometry>(odom_topic, 1));

        m_corrected_sensor_frame = (std::string)dyn_config["corrected_sensor_frame"];
        m_corrected_baselink_frame = (std::string)dyn_config["corrected_baselink_frame"];
        m_baselink_gps_frame = (std::string)dyn_config["baselink_gps_frame"];
        m_baselink_frame = (std::string)dyn_config["baselink_frame"];
        m_sensor_frame = (std::string)dyn_config["sensor_frame"];
        m_world_frame = (std::string)dyn_config["world_frame"];
        m_send_map_baselink_tf = dyn_config["send_map_baselink_tf"];
        m_tf_buffer = std::make_shared<tf2_ros::Buffer>();
        m_tf_listener = std::make_shared<tf2_ros::TransformListener>( *m_tf_buffer );
        m_tf_broadcaster = std::make_shared<tf2_ros::TransformBroadcaster>();
        m_tf_broadcaster->sendTransform(toTransformStamped(sophusToTransform(Sophus::SE3d()), ros::Time::now(), m_world_frame, m_map_frame));

        const std::string reset_first_pose_topic = (std::string)dyn_config["reset_first_pose_topic"];
        m_reset_first_pose_trigger_service = std::make_shared<ros::ServiceServer>(
                    pnh.advertiseService<FirstPoseStatus,std_srvs::Trigger::Request,std_srvs::Trigger::Response>
                    (reset_first_pose_topic, &FirstPoseStatus::trigger_reset_cb, g_first_pose_status.get()));

        {
            m_vis_thread = std::make_shared<std::thread>( &MarsSplineRegistrationAdaptor::runVisMap, this );
            m_run_vis_run = true;
        }
    }
    LOG(INFO) << "Initialized. Waiting for data...";
}

bool MarsSplineRegistrationAdaptor::was_keyframe ( ) const
{
    return m_marsRegistrator->was_keyframe();
}

std::mutex imuMsgsMutex, gpsMsgsMutex;
void MarsSplineRegistrationAdaptor::imu_msgs_ros ( const sensor_msgs::Imu::ConstPtr & imu_msg )
{
    if ( imu_msg == nullptr ) { return; }
    std::unique_lock<std::mutex> lock ( imuMsgsMutex );
    if ( ! m_imu_msgs )
        m_imu_msgs = std::make_shared<ImuMsgs>();
    (*m_imu_msgs)[imu_msg->header.stamp.toNSec()]=*imu_msg;
//    LOG(1) << "imu msgs: " <<m_imu_msgs->size();
}
void MarsSplineRegistrationAdaptor::gps_msgs_ros ( const nav_msgs::Odometry::ConstPtr& gps_msg)
{
    if ( gps_msg == nullptr ) { return; }
    std::unique_lock<std::mutex> lock ( gpsMsgsMutex );
    if ( ! m_gps_msgs )
        m_gps_msgs = std::make_shared<GpsMsgs>();
    //m_gps_msgs->push_back(*gps_msg);
    (*m_gps_msgs)[gps_msg->header.stamp.toNSec()]=*gps_msg;
//    LOG(1) << "gps msgs: " <<m_gps_msgs->size();
}

template <typename Msg, typename MsgsPtr>
MsgsPtr getCurrentMsgs( const MsgsPtr & m_msgs, std::mutex & msgMutex, const int64_t & last_cloud_stamp, const int64_t & current_cloud_stamp )
{
    MsgsPtr msgs = std::make_shared<Msg>();
    {
        std::unique_lock<std::mutex> lock ( msgMutex );
        if ( m_msgs )
        {
            const auto it_beg = m_msgs->lower_bound(last_cloud_stamp);
            const auto it_end = m_msgs->upper_bound(current_cloud_stamp);
            //LOG(1) << "grabbing msgs: " << std::distance(it_beg, it_end)<< " end? " << (it_beg == m_msgs->end()) << " " << (it_end == m_msgs->end()) << " s: " << m_msgs->size();
            if ( it_beg != m_msgs->end() )
            {
                msgs->insert(it_beg,it_end);
                //LOG(1) << "inserted mgs: " << msgs->size() << " db: " << std::distance(m_msgs->begin(),it_beg) << " db: " << std::distance(it_beg,m_msgs->end());
                m_msgs->erase(m_msgs->begin(),it_beg);
                //LOG(1) << "new msgs size: " << m_msgs->size();
            }
        }
    }
    if ( msgs->empty() ) msgs = nullptr;
    return msgs;
}


int64_t get_front_stamp_offset ( const sensor_msgs::PointCloud2::ConstPtr & cloud_msg )
{
    uint32_t offset_time = 0;
    bool has_relative_time = false;
    for(size_t i=0; i<cloud_msg->fields.size(); ++i)
    {
        if  ( cloud_msg->fields[i].name=="t" ) // ouster OS
        {
            has_relative_time = true;
            offset_time = cloud_msg->fields[i].offset;
            break;
        }
    }
    if ( ! has_relative_time ) return 0;
    // iterate backwards to get "oldest" point offset.
    const int numPoints = cloud_msg->height * cloud_msg->width;
    for ( int point_idx = numPoints-1; point_idx >= 0; --point_idx )
    {
        const uint32_t point_start = cloud_msg->point_step * point_idx;
        const uint32_t & t_ns = *reinterpret_cast<const uint32_t*>(&cloud_msg->data[point_start + offset_time]);
        if ( t_ns > 0 ) return t_ns;
    }
    return 0;
}


void MarsSplineRegistrationAdaptor::cloud_msgs_ros ( const sensor_msgs::PointCloud2::ConstPtr& inputCloud )
{
    if ( ! inputCloud ) return;

    const int64_t front_stamp_offset_ns = m_stamp_at_front ? get_front_stamp_offset ( inputCloud ) : 0 ;

    static int64_t last_cloud_stamp = inputCloud->header.stamp.toNSec() + m_scan_time_offset_ns + front_stamp_offset_ns;
    const int64_t current_cloud_stamp = inputCloud->header.stamp.toNSec() + m_scan_time_offset_ns + front_stamp_offset_ns;
    static int64_t last_seq_id = inputCloud->header.seq;
    const int64_t cur_seq_id = inputCloud->header.seq;

    const int64_t dt = TimeConversion::to_ns(1./(m_num_scans_per_second));
    const int64_t last_cloud_stamp_dt = last_cloud_stamp + 1.05*dt;
    //LOG(1) << "lcs: " << last_cloud_stamp << " dt: " << dt << " lcs+dt: " << last_cloud_stamp_dt;
    if ( last_cloud_stamp_dt < current_cloud_stamp ) // for orientation compensation of os_lidars
    {
       last_cloud_stamp = current_cloud_stamp - dt;
    }

    ImuMsgsPtr imu_msgs = getCurrentMsgs<ImuMsgs,ImuMsgsPtr>( m_imu_msgs, imuMsgsMutex, last_cloud_stamp, current_cloud_stamp );
    GpsMsgsPtr gps_msgs = getCurrentMsgs<GpsMsgs,GpsMsgsPtr>( m_gps_msgs, gpsMsgsMutex, last_cloud_stamp, current_cloud_stamp );

    //LOG(1) << "last: " << last_cloud_stamp << " cur: " << current_cloud_stamp << " imu: " << ( imu_msgs ? imu_msgs->size() : 0 );

    register_cloud_ros(inputCloud, imu_msgs, gps_msgs );
    last_cloud_stamp = current_cloud_stamp;
    last_seq_id = cur_seq_id;
}

void MarsSplineRegistrationAdaptor::imu_msgs ( const std::shared_ptr<sensor_msgs::Imu> & imu_msg )
{
    if ( imu_msg == nullptr ) { return; }
    std::unique_lock<std::mutex> lock ( imuMsgsMutex );
    if ( ! m_imu_msgs )
        m_imu_msgs = std::make_shared<ImuMsgs>();
    //m_imu_msgs->emplace_back(*imu_msg);
    (*m_imu_msgs)[imu_msg->header.stamp.toNSec()]=*imu_msg;
}

void MarsSplineRegistrationAdaptor::gps_msgs ( const std::shared_ptr<nav_msgs::Odometry> & gps_msg )
{
    if ( gps_msg == nullptr ) { return; }
    std::unique_lock<std::mutex> lock ( gpsMsgsMutex );
    if ( ! m_gps_msgs )
        m_gps_msgs = std::make_shared<GpsMsgs>();
    //m_gps_msgs->emplace_back(*gps_msg);
    (*m_gps_msgs)[gps_msg->header.stamp.toNSec()]=*gps_msg;
}

void MarsSplineRegistrationAdaptor::cloud_msg ( CloudPtr mesh )
{
    const int64_t front_stamp_offset_ns = (m_stamp_at_front && mesh->T.size() > 0) ? mesh->T.maxCoeff() : 0 ;
    static int64_t last_cloud_stamp = mesh->t + m_scan_time_offset_ns + front_stamp_offset_ns;
    const int64_t current_cloud_stamp = mesh->t + m_scan_time_offset_ns + front_stamp_offset_ns;
    const int64_t dt = TimeConversion::to_ns(1./(m_num_scans_per_second));
    const int64_t last_cloud_stamp_dt = last_cloud_stamp + 1.05*dt;
    //LOG(1) << "lcs: " << last_cloud_stamp << " dt: " << dt << " lcs+dt: " << last_cloud_stamp_dt;
    if ( last_cloud_stamp_dt < current_cloud_stamp ) // for orientation compensation of os_lidars
    {
       last_cloud_stamp = current_cloud_stamp - dt;
    }
    //LOG(1) << "last: " << last_cloud_stamp << " cur: " << current_cloud_stamp << " beg: " << mesh->t << " front_off: " << front_stamp_offset_ns;

    if constexpr ( false )
    {
        std::unique_lock<std::mutex> lock ( imuMsgsMutex );
        if ( m_imu_msgs )
        {
        for ( const auto & a : *m_imu_msgs )
            LOG(1) << a.first << " " << a.second.header.stamp.toNSec() << " " << (a.second.header.stamp.toNSec() < last_cloud_stamp) << " " << (a.second.header.stamp.toNSec() < current_cloud_stamp) ;
        LOG(1) << "imu size: " << m_imu_msgs->size();
        }
    }

    ImuMsgsPtr imu_msgs = getCurrentMsgs<ImuMsgs,ImuMsgsPtr>( m_imu_msgs, imuMsgsMutex, last_cloud_stamp, current_cloud_stamp );
    GpsMsgsPtr gps_msgs = getCurrentMsgs<GpsMsgs,GpsMsgsPtr>( m_gps_msgs, gpsMsgsMutex, last_cloud_stamp, current_cloud_stamp );

    register_cloud (mesh, imu_msgs, gps_msgs );
}

#define evalCloudType(a,b) [&]{ if constexpr( is_mars_cloud_v<InputPointCloud> ) return (a); else return (b);}()
#define evalSemanticType(a,b) [&]{ if constexpr( is_mars_cloud_v<InputPointCloud> ) { if constexpr( has_semantics_v<InputPointCloud> ) { return (a); } else {return 0;} } else return (b);}()

template<typename InputPointCloud,typename MarsPointCloud>
void MarsSplineRegistrationAdaptor::convertToMapCloud( const std::shared_ptr<InputPointCloud>& inputCloud,
                                                       const typename MarsPointCloud::Ptr& outputCloud,
                                                       const uint16_t & scan_id )
{
    if ( ! inputCloud ) LOG(FATAL) << "empty inputCloud?";
    const int num_pts = evalCloudType(inputCloud->size(),inputCloud->V.rows());
    const bool has_time = evalCloudType(inputCloud->m_scan_time.size(),inputCloud->T.size()) > 0;
    const bool has_intensity = evalCloudType(inputCloud->m_intensity.size(),inputCloud->I.rows()) > 0;
    const bool has_reflectivity = evalCloudType(inputCloud->m_reflectivity.size(),inputCloud->C.rows()) > 0;
    const bool has_semantic = evalSemanticType(inputCloud->size(),inputCloud->S_pred.rows()) > 0;

    outputCloud->resize( num_pts,  true );
    Eigen::VectorXi & scan_time = outputCloud->m_scan_time;
    Eigen::VectorXf & intensity = outputCloud->m_intensity;
    Eigen::VectorXu & reflectivity = outputCloud->m_reflectivity;
    Eigen::VectorXu & scan_line_ids = outputCloud->m_scan_line_id;

    for ( int i = 0; i < num_pts; ++i )
    {
        const Eigen::Vector3f v = evalCloudType( inputCloud->m_points.col(i), inputCloud->V.row(i).transpose().template cast<float>() );
        if ( !v.allFinite() ) continue;

        const float range = evalCloudType( v.norm(), inputCloud->D(i,0) );

        const int scan_line_id = m_organized_scans ? i % m_num_scan_lines : ((M_PI/2-acos(v(2)/range)) * 180./M_PI - m_fov_down) * m_inv_vert_sensor_res;

        //if ( scan_line_id > max_scan_line_id ) max_scan_line_id = scan_line_id;
        //if ( scan_line_id < min_scan_line_id ) min_scan_line_id = scan_line_id;
        //if ( scan_column > max_scan_column ) max_scan_column = scan_column;
        //if ( scan_column < min_scan_column ) min_scan_column = scan_column;
        //if ( ! m_organized_scans && i == inputCloud->V.rows()-1 ) LOG(1) << "pt: " << inputCloud->V.row(i) << " row: "  << scan_line_id << " ( " << min_scan_line_id << " - " << max_scan_line_id << ")" << " col: "<< scan_column << " ( " << min_scan_column << " - " << max_scan_column << ")" ;

        if ( scan_line_id >= m_max_scan_line ) continue;
        if ( range < m_min_range ) continue;
        if ( m_use_box_self_filter && (m_self_filter_lower_bound.array() < v.array()).all() && (v.array() < m_self_filter_upper_bound.array()).all() ) continue;

        if ( m_use_angular_self_filter )
        {
            const int scan_column = m_organized_scans ? i / m_num_scan_lines : atan2(v(1),v(0)) * 180./M_PI * m_inv_hori_sensor_res + m_num_scan_columns/2;
            if ( (scan_line_id > m_max_angular_scan_line || scan_line_id < m_min_angular_scan_line) && m_angular_self_filter_lower_bound < scan_column && scan_column < m_angular_self_filter_upper_bound ) continue;
        }

        const int idx = outputCloud->size();
        outputCloud->addPoint(v);

        if ( has_time ) scan_time(idx) = evalCloudType( inputCloud->m_scan_time(i), inputCloud->T.row(i)(0) );
        if ( has_intensity ) intensity(idx) = evalCloudType( inputCloud->m_intensity(i), inputCloud->I.row(i)(0) );
        if ( has_reflectivity ) reflectivity(idx) = evalCloudType( inputCloud->m_reflectivity(i), inputCloud->C.row(i)(0)*65535 );

        if constexpr (has_semantics_v<MarsPointCloud> && ( has_semantics_v<InputPointCloud> || ! is_mars_cloud_v<InputPointCloud> ))
        {
            if ( has_semantic )
            {
                outputCloud->m_semantic.col(idx) = evalCloudType ( inputCloud->m_semantic.col(i), inputCloud->S_pred.row(i).transpose() );
                outputCloud->m_class(idx) = evalCloudType ( inputCloud->m_class(i), 0 );
            }
        }
        scan_line_ids(idx) = scan_line_id;
    }
    outputCloud->m_pose = evalCloudType( inputCloud->m_pose, Sophus::SE3d(inputCloud->cur_pose().matrix()) );
    outputCloud->m_stamp = evalCloudType( inputCloud->m_stamp, inputCloud->t );
    outputCloud->m_scan_id = Eigen::VectorXi::Constant(1,1,scan_id);
}

void MarsSplineRegistrationAdaptor::showPoses ( const Sophus::SE3d & cur_pose, const Sophus::SE3d & reg_pose, const int64_t & cur_time )
{
#ifdef USE_EASY_PBR
    static VisMesh::Ptr gm = nullptr;
    static VisMesh::Ptr rm = nullptr;
    if ( ! gm  )
    {
        gm = VisMesh::create();
        rm = VisMesh::create();
    }
    rm->resetForAggregation();
    rm->addPoint( reg_pose.translation().cast<float>(), VisMesh::blue );
    rm->showAggregatedPointMesh(Sophus::SE3d(),"reg_pose");

    if (cur_pose.translation().norm() < 10e3 )
    {
        gm->resetForAggregation();
        gm->addPoint( cur_pose.translation().cast<float>(), VisMesh::red);
        gm->showAggregatedPointMesh(Sophus::SE3d(),"gps_pose");
    }
    LOG(1) << "GPS: " << cur_pose.params().transpose();
#endif

    //static std::ofstream regPosesFile ("./reg_poses.txt");
    static std::ofstream regPosesFile ("./lidar_mars_registration_after_map_poses.txt");
    static std::ofstream gpsPosesFile ("./gps_poses.txt");

    if ( gpsPosesFile.is_open() )
    {
        Eigen::Vector3d t = cur_pose.translation();
        Eigen::Quaterniond q = cur_pose.unit_quaternion();
        gpsPosesFile << cur_time << " " << t.x() << " " << t.y() << " " << t.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<"\n";
    }
    if ( regPosesFile.is_open() )
    {
        Eigen::Vector3d t = reg_pose.translation();
        Eigen::Quaterniond q = reg_pose.unit_quaternion();
        regPosesFile << cur_time << " " << t.x() << " " << t.y() << " " << t.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<"\n";
    }

}

template <typename MarsPointCloud>
sensor_msgs::PointCloud2Ptr cloud2Msg ( typename MarsPointCloud::Ptr cloud, const std::string & cloud_frame, const int & intensity = 0, const bool & intensity_as_flag = false )
{
    sensor_msgs::PointCloud2Ptr msg ( new sensor_msgs::PointCloud2() );

    // add fields
    msg->fields.resize(5); // x,y,z,rgb,intensity
    msg->fields[0].name = "x";
    msg->fields[1].name = "y";
    msg->fields[2].name = "z";
    msg->fields[3].name = "rgb";
    msg->fields[4].name = "intensity";
    for ( size_t i = 0; i < 5; ++i )
    {
        msg->fields[i].datatype = sensor_msgs::PointField::FLOAT32;
        msg->fields[i].count = 1;
        msg->fields[i].offset = i*sizeof(float);
    }

    const int num_points = cloud->size();
    const size_t point_step = 5*sizeof(float);
    const size_t offset_rgb = msg->fields[3].offset;
    const size_t offset_intensity = msg->fields[4].offset;

    msg->header.frame_id = cloud_frame;
    msg->header.stamp.fromNSec(cloud->m_stamp);
    msg->data.resize(num_points * point_step, 0);
    msg->width = num_points;
    msg->height = 1;

    //msg->is_bigendian = cloud_msg->is_bigendian;
    msg->point_step = point_step;
    msg->row_step = point_step * msg->width;
    msg->is_dense = false;

    struct PointRGB{
      union{
        struct{
          uint8_t b;
          uint8_t g;
          uint8_t r;
          uint8_t a;
        };
        float rgb;
      };
      PointRGB(const Eigen::Vector3f & c) : b(std::min<uint8_t>(255,std::min<uint8_t>(c(2)*255,0))), g(std::min<uint8_t>(255,std::min<uint8_t>(c(1)*255,0))), r(std::min<uint8_t>(255,std::min<uint8_t>(c(0)*255,0))), a(255){}
    };
    const bool has_intensity = cloud->m_intensity.size() > 0;
    bool has_rgb = false;
    if constexpr (has_color_v<MarsPointCloud>) has_rgb = cloud->m_colors.size() == cloud->m_points.size();

    for ( int idx = 0; idx < num_points; ++idx )
    {
        const size_t point_offset = idx * point_step;
        Eigen::Map<Eigen::Vector3f>((float*)&msg->data[point_offset]) = cloud->m_points.col(idx);
        if ( !intensity_as_flag && has_intensity )
            *reinterpret_cast<float*>(&msg->data[point_offset + offset_intensity]) = cloud->m_intensity(idx);
        else
            *reinterpret_cast<float*>(&msg->data[point_offset + offset_intensity]) = intensity;

        if constexpr (has_color_v<MarsPointCloud>)
                if ( has_rgb ) *reinterpret_cast<PointRGB*>(&msg->data[point_offset + offset_rgb]) = PointRGB(cloud->m_colors.col(idx));
    }
    return msg;
}

void publishCloud ( PublisherPtr publ, sensor_msgs::PointCloud2Ptr msg )
{
    publ->publish(msg);
}

template <typename MarsPointCloud>
sensor_msgs::PointCloud2Ptr publishCloud ( PublisherPtr publ, typename MarsPointCloud::Ptr cloud, const std::string & cloud_frame, const int & intensity = 0, const bool & intensity_as_flag = false )
{
    if ( publ->getNumSubscribers() == 0 ) return nullptr;
    sensor_msgs::PointCloud2Ptr msg = cloud2Msg<MarsMapPointCloud>( cloud, cloud_frame, intensity, intensity_as_flag );
    publishCloud( publ, msg );
    return msg;
}

void publishOdom ( PublisherPtr publ, const uint64_t & stamp, const std::string & header_frame, const std::string & child_frame, const Sophus::SE3d & pose, const Eigen::Vector6d & velocity, const Eigen::Matrix6d & cov )
{
    const Eigen::Quaterniond q = pose.unit_quaternion();
    nav_msgs::Odometry odom;
    odom.header.stamp.fromNSec(stamp);
    odom.header.frame_id = header_frame;
    odom.child_frame_id = child_frame;
    odom.pose.pose.position.x = pose.translation()[0];
    odom.pose.pose.position.y = pose.translation()[1];
    odom.pose.pose.position.z = pose.translation()[2];
    odom.pose.pose.orientation.x = q.x();
    odom.pose.pose.orientation.y = q.y();
    odom.pose.pose.orientation.z = q.z();
    odom.pose.pose.orientation.w = q.w();
    memcpy(odom.pose.covariance.data(),cov.data(),sizeof(double)*cov.size());
    publ->publish(odom);
}

template <typename T>
int getClosestIndex ( const int64_t & cur_ts, const T & msgs )
{
    int firstPos = -1;
    if ( msgs != nullptr && !msgs->empty() )
    {
        int64_t min_diff = std::numeric_limits<int64_t>::max();
        for ( auto msg_it = msgs->begin(); msg_it != msgs->end(); ++msg_it )
        {
            const int64_t msg_ts = msg_it->first;
            const int64_t cur_diff = std::abs<int64_t>(msg_ts-cur_ts);
            if ( cur_diff < min_diff )
            {
                firstPos = std::distance(msgs->begin(),msg_it);
                min_diff = cur_diff;
            }
        }
    }
    return firstPos;
}

void MarsSplineRegistrationAdaptor::register_cloud_ros ( const sensor_msgs::PointCloud2::ConstPtr& inputCloud, const ImuMsgsPtr & imu_msgs, const GpsMsgsPtr & gps_msgs )
{
    if ( ! inputCloud ) return;
    ZoneScopedN("MarsSplineRegistrationAdaptor::Register");

    // tf lookup for pose.
    geometry_msgs::TransformStamped gps_transform_world_baselink;
    geometry_msgs::TransformStamped transform_baselink_sensor;
    if ( ! m_use_tf_for_field_baselink )
        gps_transform_world_baselink.transform = sophusToTransform(Sophus::SE3d());
    if ( ! m_use_tf_for_baselink_sensor )
        transform_baselink_sensor.transform = sophusToTransform(Sophus::SE3d());
    try{
        if ( m_use_tf_for_field_baselink )
            gps_transform_world_baselink = m_tf_buffer->lookupTransform( m_world_frame, m_baselink_gps_frame, ros::Time(0) );
        if ( m_use_tf_for_baselink_sensor )
            transform_baselink_sensor = m_tf_buffer->lookupTransform( m_baselink_frame, m_sensor_frame, ros::Time(0) );
    }
    catch (const tf2::TransformException & ex){
        ROS_ERROR_STREAM_THROTTLE(1,ex.what());
        return;
    }

    Sophus::SE3d cur_pose_gps_world_baselink = transformToSophus(gps_transform_world_baselink.transform);
    if ( m_center_transform )
    {
        cur_pose_gps_world_baselink.translation().setZero();
    }

    m_cur_pose_baselink_sensor = transformToSophus(transform_baselink_sensor.transform); //Sophus::SE3d ( cur_pose_baselink_sensor_eigen.matrix() );
    //LOG(1) << "baselink_sensor " << m_cur_pose_baselink_sensor.params().transpose() << " R:\n" << m_cur_pose_baselink_sensor.so3().matrix();

    static ros::Time last_stamp = ros::Time().fromNSec(inputCloud->header.stamp.toNSec()+m_scan_time_offset_ns); // first time without the front offset
    MarsMapPointCloud::Ptr sceneCloud = MarsMapPointCloud::create();
    ros::Time input_stamp;
    Sophus::SO3d rotation_prior;
    bool rotation_prior_valid = false;
    {
    ZoneScopedN("MarsSplineRegistrationAdaptor::Conversion");
    MarsMapPointCloud::Ptr sceneCloudTmp = MarsMapPointCloud::create();
    copyCloud<MarsMapPointCloud> ( inputCloud, sceneCloudTmp, nullptr, m_stamp_at_front );
    convertToMapCloud<MarsMapPointCloud,MarsMapPointCloud>( sceneCloudTmp, sceneCloud, m_scan_id );
    input_stamp = ros::Time().fromNSec(sceneCloud->m_stamp+m_scan_time_offset_ns); // includes the front offset !
    if ( m_compensate_orientation && imu_msgs && input_stamp > last_stamp )
    {
        rotation_prior = compensateOrientation ( sceneCloud, sceneCloud->m_scan_time, input_stamp.toNSec(), last_stamp.toNSec(), *imu_msgs, m_cur_pose_baselink_sensor.so3(), m_min_range*m_min_range, m_use_gyro_directly, false ); // m_stamp_at_front unnecessary, as it has been copied to match at end
        if ( rotation_prior.params().squaredNorm() > 0.1 ) // check that it is not given as invalid quat.
            rotation_prior_valid = true;
        //LOG(1) << "rot prior: " << rotation_prior.params().transpose();
    }
//    else {
//        LOG(1) << "not compensating: " << m_compensate_orientation << " " << (imu_msgs != nullptr) << " "<< input_stamp.toNSec() << " " << last_stamp.toNSec() << " " << (input_stamp > last_stamp);
//    }
    }

    LOG(1) << "sceneCloud: " << sceneCloud->size() << " id:" << m_scan_id << " gps: " << (gps_msgs==nullptr?0:gps_msgs->size()) << " imu: " <<(imu_msgs==nullptr?0:imu_msgs->size());

    Sophus::SE3d cur_gps_world_base;
    bool gps_valid = false;
    const int gps_idx = getClosestIndex<GpsMsgsPtr>( input_stamp.toNSec(), gps_msgs );
    if ( gps_idx >= 0 )
    {
        auto gps_it = gps_msgs->begin(); std::advance(gps_it,gps_idx);
        const nav_msgs::Odometry & gps_msg = gps_it->second;
        cur_gps_world_base.translation() << gps_msg.pose.pose.position.x,gps_msg.pose.pose.position.y,gps_msg.pose.pose.position.z;
        cur_gps_world_base.setQuaternion(Eigen::Quaterniond(gps_msg.pose.pose.orientation.w,gps_msg.pose.pose.orientation.x,gps_msg.pose.pose.orientation.y,gps_msg.pose.pose.orientation.z));
        gps_valid = gps_msg.pose.covariance[0] < m_gps_invalid_threshold;
    }

    static Sophus::SE3d m_pose_gravity_aligned_map;
    static bool gravity_aligned_set = false;
    if ( gravity_aligned_set )
    {
        const int imu_idx = getClosestIndex<ImuMsgsPtr>( input_stamp.toNSec(), imu_msgs );
        if ( imu_idx >= 0 )
        {
            auto imu_it = imu_msgs->begin(); std::advance(imu_it,imu_idx);
            const sensor_msgs::Imu & imu_msg = imu_it->second;
            m_pose_gravity_aligned_map.setQuaternion(Eigen::Quaterniond(imu_msg.orientation.w,imu_msg.orientation.x,imu_msg.orientation.y,imu_msg.orientation.z).inverse());
            gravity_aligned_set = true;
        }
    }

    if ( !g_first_pose_status->is_initialized() )
    {
        g_first_pose_status->set_initialized();
        m_pose_field_map_init = cur_pose_gps_world_baselink; //Sophus::SE3d( cur_pose_world_baselink_gps.matrix() );
        m_gps_world_base_init = cur_gps_world_base;
        if ( m_center_transform )
        {
            m_pose_field_map_init.translation().setZero();
            m_gps_world_base_init.translation().setZero();
        }
        LOG(1) << "initialized Field to Map: " << m_pose_field_map_init.params().transpose();
    }

    const Sophus::SE3d firstPose_lidar_imu = (m_cur_pose_baselink_sensor*m_gps_world_base_init.inverse());
    const Sophus::SE3d curPose_imu_lidar = (cur_gps_world_base * m_cur_pose_baselink_sensor.inverse());
    const Sophus::SE3d cur_pose_gps_world_sensor = (firstPose_lidar_imu * curPose_imu_lidar);

    // assuming the cloud is centered at Identity()
    m_marsRegistrator->setCurRot ( rotation_prior, rotation_prior_valid );
    m_marsRegistrator->setCurPose ( cur_pose_gps_world_sensor, gps_valid );
    const Sophus::SE3d pose_map_sensor = m_marsRegistrator->register_cloud ( sceneCloud );
    const Sophus::SE3d pose_world_map = m_marsRegistrator->getMapPose();
    const Sophus::SE3d interp_pose_map_sensor = pose_world_map * m_marsRegistrator->getInterpolatedSplinePose();
    const Sophus::SE3d interp_pose_map_baselink = m_cur_pose_baselink_sensor * interp_pose_map_sensor * m_cur_pose_baselink_sensor.inverse(); // change of basis necessary?

#ifdef USE_EASY_PBR
    showPoses( cur_pose_gps_world_sensor, pose_map_sensor, sceneCloud->m_stamp );
#endif
    // send tf
    const ros::Duration validityTimeOffset (m_validityTimeOffsetTF);

    const geometry_msgs::Transform field_map = sophusToTransform( m_pose_field_map_init * m_pose_gravity_aligned_map.inverse() );

    const geometry_msgs::TransformStamped stampedTransformWorldMap = toTransformStamped( field_map, input_stamp+validityTimeOffset, m_world_frame, m_map_frame );
    m_tf_broadcaster->sendTransform(stampedTransformWorldMap);

    //Eigen::Vector6d velocity_sensor = m_marsRegistrator->getOdomVelocity();
    const Eigen::Matrix6d cov = m_marsRegistrator->getOdomCovariance();
    const Eigen::Vector6d velocity_baselink = Eigen::Vector6d::Zero();
    const Sophus::SE3d interp_pose_firstMap_baselink = //m_pose_field_map_init *
            m_pose_gravity_aligned_map *
            interp_pose_map_baselink;
    // TODO: add interpolation offset on the time stamp
    publishOdom ( m_odom_publisher, sceneCloud->m_stamp, m_map_frame, m_corrected_baselink_frame, interp_pose_firstMap_baselink, velocity_baselink, cov );

    if ( m_send_map_baselink_tf )
    {
        const geometry_msgs::Transform map_baselink = sophusToTransform( interp_pose_firstMap_baselink );
        const geometry_msgs::TransformStamped stampedTransformMapBaselink =
        toTransformStamped( map_baselink, input_stamp+validityTimeOffset, m_map_frame, m_corrected_baselink_frame );
        m_tf_broadcaster->sendTransform(stampedTransformMapBaselink);
    }

    {
        static std::ofstream regPosesFile ("./lidar_mars_registration_after_map_poses.txt");
        if ( regPosesFile.is_open() )
        {
            const Eigen::Vector3d t = interp_pose_firstMap_baselink.translation();
            const Eigen::Quaterniond q = interp_pose_firstMap_baselink.unit_quaternion();
            regPosesFile << (inputCloud->header.stamp.toNSec()) << " " << t.x() << " " << t.y() << " " << t.z() << " " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << std::endl;
        }
    }


    triggerVis( sceneCloud );

    ++m_scan_id;
    last_stamp = input_stamp;
}

//#define USE_EASY_PBR
void MarsSplineRegistrationAdaptor::register_cloud ( CloudPtr mesh, const ImuMsgsPtr & imu_msgs, const GpsMsgsPtr & gps_msgs )
{
    ZoneScopedN("MarsSplineRegistrationAdaptor::Register");
#ifdef USE_EASY_PBR
    // back to ros coord system
    mesh->worldGL2worldROS();
#endif

    static Sophus::SE3d first_gps_world_base;
    static Eigen::Quaterniond first_orientation_world_sensor = Eigen::Quaterniond::Identity();
    static Eigen::Quaterniond first_quat = Eigen::Quaterniond::Identity();
    static bool firstSet = false;
    if ( ! firstSet )
    {
        firstSet = true;
#ifdef USE_EASY_PBR
        if ( m_bag_adaptor != nullptr )
            m_cur_pose_baselink_sensor = m_bag_adaptor->get_pose_lidar_imu().inverse();
        const int64_t stamp = mesh->t;
#else
        const int64_t stamp = mesh->m_stamp;
#endif
        const int imu_idx = getClosestIndex<ImuMsgsPtr>( stamp, imu_msgs );
        if ( imu_idx >= 0 )
        {
            auto imu_it = imu_msgs->begin(); std::advance(imu_it,imu_idx);
            const sensor_msgs::Imu & msg = imu_it->second;
            first_orientation_world_sensor = m_cur_pose_baselink_sensor.so3().inverse().unit_quaternion() * //Eigen::Quaterniond(0,0,0,1).inverse() *
                    Eigen::Quaterniond(msg.orientation.w, msg.orientation.x, msg.orientation.y,  msg.orientation.z );
            first_quat = Eigen::Quaterniond(msg.orientation.w, msg.orientation.x, msg.orientation.y, msg.orientation.z );
        }

        const int gps_idx = getClosestIndex<GpsMsgsPtr>( stamp, gps_msgs );
        if ( gps_idx >= 0 )
        {
            auto gps_it = gps_msgs->begin(); std::advance(gps_it,gps_idx);
            const nav_msgs::Odometry & gps_msg = gps_it->second;
            first_gps_world_base.translation() << gps_msg.pose.pose.position.x,gps_msg.pose.pose.position.y,gps_msg.pose.pose.position.z;
            first_gps_world_base.setQuaternion(Eigen::Quaterniond(gps_msg.pose.pose.orientation.w,gps_msg.pose.pose.orientation.x,gps_msg.pose.pose.orientation.y,gps_msg.pose.pose.orientation.z));
            const bool invalid_pos = first_gps_world_base.translation().norm() > 10e4;
            if ( gps_msg.pose.covariance[0] < m_gps_invalid_threshold || invalid_pos )
            {
                first_gps_world_base.translation().setZero();
                LOG(1) << "reset first pose.";
            }
        }

        m_marsRegistrator->setFirstOrientation( first_orientation_world_sensor );
    }
    //m_marsRegistrator->setFirstPose ( first_pose_world_sensor );
#ifdef USE_EASY_PBR
    static Sophus::SE3d first_pose_world_sensor = Sophus::SE3d(mesh->cur_pose().matrix());
    Sophus::SE3d cur_pose_world_sensor = Sophus::SE3d(mesh->cur_pose().matrix());
#else
    static Sophus::SE3d first_pose_world_sensor = first_gps_world_base;
    Sophus::SE3d cur_pose_world_sensor = mesh->m_pose;
#endif

    MarsMapPointCloud::Ptr sceneCloud = MarsMapPointCloud::create();
    convertToMapCloud<Cloud,MarsMapPointCloud>(mesh, sceneCloud, m_scan_id);
    if ( m_stamp_at_front )
    {
        sceneCloud->m_stamp += sceneCloud->m_scan_time.head(sceneCloud->size()).maxCoeff();
    }
    // from here on, m_stamp_at_front is not important anymore!

    Sophus::SO3d rotation_prior;
    bool rotation_prior_valid = false;
    static ros::Time last_stamp = ros::Time().fromNSec(sceneCloud->m_stamp+m_scan_time_offset_ns); // first time without the front offset
    const ros::Time input_stamp = ros::Time().fromNSec(sceneCloud->m_stamp+m_scan_time_offset_ns); // includes the front offset !
    if ( m_compensate_orientation && imu_msgs  && input_stamp > last_stamp )
    {
        rotation_prior = compensateOrientation ( sceneCloud, sceneCloud->m_scan_time, input_stamp.toNSec(), last_stamp.toNSec(), *imu_msgs, m_cur_pose_baselink_sensor.so3(), m_min_range*m_min_range, m_use_gyro_directly );
        if ( rotation_prior.params().squaredNorm() > 0.1 ) // check that it is not given as invalid quat.
            rotation_prior_valid = true;
        //LOG(1) << "rot prior: " << rotation_prior.params().transpose();
    }
    //rotation_prior = rotation_prior.inverse();
    //rotation_prior_valid = false;

    LOG(1) << "sceneCloud: " << sceneCloud->size() << " id:" << m_scan_id << " gps: " << (gps_msgs==nullptr?0:gps_msgs->size()) << " imu: " <<(imu_msgs==nullptr?0:imu_msgs->size());

#ifdef USE_EASY_PBR
    sceneCloud->transform_points(cur_pose_world_sensor.inverse()); // center cloud at Identity
#else
    // assuming the cloud is centered at Identity()
#endif

    {
        const int gps_idx = getClosestIndex<GpsMsgsPtr>( sceneCloud->m_stamp, gps_msgs );
        if ( gps_idx >= 0 )
        {
            auto gps_it = gps_msgs->begin(); std::advance(gps_it,gps_idx);
            const nav_msgs::Odometry & gps_msg = gps_it->second;
            Sophus::SE3d raw_gps_pose;
            raw_gps_pose.translation() << gps_msg.pose.pose.position.x,gps_msg.pose.pose.position.y,gps_msg.pose.pose.position.z;
            raw_gps_pose.setQuaternion(Eigen::Quaterniond(gps_msg.pose.pose.orientation.w,gps_msg.pose.pose.orientation.x,gps_msg.pose.pose.orientation.y,gps_msg.pose.pose.orientation.z));

            const Sophus::SE3d pose_lidar_imu = m_cur_pose_baselink_sensor;
            const Sophus::SE3d firstPose_lidar_imu = (pose_lidar_imu*first_gps_world_base.inverse());
            const Sophus::SE3d curPose_imu_lidar = (raw_gps_pose * pose_lidar_imu.inverse());
            const Sophus::SE3d diffPose_lidar = (firstPose_lidar_imu * curPose_imu_lidar);
            cur_pose_world_sensor = diffPose_lidar;
        }
    }
    Sophus::SE3d cur_pose_vo_world_sensor;
    bool cur_pose_vo_valid = false;

    if ( m_center_transform )
        cur_pose_world_sensor.translation().setZero();

    m_marsRegistrator->setCurPose ( cur_pose_world_sensor, false );
    m_marsRegistrator->setCurPos ( cur_pose_vo_world_sensor, cur_pose_vo_valid );
    m_marsRegistrator->setCurRot ( rotation_prior, rotation_prior_valid );
    const Sophus::SE3d pose_map_sensor = m_marsRegistrator->register_cloud ( sceneCloud );
    const Sophus::SE3d pose_world_map = m_marsRegistrator->getMapPose();

#ifdef USE_EASY_PBR
    const Sophus::SE3d interp_pose_map_sensor = pose_world_map * m_marsRegistrator->getInterpolatedSplinePose();
    std::vector<Sophus::SE3d> spline_poses(3,Sophus::SE3d());

    for ( Sophus::SE3d & s : spline_poses )
        s = pose_world_map * s;

    const Sophus::SE3d cur_pose_first_sensor = first_pose_world_sensor.inverse() * cur_pose_world_sensor;
    Eigen::Quaterniond imu_ori = cur_pose_world_sensor.so3().unit_quaternion();

    const int imu_idx = getClosestIndex<ImuMsgsPtr>( mesh->t, imu_msgs );
    if ( imu_idx >= 0 )
    {
        auto imu_it = imu_msgs->begin(); std::advance(imu_it,imu_idx);
        const sensor_msgs::Imu & msg = imu_it->second;
        const Eigen::Quaterniond cur_ori(msg.orientation.w, msg.orientation.x, msg.orientation.y, msg.orientation.z );

        constexpr bool store_ori = false;
        if ( store_ori )
        {
            static std::ofstream imuPosesFile ("./imu_ori.txt");
            const Eigen::Quaterniond q = cur_ori;
            imuPosesFile << msg.header.stamp.toNSec() << " 0.0 0.0 0.0 " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<"\n";
        }
        if ( store_ori )
        {
            static std::ofstream regPosesFile ("./reg_ori.txt");
            const Eigen::Quaterniond q = interp_pose_map_sensor.so3().unit_quaternion();
            regPosesFile << msg.header.stamp.toNSec() << " 0.0 0.0 0.0 " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() <<"\n";
        }
        const Eigen::Quaterniond q_lidar_imu (0.0247437, 0.706674, -0.706674, 0.0247437); // initial estim
        //const Eigen::Quaterniond q_lidar_imu (-0.00355035, 0.701345, -0.712781, 0.00679966); // huber 0.35
        //const Eigen::Quaterniond q_lidar_imu (-0.00277442, 0.70259, -0.711551, 0.00742735); // l2-loss
        //const Eigen::Quaterniond q_lidar_imu = m_cur_pose_baselink_sensor.so3().inverse().unit_quaternion();
        const Eigen::Quaterniond firstOrientation_lidar_imu = (q_lidar_imu*first_quat.inverse()).normalized();
        const Eigen::Quaterniond curOrientation_imu_lidar = (cur_ori * q_lidar_imu.inverse()).normalized();
        const Eigen::Quaterniond diffOrientation_lidar = (firstOrientation_lidar_imu * curOrientation_imu_lidar).normalized();
        imu_ori = diffOrientation_lidar;


        constexpr bool show_with_both_ori = false;
        if ( show_with_both_ori )
        {
            static VisMesh::Ptr om = VisMesh::create();
            static VisMesh::Ptr im = VisMesh::create();
            if ( was_keyframe() )
            {
                Sophus::SE3d imu_pose = interp_pose_map_sensor;
                imu_pose.so3().setQuaternion(imu_ori);
                om->resetForAggregation();
                im->resetForAggregation();
                const Eigen::Matrix3Xf & pts = sceneCloud->toMatrix3Xf();
                const int num_pts = sceneCloud->size();
                om->reservePoints(num_pts);
                im->reservePoints(num_pts);
                for ( int i = 0; i < num_pts; ++i )
                {
                    om->addPoint( interp_pose_map_sensor.cast<float>() * pts.col(i), VisMesh::blue);
                    im->addPoint( imu_pose.cast<float>() * pts.col(i), VisMesh::red);
                }
                om->showAggregatedPointMesh(Sophus::SE3d(),"lidar_cloud");
                im->showAggregatedPointMesh(Sophus::SE3d(),"imu_cloud");
            }
        }
    }
    if ( gps_msgs != nullptr && !gps_msgs->empty() )
    {
        constexpr bool show_with_gps = false;
        if ( show_with_gps )
        {
            const Eigen::Quaterniond q_lidar_imu (0.0247437, 0.706674, -0.706674, 0.0247437); // initial estim
            //const Eigen::Quaterniond q_lidar_imu (-0.00355035, 0.701345, -0.712781, 0.00679966); // huber 0.35
            //const Eigen::Quaterniond q_lidar_imu (-0.00277442, 0.70259, -0.711551, 0.00742735); // l2-loss
            //const Eigen::Quaterniond q_lidar_imu = m_cur_pose_baselink_sensor.so3().inverse().unit_quaternion();

            const Sophus::SE3d pose_lidar_imu = m_cur_pose_baselink_sensor;
            //pose_lidar_imu.setQuaternion( q_lidar_imu );
            const Sophus::SE3d firstPose_lidar_imu = (pose_lidar_imu*first_gps_world_base.inverse());

            static VisMesh::Ptr gm = VisMesh::create();
            static VisMesh::Ptr im = VisMesh::create();
            gm->resetForAggregation();
            im->resetForAggregation();
            constexpr bool show_cloud_with_gps = false;
            if ( show_cloud_with_gps )
            {
                static VisMesh::Ptr im = VisMesh::create();
                {
                    const int gps_idx = getClosestIndex<GpsMsgsPtr>( mesh->t, gps_msgs );
                    if ( gps_idx >= 0 )
                    {
                        auto gps_it = gps_msgs->begin(); std::advance(gps_it,gps_idx);
                        const nav_msgs::Odometry & gps_msg = gps_it->second;
                        Sophus::SE3d raw_gps_pose;
                        raw_gps_pose.translation() << gps_msg.pose.pose.position.x,gps_msg.pose.pose.position.y,gps_msg.pose.pose.position.z;
                        raw_gps_pose.setQuaternion(Eigen::Quaterniond(gps_msg.pose.pose.orientation.w,gps_msg.pose.pose.orientation.x,gps_msg.pose.pose.orientation.y,gps_msg.pose.pose.orientation.z));

                        const Sophus::SE3d curPose_imu_lidar = (raw_gps_pose * pose_lidar_imu.inverse());
                        const Sophus::SE3d diffPose_lidar = (firstPose_lidar_imu * curPose_imu_lidar);
                        const Sophus::SE3d gps_pose = diffPose_lidar;

                        static Sophus::SE3d last_gps_kf_pose = gps_pose;
                        const bool should_be_keyframe = (last_gps_kf_pose.translation() - gps_pose.translation()).norm() > 2;

                        if ( was_keyframe() || should_be_keyframe )
                        {
                            last_gps_kf_pose = gps_pose;
                            Sophus::SE3d imu_pose = gps_pose;
                            im->resetForAggregation();
                            const Eigen::Matrix3Xf & pts = sceneCloud->toMatrix3Xf();
                            const int num_pts = sceneCloud->size();
                            im->reservePoints(num_pts);
                            for ( int i = 0; i < num_pts; ++i )
                            {
                                im->addPoint( imu_pose.cast<float>() * pts.col(i), VisMesh::red);
                            }
                            im->showAggregatedPointMesh(Sophus::SE3d(),"agg_cloud");
                        }
                    }
                }
            }


            for ( const auto & ts_gps : *gps_msgs )
            {
                const nav_msgs::Odometry & gps_msg = ts_gps.second;
                Sophus::SE3d raw_gps_pose;
                raw_gps_pose.translation() << gps_msg.pose.pose.position.x,gps_msg.pose.pose.position.y,gps_msg.pose.pose.position.z;
                raw_gps_pose.setQuaternion(Eigen::Quaterniond(gps_msg.pose.pose.orientation.w,gps_msg.pose.pose.orientation.x,gps_msg.pose.pose.orientation.y,gps_msg.pose.pose.orientation.z));
                const Sophus::SE3d curPose_imu_lidar = (raw_gps_pose * pose_lidar_imu.inverse());
                const Sophus::SE3d diffPose_lidar = (firstPose_lidar_imu * curPose_imu_lidar);
                Sophus::SE3d gps_pose = diffPose_lidar;
                const bool invalid_pos = gps_pose.translation().norm() > 10e4;
                if ( invalid_pos ) gps_pose.translation().setZero();
                if ( gps_msg.pose.covariance[0] < m_gps_invalid_threshold || invalid_pos )
                {
                    gm->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitX()*0.1)).cast<float>(), VisMesh::red );
                    gm->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitY()*0.1)).cast<float>(), VisMesh::blue );
                    gm->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitZ()*0.1)).cast<float>(), VisMesh::green );
                }
                else
                {
                    im->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitX()*0.1)).cast<float>(), VisMesh::red );
                    im->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitY()*0.1)).cast<float>(), VisMesh::blue );
                    im->addEdge ( gps_pose.translation().cast<float>(), (gps_pose * (Eigen::Vector3d::UnitZ()*0.1)).cast<float>(), VisMesh::green );
                }
            }
            gm->showAggregatedEdgeMesh(Sophus::SE3d(),"gps_valid");
            im->showAggregatedEdgeMesh(Sophus::SE3d(),"gps_invalid");
        }
    }

    showPoses( cur_pose_first_sensor, interp_pose_map_sensor, sceneCloud->m_stamp );
#endif
    ++m_scan_id;
    last_stamp = input_stamp;
}

void MarsSplineRegistrationAdaptor::triggerVis ( MarsMapPointCloud::Ptr sceneCloud )
{
    constexpr int maxAllowedSize = 1;
    std::unique_lock<std::mutex> lock(m_visCloudDequeMutex);
    m_visCloudDeque.emplace_back(sceneCloud);
    while ( m_visCloudDeque.size() > maxAllowedSize ) { m_visCloudDeque.pop_front(); }
    m_visCloudDequeCv.notify_one();
}

// visualization:
void MarsSplineRegistrationAdaptor::runVisMap ( )
{
#ifdef USE_EASY_PBR
    while ( m_run_vis_run )
#else
    while( m_run_vis_run && ros::ok() )
#endif
    {
        MarsMapPointCloud::Ptr sceneCloud = nullptr;
        // Wait until callback() sends data
        {
            std::unique_lock<std::mutex> lock ( m_visCloudDequeMutex );
            if ( !m_visCloudDeque.empty() )
            {
                sceneCloud = m_visCloudDeque.front();
                m_visCloudDeque.pop_front();
            }
        }

        // do the job
        if ( sceneCloud != nullptr )
        {
            ZoneScopedN("MarsSplineRegistrationAdaptor::Visualization");
            const int64_t now = sceneCloud->m_stamp;
            //LOG(1) << "now: " << now << " map_sub: " << m_map_publisher->getNumSubscribers() << " topic: " << m_map_publisher->getTopic() ;
            if ( m_map_publisher->getNumSubscribers() > 0 )
            {
                constexpr int64_t time_between_map_publish_ns = 1e9; // 1 sec
                static int64_t last_map_publishing_time = 0;
                static sensor_msgs::PointCloud2Ptr msg = nullptr;
                static int64_t prev_last_keyframe_stamp = m_marsRegistrator->getLastKeyFrameStamp();
                const int64_t last_keyframe_stamp = m_marsRegistrator->getLastKeyFrameStamp();
                //LOG(1) << "ts: " << prev_last_keyframe_stamp << " " << last_keyframe_stamp << " " << (msg==nullptr);
                if ( last_keyframe_stamp != prev_last_keyframe_stamp || msg == nullptr )
                {
                    MarsMapPointCloud::Ptr local_map_centers = m_marsRegistrator->getLocalMapSurfelCenters();
                    local_map_centers->transform_points( /*m_pose_field_map_init */ m_cur_pose_baselink_sensor);
                    local_map_centers->m_stamp = now;
                    msg = publishCloud<MarsMapPointCloud> ( m_map_publisher, local_map_centers, m_map_frame );
                    prev_last_keyframe_stamp = last_keyframe_stamp;
                }
                else
                {
                    if ( msg != nullptr && now > (last_map_publishing_time+time_between_map_publish_ns) )
                    {
                        msg->header.stamp.fromNSec(now);
                        publishCloud( m_map_publisher, msg );
                        last_map_publishing_time = now;
                    }
                }
            }

            if ( m_scene_publisher->getNumSubscribers() > 0 )
            {
                MarsMapPointCloud::Ptr scene_map_centers = m_marsRegistrator->getSceneMapSurfelCenters();
                scene_map_centers->transform_points( /*m_pose_field_map_init */ m_cur_pose_baselink_sensor);
                scene_map_centers->m_stamp = now;
                publishCloud<MarsMapPointCloud> ( m_scene_publisher, scene_map_centers, m_map_frame );
            }

            if ( m_scan_publisher->getNumSubscribers() > 0 )
            {
                MarsMapPointCloud::Ptr transformed_scene_cloud = MarsMapPointCloud::create( sceneCloud );
                //transformed_scene_cloud->transform_points(pose_world_map);
                transformed_scene_cloud->m_stamp = now;
                publishCloud<MarsMapPointCloud> ( m_scan_publisher, transformed_scene_cloud, m_corrected_sensor_frame ); // should be corrected_sensor_frame
            }

            if ( m_scan_downsampled_publisher->getNumSubscribers() > 0 )
            {
                MarsMapPointCloud::Ptr downsampled_scene_cloud = DownsampleCloud<MarsMapPointCloud>( sceneCloud, m_downsampled_scan_res );
                downsampled_scene_cloud->m_stamp = now;
                publishCloud<MarsMapPointCloud> ( m_scan_downsampled_publisher, downsampled_scene_cloud, m_sensor_frame );
            }
        }

        std::unique_lock<std::mutex> lk ( m_visCloudDequeMutex );
        if ( m_visCloudDeque.empty() )
            m_visCloudDequeCv.wait_for( lk, std::chrono::milliseconds(50) );
    }
}
