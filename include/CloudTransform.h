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
#include "MarsMapParameters.h"
#ifdef USE_EASY_PBR
#include "easy_pbr/Mesh.h"
namespace easy_pbr{
    class Mesh;
    typedef std::shared_ptr<Mesh> MeshSharedPtr;
}
typedef easy_pbr::Mesh MeshCloud;
typedef easy_pbr::MeshSharedPtr MeshCloudPtr;
#endif

#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <absl/container/flat_hash_map.h>

#define LOGURU_REPLACE_GLOG 1
#include "loguru.hpp"

inline bool is_imu_ori_valid ( const sensor_msgs::Imu & msg )
{
    return ! ((std::abs<double>(msg.orientation_covariance[0]-(-1)) < 1e-12) || Eigen::Quaterniond(msg.orientation.w,msg.orientation.x,msg.orientation.y,msg.orientation.z).squaredNorm() < 1e-4);
}
inline Sophus::SO3d ori_from_imu_msg ( const sensor_msgs::Imu & msg )
{
    return Sophus::SO3d ( Eigen::Quaterniond(msg.orientation.w,msg.orientation.x,msg.orientation.y,msg.orientation.z).normalized() );
}
inline Eigen::Vector3d get_gyr_from_msg ( const sensor_msgs::Imu & msg )
{
    return {msg.angular_velocity.x, msg.angular_velocity.y, msg.angular_velocity.z}; // already in FLU
}

#ifdef USE_EASY_PBR
#include "easy_pbr/Scene.h"

static Sophus::SO3d compensateOrientation( MeshCloudPtr cloud, const Eigen::VectorXi & times, const int64_t & cur_scan_ns, const int64_t & last_scan_ns, const std::map<int64_t,sensor_msgs::Imu> & imu_msgs, const Sophus::SO3d & q_lidar_imu, const float & min_range2 = 0.1, const bool & should_use_gyro_directly = false, const bool & stamp_at_beginning = false )
{
    //constexpr bool print_info = false;
    if ( last_scan_ns == 0 ) return;
    if ( imu_msgs.size() < 2 ) return;
    //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation");
    //StopWatch watch;
    int64_t first_imu_ns = imu_msgs.cbegin()->first;
    int64_t last_imu_ns = imu_msgs.crbegin()->first;

    for ( const std::pair<const int64_t, sensor_msgs::Imu> & it : imu_msgs )
    {
        const int64_t & imu_msg_ts = it.first;
        if ( imu_msg_ts < last_scan_ns )
            first_imu_ns = imu_msg_ts;
        if ( imu_msg_ts < cur_scan_ns )
            last_imu_ns = imu_msg_ts;
    }
    const int64_t ref_imu_ns = last_imu_ns; // since cur_scan_ns is at the end of the scan.
    const bool firstOneIsOlder = first_imu_ns < last_scan_ns;
    // if first msg is older than last time: it should only be used partially.
    // TODO: handle this correctly!

    const Eigen::VectorXi & pt_times_ns = times;

    const int num_pts = cloud->V.rows();
    const int num_times = pt_times_ns.rows();
    const bool t_cont = (pt_times_ns.head(num_pts-1).array() <= pt_times_ns.segment(1,num_pts-1).array()).all();
    const int min_time_ns = pt_times_ns.head(num_pts).minCoeff();
    const int max_time_ns = pt_times_ns.head(num_pts).maxCoeff();

    bool use_gyro_directly = should_use_gyro_directly;
    if ( !should_use_gyro_directly && !is_imu_ori_valid(imu_msgs.find(ref_imu_ns)->second) )
    {
        //LOG(WARNING) << "Reference orientation is invalid, but I should not use gyro directly? I still have to!";
        use_gyro_directly = true;
    }
    std::map<int64_t,std::pair<Sophus::SO3d, Eigen::Vector3d>> oris;

    const auto bit = imu_msgs.find(first_imu_ns);
    const auto lit = imu_msgs.find(last_imu_ns);
    if ( bit == imu_msgs.end() || lit == imu_msgs.end() )
    {
        //LOG(WARNING) << "could not find iterators. t1: " << first_imu_ns << " tn: " << last_imu_ns;
        return;
    }
    const auto nit = std::next(lit);
    {
        //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation::compute_log");
        Sophus::SO3d prev_estim;
        Sophus::SO3d prev_estim_lidar;
        int64_t prev_imu_ns = bit->first;
        for ( auto it = bit; it != nit; ++it )
        {
            const int64_t & cur_imu_ns = it->first;
            const sensor_msgs::Imu & msg = it->second;
            Sophus::SO3d cur_estim;
            if ( use_gyro_directly )
                cur_estim = prev_estim * Sophus::SO3d::exp( TimeConversion::to_s(cur_imu_ns - prev_imu_ns) * get_gyr_from_msg ( msg ));
            else
                cur_estim = ori_from_imu_msg ( msg );
            const Sophus::SO3d cur_estim_lidar = cur_estim * q_lidar_imu.inverse();
            const Eigen::Vector3d log_vec = (prev_estim_lidar.inverse()*cur_estim_lidar).log();
            oris[cur_imu_ns] = {cur_estim_lidar,log_vec};
            prev_estim_lidar = cur_estim_lidar;
            prev_estim = cur_estim;
            prev_imu_ns = cur_imu_ns;
        }
    }

    const Sophus::SO3d ref_ori_lidar = oris.crbegin()->second.first.inverse(); // q_lidar_imu * q_imu0_world = q_lidar0_world
    for ( std::pair<const int64_t, std::pair<Sophus::SO3d,Eigen::Vector3d>> & it : oris )
    {
        it.second.first = ref_ori_lidar * it.second.first; // premultiply ref_ori to get relative to last one.
    }
    auto it_cont = oris.cbegin();

    const int64_t offset = stamp_at_beginning ? 0 : -max_time_ns;  // max_time_ns is at cloud_stamp => go to first time.

    Sophus::SO3d diff_ori_lidar;
    //int num_slerped = 0;
    int prev_pt_ns = std::numeric_limits<int>::lowest(); // force computation, since times_ns[i] >= 0
    const int64_t first_scan_pt_ns = cur_scan_ns + offset;
    //{
    //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation::points");
    for ( int idx = 0; idx < num_pts; ++idx )
    {
        if ( cloud->V.row(idx).squaredNorm() < min_range2 ) continue;
        if ( prev_pt_ns != pt_times_ns(idx) )
        {
            prev_pt_ns = pt_times_ns(idx);
            const int64_t cur_pt_ns = first_scan_pt_ns + pt_times_ns(idx);

            auto it_next = it_cont;
            if ( t_cont )
            {
                if ( (it_next != oris.cend()) && (it_next->first <= cur_pt_ns) )
                {
                    it_next = std::next(it_next);
                    it_cont = it_next;
                }
            }
            else
            {
                it_next = oris.upper_bound(cur_pt_ns); // => greater than
            }

            if ( it_next == oris.cend() ) // none larger: get last one, we will cap time
            {
                diff_ori_lidar = Sophus::SO3d();
            }
            else
            {
                if ( it_next == oris.cbegin() ) // first one is larger:
                {
                    diff_ori_lidar = oris.cbegin()->second.first;// due to time cap, we use first
                }
                else
                {
                    // somewhere inbetween
                    const auto it_prev = std::prev(it_next);
                    const int64_t next_imu_ns = it_next->first;
                    const int64_t prev_imu_ns = it_prev->first;
                    //LOG(INFO) << "idx: " << idx << " oi: " << prevOriIdx << " i: " << imu_msgs.size() << " s: " << oris.size() << " dt: " << dt<< " t: " << cur_t_ns << " o: " << oris[prevOriIdx].first << " n: " << oris[prevOriIdx+1].first << " t-o: " << (cur_t_ns - oris[prevOriIdx].first) << " n-t: " << ( oris[prevOriIdx+1].first - cur_t_ns);
                    if( cur_pt_ns == prev_imu_ns ) // no need to check for next_imu_ns, due to upper_bound
                        diff_ori_lidar = it_prev->second.first;
                    else {
                        //SLERP
                        const int64_t cur_dt_ns = cur_pt_ns - prev_imu_ns;
                        const int64_t imu_dt_ns = next_imu_ns - prev_imu_ns;
                        const float dt = std::min(1.f, std::max( 0.f, float(TimeConversion::to_s(cur_dt_ns) / TimeConversion::to_s(imu_dt_ns))));
                        diff_ori_lidar = (it_prev->second.first * Sophus::SO3d::exp( dt * it_prev->second.second ));
                        //++num_slerped;
                    }
                }
            }
        }
        cloud->V.row(idx) = (diff_ori_lidar.unit_quaternion() * cloud->V.row(idx).transpose()).transpose();
    }
    //}
    //if constexpr ( print_info )
    //LOG(1) << "orientation comp took: " << watch.getTime() << " i: " << imu_msgs.size() << " o: " << oris.size() << " s: " << num_slerped;
    if ( !oris.empty() )
        return oris.crbegin()->second.first;
    else
        return Sophus::SO3d();
}
#endif

static Sophus::SO3d compensateOrientation( MarsMapPointCloud::Ptr cloud, const Eigen::VectorXi & times, const int64_t & cur_scan_ns, const int64_t & last_scan_ns, const std::map<int64_t,sensor_msgs::Imu> & imu_msgs, const Sophus::SO3d & q_lidar_imu, const float & min_range2 = 0.1, const bool & should_use_gyro_directly = false, const bool & stamp_at_beginning = false )
{
    static Sophus::SO3d invalid_quaternion;
    if ( invalid_quaternion.data()[3] > 0.1 ) invalid_quaternion.data()[3] = 0;

    //constexpr bool print_info = false;
    if ( last_scan_ns == 0 ) return invalid_quaternion;
    if ( imu_msgs.size() < 2 ) return invalid_quaternion;
    //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation");
    //StopWatch watch;
    int64_t first_imu_ns = imu_msgs.cbegin()->first;
    int64_t last_imu_ns = imu_msgs.crbegin()->first;

    //LOG(1) << "imu("<<imu_msgs.size()<<" times: " << first_imu_ns << " last: " << last_imu_ns << " scan: " << cur_scan_ns << " last: " << last_scan_ns << " first<scan? " << (first_imu_ns < cur_scan_ns) << " last<scan? " << (last_imu_ns < cur_scan_ns) << " Last: first<scan? " << (first_imu_ns < last_scan_ns) << " last<scan? " << (last_imu_ns < last_scan_ns) ;

    for ( const std::pair<const int64_t, sensor_msgs::Imu> & it : imu_msgs )
    {
        const int64_t & imu_msg_ts = it.first;
        if ( imu_msg_ts < last_scan_ns )
            first_imu_ns = imu_msg_ts;
        if ( imu_msg_ts < cur_scan_ns )
            last_imu_ns = imu_msg_ts;
    }
    const int64_t ref_imu_ns = last_imu_ns; // since cur_scan_ns is at the end of the scan.
    const bool firstOneIsOlder = first_imu_ns < last_scan_ns;
    // if first msg is older than last time: it should only be used partially.
    // TODO: handle this correctly!

    if ( firstOneIsOlder )
    {
        LOG(WARNING) << "OriComp: First imu is older than last scan: " << first_imu_ns << " < " << last_scan_ns << " skipping and no compensation done.";
        return invalid_quaternion;
    }

    //LOG(1) << "After start: imu("<<imu_msgs.size()<<" times: " << first_imu_ns << " last: " << last_imu_ns << " scan: " << cur_scan_ns << " last: " << last_scan_ns << " first<scan? " << (first_imu_ns < cur_scan_ns) << " last<scan? " << (last_imu_ns < cur_scan_ns) << " Last: first<scan? " << (first_imu_ns < last_scan_ns) << " last<scan? " << (last_imu_ns < last_scan_ns) ;


    Eigen::Matrix3Xf & pts = cloud->m_points;
    const Eigen::VectorXi & pt_times_ns = times;

    const int num_pts = cloud->size();
    const int num_times = pt_times_ns.rows();
    const bool t_cont = (pt_times_ns.head(num_pts-1).array() <= pt_times_ns.segment(1,num_pts-1).array()).all();
    const int min_time_ns = pt_times_ns.head(num_pts).minCoeff();
    const int max_time_ns = pt_times_ns.head(num_pts).maxCoeff();

    bool use_gyro_directly = should_use_gyro_directly;
    if ( !should_use_gyro_directly && !is_imu_ori_valid(imu_msgs.find(ref_imu_ns)->second) )
    {
        //LOG(WARNING) << "Reference orientation is invalid, but I should not use gyro directly? I still have to!";
        use_gyro_directly = true;
    }
    std::map<int64_t,std::pair<Sophus::SO3d, Eigen::Vector3d>> oris;

    const auto bit = imu_msgs.find(first_imu_ns);
    const auto lit = imu_msgs.find(last_imu_ns);
    if ( bit == imu_msgs.end() || lit == imu_msgs.end() )
    {
        LOG(WARNING) << "OriComp: could not find iterators. t1: " << first_imu_ns << " tn: " << last_imu_ns;
        return invalid_quaternion;
    }
    const auto nit = std::next(lit);
    {
        //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation::compute_log");
        Sophus::SO3d prev_estim;
        Sophus::SO3d prev_estim_lidar;
        int64_t prev_imu_ns = bit->first;
        for ( auto it = bit; it != nit; ++it )
        {
            const int64_t & cur_imu_ns = it->first;
            const sensor_msgs::Imu & msg = it->second;
            Sophus::SO3d cur_estim;
            if ( use_gyro_directly )
                cur_estim = prev_estim * Sophus::SO3d::exp( TimeConversion::to_s(cur_imu_ns - prev_imu_ns) * get_gyr_from_msg ( msg ));
            else
                cur_estim = ori_from_imu_msg ( msg );
            //LOG(1) << "estim: " << cur_estim.params().transpose()<< " prev: " << prev_estim.params().transpose() <<  " at t: " << cur_imu_ns  << " ( " << msg.header.stamp.toNSec() << " ) p: " << prev_imu_ns  << " d: " << std::distance(bit,nit);
            const Sophus::SO3d cur_estim_lidar = cur_estim * q_lidar_imu.inverse();
            const Eigen::Vector3d log_vec = (prev_estim_lidar.inverse()*cur_estim_lidar).log();
            oris[cur_imu_ns] = {cur_estim_lidar,log_vec};
            prev_estim_lidar = cur_estim_lidar;
            prev_estim = cur_estim;
            prev_imu_ns = cur_imu_ns;
        }
    }

    const Sophus::SO3d ref_ori_lidar = oris.crbegin()->second.first.inverse(); // q_lidar_imu * q_imu0_world = q_lidar0_world
    for ( std::pair<const int64_t, std::pair<Sophus::SO3d,Eigen::Vector3d>> & it : oris )
    {
        it.second.first = ref_ori_lidar * it.second.first; // premultiply ref_ori to get relative to last one.
    }
    auto it_cont = oris.cbegin();

    const int64_t offset = stamp_at_beginning ? 0 : -max_time_ns;  // max_time_ns is at cloud_stamp => go to first time.

    Sophus::SO3f diff_ori_lidar;
    //int num_slerped = 0;
    int prev_pt_ns = std::numeric_limits<int>::lowest(); // force computation, since times_ns[i] >= 0
    const int64_t first_scan_pt_ns = cur_scan_ns + offset;
    //{
    //ZoneScopedN("MotionCompensation::Cloud::compensateOrientation::points");
    for ( int idx = 0; idx < num_pts; ++idx )
    {
        if ( pts.col(idx).squaredNorm() < min_range2 ) continue;
        if ( prev_pt_ns != pt_times_ns(idx) )
        {
            prev_pt_ns = pt_times_ns(idx);
            const int64_t cur_pt_ns = first_scan_pt_ns + pt_times_ns(idx);

            auto it_next = it_cont;
            if ( t_cont )
            {
                if ( (it_next != oris.cend()) && (it_next->first <= cur_pt_ns) )
                {
                    it_next = std::next(it_next);
                    it_cont = it_next;
                }
            }
            else
            {
                it_next = oris.upper_bound(cur_pt_ns); // => greater than
            }

            if ( it_next == oris.cend() ) // none larger: get last one, we will cap time
            {
                diff_ori_lidar = Sophus::SO3f();
            }
            else
            {
                if ( it_next == oris.cbegin() ) // first one is larger:
                {
                    diff_ori_lidar = oris.cbegin()->second.first.template cast<float>();// due to time cap, we use first
                }
                else
                {
                    // somewhere inbetween
                    const auto it_prev = std::prev(it_next);
                    const int64_t next_imu_ns = it_next->first;
                    const int64_t prev_imu_ns = it_prev->first;
                    //LOG(INFO) << "idx: " << idx << " oi: " << prevOriIdx << " i: " << imu_msgs.size() << " s: " << oris.size() << " dt: " << dt<< " t: " << cur_t_ns << " o: " << oris[prevOriIdx].first << " n: " << oris[prevOriIdx+1].first << " t-o: " << (cur_t_ns - oris[prevOriIdx].first) << " n-t: " << ( oris[prevOriIdx+1].first - cur_t_ns);
                    if( cur_pt_ns == prev_imu_ns ) // no need to check for next_imu_ns, due to upper_bound
                        diff_ori_lidar = it_prev->second.first.template cast<float>();
                    else {
                        //SLERP
                        const int64_t cur_dt_ns = cur_pt_ns - prev_imu_ns;
                        const int64_t imu_dt_ns = next_imu_ns - prev_imu_ns;
                        const float dt = std::min(1.f, std::max( 0.f, float(TimeConversion::to_s(cur_dt_ns) / TimeConversion::to_s(imu_dt_ns))));
                        diff_ori_lidar = (it_prev->second.first * Sophus::SO3d::exp( dt * it_prev->second.second )).template cast<float>();
                        //++num_slerped;
                    }
                }
            }
        }
        pts.col(idx) = diff_ori_lidar * pts.col(idx);
    }

    //if ( ! oris.empty() ) { LOG(1) << "ori beg: " << oris.cbegin()->second.first.params().transpose() << " ori end: " << oris.crbegin()->second.first.params().transpose(); }
    //}
    //if constexpr ( print_info )
    //LOG(1) << "orientation comp took: " << watch.getTime() << " i: " << imu_msgs.size() << " o: " << oris.size() << " s: " << num_slerped;
    const Sophus::SO3d rel_ori = oris.empty() ? invalid_quaternion : oris.cbegin()->second.first.inverse();
    //LOG(1) << "rel ori: " << rel_ori.params().transpose() << " empty? " << oris.empty() << " ref: " << ref_ori_lidar.params().transpose() << " gyro? " << use_gyro_directly;
    return rel_ori;
}

#ifdef USE_EASY_PBR
static void copyCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, MeshCloudPtr cloud, Eigen::VectorXi * times = nullptr, const bool & stamp_at_front = false )
{
    if ( ! cloudMsg ) return;
    bool has_rgb = false;
    bool has_x = false;
    bool has_y = false;
    bool has_z = false;
    bool has_intensity = false;
    bool has_reflectivity = false;
    bool has_time = true;
    bool has_semantic = false;
    uint32_t offset_x = 0;
    uint32_t offset_y = 0;
    uint32_t offset_z = 0;
    uint32_t offset_intensity = 0;
    uint32_t offset_reflectivity = 0;
    uint32_t offset_rgb = 0;
    uint32_t offset_time = 0;
    uint32_t offset_semantic = 0;
    uint32_t width_semantic = 0;
    for(size_t i=0; i<cloudMsg->fields.size(); ++i)
    {
        if ( cloudMsg->fields[i].name=="x" )
        {
            has_x = true;
            offset_x = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="y" )
        {
            has_y = true;
            offset_y = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="z" )
        {
            has_z = true;
            offset_z = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="rgb" )
        {
            has_rgb = true;
            offset_rgb = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="intensity" )
        {
            has_intensity = true;
            offset_intensity = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="reflectivity" )
        {
            has_reflectivity = true;
            offset_reflectivity = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="t" )
        {
            has_time = true;
            offset_time = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name == "semantic")
        {
            has_semantic = true;
            offset_semantic = cloudMsg->fields[i].offset;
            width_semantic = cloudMsg->fields[i].count;
        }
    }
    if ( ! has_x || ! has_y || ! has_z ) { LOG(INFO) << "cloud does not contain x,y,z!"; return; }
    const size_t numPoints = cloudMsg->height*cloudMsg->width;

    cloud->V.resize(numPoints, 3);
    cloud->D.resize(numPoints, 1);
    if ( has_intensity ) cloud->I.resize(numPoints, 1);
    if ( has_reflectivity ) cloud->C.resize(numPoints, 3);
    if ( has_time && times != nullptr ) times->resize(numPoints,1);
    if ( has_semantic ) cloud->S_pred.resize(numPoints,width_semantic);
    if ( has_rgb ) cloud->C.resize(numPoints, 3);

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
      Eigen::Vector3d toVec3d() const { return Eigen::Vector3d(r,g,b); }
    };

#ifdef USE_OMP
    #pragma omp parallel for num_threads ( 3 )
#endif
    for (size_t i = 0; i < numPoints; ++i)
    {
        const uint32_t point_start = cloudMsg->point_step * i;
        cloud->V.row(i) << *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_x]),
                           *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_y]),
                           *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_z]);

        cloud->D.row(i) << cloud->V.row(i).norm();
        if ( has_intensity ) cloud->I.row(i) << *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_intensity])/6000.0;
        if ( has_reflectivity ) cloud->C.row(i).setConstant((*reinterpret_cast<const uint16_t*>(&cloudMsg->data[point_start + offset_reflectivity])/ float(std::numeric_limits<uint16_t>::max() ) ) );
        if ( has_time && times != nullptr ) times->row(i) << *reinterpret_cast<const uint32_t*>(&cloudMsg->data[point_start + offset_time]);
        if ( has_semantic ) cloud->S_pred.row(i) = Eigen::Map<const Eigen::VectorXf>(reinterpret_cast<const float*>(&cloudMsg->data[point_start+offset_semantic]),width_semantic,1).cast<double>();
        if ( has_rgb ) cloud->C.row(i) = reinterpret_cast<const PointRGB*>(&cloudMsg->data[point_start+offset_rgb])->toVec3d();

        if ( !cloud->D.row(i).allFinite() ){
            cloud->V.row(i).setZero();
            cloud->D.row(i).setZero();
            if ( has_intensity ) cloud->I.row(i).setZero();
            if ( has_reflectivity ) cloud->C.row(i).setZero();
            if ( has_time && times != nullptr ) times->row(i).setZero();
        }
    }
    cloud->t = cloudMsg->header.stamp.toNSec();
    if ( has_time && times != nullptr && stamp_at_front )
        cloud->t += times->maxCoeff();
}
#endif

template<typename MarsPointCloud>
static void copyCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, typename MarsPointCloud::Ptr cloud, Eigen::VectorXi * times = nullptr, const bool & stamp_at_front = false )
{
    if ( ! cloudMsg ) return;
    bool has_rgb = false;
    bool has_x = false;
    bool has_y = false;
    bool has_z = false;
    bool has_intensity = false;
    bool has_reflectivity = false;
//    bool has_ring = false;
    bool has_range = false;
    bool has_time = true;
    bool has_semantic = false;
    uint32_t offset_x = 0;
    uint32_t offset_y = 0;
    uint32_t offset_z = 0;
    uint32_t offset_rgb = 0;
    uint32_t offset_intensity = 0;
    uint32_t offset_reflectivity = 0;
//    uint32_t offset_ring = 0;
//    uint32_t offset_range = 0;
    uint32_t offset_time = 0;
    uint32_t offset_semantic = 0;
    uint32_t width_semantic = 0;
    for(size_t i=0; i<cloudMsg->fields.size(); ++i)
    {
        if ( cloudMsg->fields[i].name=="x" )
        {
            has_x = true;
            offset_x = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="y" )
        {
            has_y = true;
            offset_y = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="z" )
        {
            has_z = true;
            offset_z = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="rgb" )
        {
            has_rgb = true;
            offset_rgb = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="intensity" )
        {
            has_intensity = true;
            offset_intensity = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="reflectivity" )
        {
            has_reflectivity = true;
            offset_reflectivity = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name=="range" )
        {
            has_range = true;
        }
        if ( cloudMsg->fields[i].name=="t" )
        {
            has_time = true;
            offset_time = cloudMsg->fields[i].offset;
        }
        if ( cloudMsg->fields[i].name == "semantic")
        {
            has_semantic = true;
            offset_semantic = cloudMsg->fields[i].offset;
            width_semantic = cloudMsg->fields[i].count;
        }
    }
    if ( ! has_x || ! has_y || ! has_z ) return;

    const int numPoints = cloudMsg->height * cloudMsg->width;
    cloud->resize( numPoints, false );
    cloud->resizeAdditionals( has_intensity, has_reflectivity, has_time );

    Eigen::Matrix3Xf * colors = nullptr;
    if constexpr ( has_color_v<MarsPointCloud> )
        if ( has_rgb )
        {
            colors = &cloud->m_colors;
            cloud->m_colors.resize(3,numPoints);
        }
    Eigen::MatrixXf * semantics = nullptr;
    if  constexpr( has_semantics_v<MarsPointCloud>)
        if ( has_semantic )
        {
            semantics = &cloud->m_semantic;
            cloud->m_class.resize(numPoints,1);
            semantics->resize(width_semantic,numPoints);
        }
    Eigen::Matrix3Xf & pts = cloud->m_points;
    Eigen::VectorXi & scan_time = cloud->m_scan_time;
    cloud->m_scan_id = Eigen::VectorXi::Zero(1,1);
    Eigen::VectorXf & intensity = cloud->m_intensity;
    Eigen::VectorXu & reflectivity = cloud->m_reflectivity;
    if ( has_time ) scan_time.setZero();
    if ( has_intensity ) intensity.setZero();
    if ( has_reflectivity ) reflectivity.setZero();

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
      Eigen::Vector3d toVec3d() const { return Eigen::Vector3d(r,g,b); }
    };

    //#pragma omp parallel for num_threads ( 3 )
    for (int i = 0; i < numPoints; ++i)
    {
        const uint32_t point_start = cloudMsg->point_step * i;

        pts.col(i) << *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_x]),
                      *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_y]),
                      *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_z]);


        if ( has_intensity ) intensity(i) = *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_intensity])/6000.0;
        if ( has_reflectivity ) reflectivity(i) = *reinterpret_cast<const uint16_t*>(&cloudMsg->data[point_start + offset_reflectivity])/std::numeric_limits<uint16_t>::max();
        if ( has_time ) scan_time(i) = *reinterpret_cast<const uint32_t*>(&cloudMsg->data[point_start + offset_time]);
        if constexpr (has_semantics_v<MarsPointCloud>)
                if ( has_semantic ) semantics->col(i) = Eigen::Map<const Eigen::VectorXf>(reinterpret_cast<const float*>(&cloudMsg->data[point_start+offset_semantic]),width_semantic,1);
        if constexpr (has_color_v<MarsPointCloud>)
                if ( has_rgb ) colors->col(i) = reinterpret_cast<const PointRGB*>(&cloudMsg->data[point_start+offset_rgb])->toVec3d();

        if ( !pts.col(i).allFinite() ){
            pts.col(i).setZero();
            if ( has_intensity ) intensity(i) = 0;
            if ( has_reflectivity ) reflectivity(i) = 0;
            if ( has_time ) scan_time(i) = 0;
        }
    }
    cloud->m_last_added_point = numPoints;
    cloud->m_stamp = cloudMsg->header.stamp.toNSec();
    if ( has_time && times != nullptr ) *times = cloud->m_scan_time;
    if ( has_time && stamp_at_front ) // moves stamp to back like we need!
        cloud->m_stamp += cloud->m_scan_time.maxCoeff();
}

template<typename MarsPointCloud>
static typename MarsPointCloud::Ptr DownsampleCloud ( typename MarsPointCloud::Ptr cloud, const float & res )
{
    struct MeanComputation
    {
        int num = 0;
        Eigen::Vector3f m = Eigen::Vector3f::Zero();
        Eigen::Vector3f m_center = Eigen::Vector3f::Zero();
        MeanComputation() {}
        void add (const MeanComputation & other)
        {
            if ( num == 0 )
            {
                num = other.num;
                m = other.m;
            }
            else
            {
                if ( num < 20000 )
                {
                    m+=other.m;
                    num+=other.num;
                }
            }
        }
        void add ( const Eigen::Vector3f & pt )
        {
            ++num;
            m+=pt;
        }
        Eigen::Vector3f mean() const
        {
            return m_center + m / num;
        }
    };
    constexpr float map_side_length = 100.;
    constexpr int num_levels = 1;
    const int num_cells = std::round<int>(map_side_length / res);
    const int num_used_cells = std::min<int>(num_cells,1024);
    const float used_side_length = num_used_cells * res;
    const MapParameters params = MapParameters::getMapParameters( used_side_length, num_used_cells, num_levels);

    const Eigen::Matrix3Xf & pts = cloud->toMatrix3Xf();
    const int num_pts = cloud->size();
    typedef int CellIndexType;
    absl::flat_hash_map<CellIndexType,MeanComputation> cells;
    cells.reserve(num_pts);
    for ( int i = 0; i < num_pts; ++i )
    {
        if ( ! params.isInBounds<float>( pts.col(i), 0) ) continue;
        const Eigen::Vector3i cellIndex = params.toCellIndexVector<float>( pts.col(i), 0 );
        const CellIndexType cellIdx = params.toCellIndexFromIndexVector ( cellIndex );
        if ( cells[cellIdx].num == 0 ) cells[cellIdx].m_center = params.fromCellIdxToCenter<float>(cellIndex,0);
        cells[cellIdx].add(pts.col(i) - cells[cellIdx].m_center);
    }

    typename MarsPointCloud::Ptr cloud_downsampled = MarsPointCloud::create();
    cloud_downsampled->resize(cells.size());
    for ( const auto & cell : cells )
        cloud_downsampled->addPoint(cell.second.mean());
    return cloud_downsampled;
}

