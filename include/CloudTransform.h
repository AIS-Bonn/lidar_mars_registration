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

#ifdef USE_EASY_PBR
#include "easy_pbr/Scene.h"
static void compensateOrientation( MeshCloudPtr cloud, const Eigen::VectorXi & times, const int64_t & cur_time, const int64_t & last_time, const std::vector<sensor_msgs::Imu> & imu_msgs, const Eigen::Quaterniond & q_lidar_imu )
{
    if ( last_time == 0 ) return;
    if ( imu_msgs.empty() ) return;

    size_t firstPos = 0, lastPos = 0;
    for ( size_t idx = 0; idx < imu_msgs.size(); ++idx)
    {
        const int64_t imu_msg_ts = int64_t(imu_msgs[idx].header.stamp.toNSec());
        if ( imu_msg_ts < last_time )
            firstPos = idx;
        if ( imu_msg_ts < cur_time )
            lastPos = idx;
    }
    int32_t max_time = std::numeric_limits<int32_t>::min();
    int32_t min_time = std::numeric_limits<int32_t>::max();
    for ( int idx = 0; idx < times.rows(); ++idx )
    {
        if ( times(idx) > max_time ) max_time = times(idx);
        if ( times(idx) < min_time && cloud->D(idx) > 0.1 ) min_time = times(idx);
    }

    const double dt = max_time - min_time; // difference in nano sec
    const double dp = lastPos - firstPos;

    const Eigen::Quaterniond firstOrientation = Eigen::Quaterniond(imu_msgs[firstPos].orientation.w,imu_msgs[firstPos].orientation.x,imu_msgs[firstPos].orientation.y,imu_msgs[firstPos].orientation.z).normalized();
    const Eigen::Quaterniond firstOrientation_lidar_imu = (q_lidar_imu*firstOrientation.inverse()).normalized();

    std::vector<Eigen::Quaterniond> diffOri(lastPos-firstPos+1,Eigen::Quaterniond::Identity());
    for ( size_t idx = firstPos; idx <= lastPos; ++idx)
    {
        const Eigen::Quaterniond cur_ori = Eigen::Quaterniond(imu_msgs[idx].orientation.w,imu_msgs[idx].orientation.x,imu_msgs[idx].orientation.y,imu_msgs[idx].orientation.z).normalized();
        const Eigen::Quaterniond curOrientation_imu_lidar = (cur_ori * q_lidar_imu.inverse()).normalized();
        const Eigen::Quaterniond diffOrientation_lidar = (firstOrientation_lidar_imu * curOrientation_imu_lidar).normalized();
        diffOri[idx-firstPos] = diffOrientation_lidar;
    }


    Eigen::Quaterniond diffOrientation_lidar;
    int prevClosestIdx = -1;
    const int num_pts = cloud->V.rows();
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for ( int idx = 0; idx < num_pts; ++idx )
    {
        if ( cloud->D(idx) < 0.1 ) continue;

        const double ct = times(idx) - min_time;
        const double ctdt = ct / dt; // in [0,1]
        const int closestIdx = std::max<size_t>(firstPos,std::min<size_t>(std::round<size_t>(dp * ctdt)+firstPos, lastPos))-firstPos;
        const Eigen::Quaterniond & diffOrientation_lidar = diffOri[closestIdx];
        cloud->V.row(idx) = ( diffOrientation_lidar * cloud->V.row(idx).transpose()).transpose();
    }
    LOG(INFO) << "firstOri: "<< firstOrientation.coeffs().transpose() << " lastDiffOriL: " << diffOrientation_lidar.coeffs().transpose() << " q_li: " << q_lidar_imu.coeffs().transpose() << " #i: " << imu_msgs.size() << " dt: " << dt << " dp: " << dp << " +t: " << max_time << " -t: " << min_time << " lp: " << lastPos << " fp: " << firstPos;
}
#endif

static void compensateOrientation( MarsMapPointCloud::Ptr cloud, const Eigen::VectorXi & times, const int64_t & cur_time, const int64_t & last_time, const std::vector<sensor_msgs::Imu> & imu_msgs, const Eigen::Quaterniond & q_lidar_imu, const double & min_range = 0.1 )
{
    if ( last_time == 0 ) return;
    if ( imu_msgs.empty() ) return;

    size_t firstPos = 0, lastPos = 0;
    for ( size_t idx = 0; idx < imu_msgs.size(); ++idx)
    {
        const int64_t imu_msg_ts = int64_t(imu_msgs[idx].header.stamp.toNSec());
        if ( imu_msg_ts < last_time )
            firstPos = idx;
        if ( imu_msg_ts < cur_time )
            lastPos = idx;
    }
    int32_t max_time = std::numeric_limits<int32_t>::min();
    int32_t min_time = std::numeric_limits<int32_t>::max();

    Eigen::Matrix3Xf & pts = cloud->m_points;
    const int num_pts = cloud->size();
    const int num_times = times.rows();
    for ( int idx = 0; idx < num_times && idx < num_pts; ++idx )
    {
        if ( times(idx) > max_time ) max_time = times(idx);
        if ( times(idx) < min_time && pts.col(idx).norm() > min_range ) min_time = times(idx);
    }

    const double dt = max_time - min_time; // difference in nano sec
    const double dp = lastPos - firstPos;
    const Eigen::Quaterniond firstOrientation = Eigen::Quaterniond(imu_msgs[firstPos].orientation.w,imu_msgs[firstPos].orientation.x,imu_msgs[firstPos].orientation.y,imu_msgs[firstPos].orientation.z).normalized();
    const Eigen::Quaterniond firstOrientation_lidar_imu = (q_lidar_imu*firstOrientation.inverse()).normalized();

    std::vector<Eigen::Quaternionf> diffOri(lastPos-firstPos+1,Eigen::Quaternionf::Identity());
    for ( size_t idx = firstPos; idx <= lastPos; ++idx)
    {
        const Eigen::Quaterniond cur_ori = Eigen::Quaterniond(imu_msgs[idx].orientation.w,imu_msgs[idx].orientation.x,imu_msgs[idx].orientation.y,imu_msgs[idx].orientation.z).normalized();
        const Eigen::Quaterniond curOrientation_imu_lidar = (cur_ori * q_lidar_imu.inverse()).normalized();
        const Eigen::Quaterniond diffOrientation_lidar = (firstOrientation_lidar_imu * curOrientation_imu_lidar).normalized();
        diffOri[idx-firstPos] = diffOrientation_lidar.cast<float>();
    }

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for ( int idx = 0; idx < num_pts; ++idx )
    {
        if ( pts.col(idx).norm() < min_range ) continue;

        const int32_t ct = times(idx) - min_time;
        const double ctdt = ct / dt; // in [0,1]
        const int closestIdx = std::max<size_t>(firstPos,std::min<size_t>(std::round<size_t>(dp * ctdt)+firstPos, lastPos))-firstPos;

        if ( closestIdx < 0 || closestIdx >= diffOri.size() ) LOG(FATAL) << "ci: " << closestIdx << " firstPos: " << firstPos << " lastPos: " << lastPos << " do: " << diffOri.size();

        const Eigen::Quaternionf & diffOrientation_lidar = diffOri[closestIdx];
        pts.col(idx) = diffOrientation_lidar * pts.col(idx);
    }
}

#ifdef USE_EASY_PBR
static void copyCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, MeshCloudPtr cloud, Eigen::VectorXi * times = nullptr )
{
    if ( ! cloudMsg ) return;
    bool has_rgb = false;
    bool has_x = false;
    bool has_y = false;
    bool has_z = false;
    bool has_intensity = false;
    bool has_reflectivity = false;
    bool has_range = false;
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
    if ( ! has_x || ! has_y || ! has_z ) { LOG(INFO) << "cloud does not contain x,y,z!"; return; }
    const size_t numPoints = cloudMsg->height*cloudMsg->width;

    cloud->V.resize(numPoints, 3);
    if ( has_range ) cloud->D.resize(numPoints, 1);
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

        if ( has_range ) cloud->D.row(i) << cloud->V.row(i).norm();
        if ( has_intensity ) cloud->I.row(i) << *reinterpret_cast<const float*>(&cloudMsg->data[point_start + offset_intensity])/6000.0;
        if ( has_reflectivity ) cloud->C.row(i).setConstant((*reinterpret_cast<const uint16_t*>(&cloudMsg->data[point_start + offset_reflectivity])/ float(std::numeric_limits<uint16_t>::max() ) ) );
        if ( has_time && times != nullptr ) times->row(i) << *reinterpret_cast<const uint32_t*>(&cloudMsg->data[point_start + offset_time]);
        if ( has_semantic ) cloud->S_pred.row(i) = Eigen::Map<const Eigen::VectorXf>(reinterpret_cast<const float*>(&cloudMsg->data[point_start+offset_semantic]),width_semantic,1).cast<double>();
        if ( has_rgb ) cloud->C.row(i) = reinterpret_cast<const PointRGB*>(&cloudMsg->data[point_start+offset_rgb])->toVec3d();

        if ( !cloud->D.row(i).allFinite() ){
            cloud->V.row(i).setZero();
            if ( has_range ) cloud->D.row(i).setZero();
            if ( has_intensity ) cloud->I.row(i).setZero();
            if ( has_reflectivity ) cloud->C.row(i).setZero();
            if ( has_time && times != nullptr ) times->row(i).setZero();
        }
    }
    cloud->t = cloudMsg->header.stamp.toNSec();
}
#endif

template<typename MarsPointCloud>
static void copyCloud( sensor_msgs::PointCloud2::ConstPtr cloudMsg, typename MarsPointCloud::Ptr cloud, Eigen::VectorXi * times = nullptr )
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

