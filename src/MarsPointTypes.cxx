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
#include "MarsPointTypes.h"

void MarsPointCloudBase::add ( MarsPointCloudBase::Ptr other, const Sophus::SE3d & tf_ww )
{
    if ( ! other ) return;
    MarsPointCloudBase::Ptr tmp = MarsPointCloudBase::create(other);
    const Sophus::SE3d tf = m_pose.inverse()*tf_ww*other->m_pose;
    tmp->transform_points(tf);
    const int old_num_pts =  size();
    const int tmp_num_pts =  tmp->size();
    m_points.conservativeResize(3,old_num_pts + tmp_num_pts);
    m_points.block(0,old_num_pts,3,tmp_num_pts) = tmp->m_points.leftCols(tmp_num_pts);
    m_last_added_point = m_points.cols();
}

void MarsPointCloudBase::addPoint ( const Eigen::Vector3f & point )
{
    if ( m_last_added_point+1 > capacity() )
    {
        const int tmp = m_last_added_point;
        resize(m_last_added_point+1); // resets m_last_added_point
        m_last_added_point = tmp;
    }
    m_points.col(m_last_added_point) = point;
    ++m_last_added_point;
}

int MarsPointCloudBase::capacity() const { return m_points.cols(); }
int MarsPointCloudBase::size() const { return m_last_added_point; }

inline void MarsPointCloudBase::clear() { m_points = Eigen::Matrix3Xf(); m_pose = Sophus::SE3d(); m_last_added_point = 0; }

void MarsPointCloudBase::resize( const int & n, const bool & all )
{
    m_last_added_point = 0;
    m_points.conservativeResize(3,n);
    if ( all )
        resizeAdditionals(true, true, true, true);
}

void MarsPointCloudBase::resizeAdditionals( const bool & has_intensity, const bool & has_reflectivity, const bool & has_time, const bool & has_line_id  )
{
    const int n = m_points.cols();
    if ( has_intensity ) m_intensity.conservativeResize(n,1);
    if ( has_reflectivity ) m_reflectivity.conservativeResize(n,1);
    if ( has_time ) m_scan_time.conservativeResize(n,1);
    if ( has_line_id ) m_scan_line_id.conservativeResize(n,1);
}

Sophus::SE3d MarsPointCloudBase::getPose() const { return m_pose; }
void MarsPointCloudBase::setPose( const Sophus::SE3d & tf ) { m_pose = tf; }

void MarsPointCloudBase::transform_points( const Sophus::SE3d & tf )
{
    m_points = tf.so3().matrix().cast<float>() * m_points;
    m_points.colwise() += tf.translation().cast<float>();
    transform_pose(tf);
}

inline void MarsPointCloudBase::transform_origin( const Eigen::Vector3d & orig )
{
    m_points.colwise() += orig.cast<float>();
    m_pose.translation() += orig;
}
inline void MarsPointCloudBase::transform_pose ( const Sophus::SE3d & tf )
{
    m_pose = tf * m_pose;
}

typename MarsPointCloudBase::Ptr MarsPointCloudBase::create( MarsPointCloudBase::Ptr other )
{
    Ptr ptr = std::make_shared<MarsPointCloudBase>();
    if ( other != nullptr )
    {
        ptr->m_points = other->m_points;
        ptr->m_scan_line_id = other->m_scan_line_id;
        ptr->m_scan_time = other->m_scan_time;
        ptr->m_scan_id = other->m_scan_id;
        ptr->m_last_added_point = other->m_last_added_point;

        ptr->m_intensity = other->m_intensity;
        ptr->m_reflectivity = other->m_reflectivity;

        ptr->m_pose = other->m_pose;
        ptr->m_stamp = other->m_stamp;
        ptr->m_frame_id = other->m_frame_id;
    }
    return ptr;
}
typename MarsColorPointCloud::Ptr MarsColorPointCloud::create( MarsColorPointCloud::Ptr other )
{
    Ptr ptr = std::make_shared<MarsColorPointCloud>();
    if ( other != nullptr )
    {
        ptr->m_points = other->m_points;
        ptr->m_scan_line_id = other->m_scan_line_id;
        ptr->m_scan_time = other->m_scan_time;
        ptr->m_scan_id = other->m_scan_id;
        ptr->m_last_added_point = other->m_last_added_point;

        ptr->m_color = other->m_color;
        ptr->m_intensity = other->m_intensity;
        ptr->m_reflectivity = other->m_reflectivity;

        ptr->m_pose = other->m_pose;
        ptr->m_stamp = other->m_stamp;
        ptr->m_frame_id = other->m_frame_id;
    }
    return ptr;
}
typename MarsSemanticPointCloud::Ptr MarsSemanticPointCloud::create( MarsSemanticPointCloud::Ptr other )
{
    Ptr ptr = std::make_shared<MarsSemanticPointCloud>();
    if ( other != nullptr )
    {
        ptr->m_points = other->m_points;
        ptr->m_scan_line_id = other->m_scan_line_id;
        ptr->m_scan_time = other->m_scan_time;
        ptr->m_scan_id = other->m_scan_id;
        ptr->m_last_added_point = other->m_last_added_point;

        ptr->m_class = other->m_class;
        ptr->m_semantic = other->m_semantic;
        ptr->m_intensity = other->m_intensity;
        ptr->m_reflectivity = other->m_reflectivity;

        ptr->m_pose = other->m_pose;
        ptr->m_stamp = other->m_stamp;
        ptr->m_frame_id = other->m_frame_id;
    }
    return ptr;
}

void MarsColorPointCloud::resize( const int & n, const bool & all )
{
    MarsPointCloudBase::resize(n,all);
    m_color.resize(3,n);
}
void MarsSemanticPointCloud::resize( const int & n, const bool & all )
{
    constexpr int m = 15;
    MarsPointCloudBase::resize(n,all);
    m_class.resize(n);
    m_semantic.resize(m, n);
}
