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
#include <memory>
#include <vector>
#include <limits>
#include "sophus/se3.hpp"

struct MarsPointCloudBase
{
    typedef std::shared_ptr<MarsPointCloudBase> Ptr;
    Eigen::Matrix3Xf m_points;
    Eigen::VectorXu m_scan_line_id;
    Eigen::VectorXi m_scan_time;
    Eigen::VectorXi m_scan_id;
    Eigen::VectorXf m_intensity;
    Eigen::VectorXu m_reflectivity;
    int m_last_added_point = 0;

    const Eigen::VectorXf & intensity() const{ return m_intensity; }
    const Eigen::VectorXu & reflectivity() const{ return m_reflectivity; }

    Sophus::SE3d m_pose;
    uint64_t m_stamp = 0;
    std::string m_frame_id = "";

    void resize( const int & n, const bool & all = true);
    void resizeAdditionals( const bool & has_intensity = false, const bool & reflectivity = false, const bool & time = false, const bool & has_line_id = true );
    int size() const;
    int capacity() const;
    void clear();

    void addPoint ( const Eigen::Vector3f & point );

    void add ( Ptr other, const Sophus::SE3d & tf_ww = Sophus::SE3d() );

    bool single_scan() const{ return m_scan_id.size() == 1; }

    const Eigen::Matrix3Xf & toMatrix3Xf() const { return m_points; }

    Sophus::SE3d getPose() const;
    void setPose( const Sophus::SE3d & tf );
    void transform_points( const Sophus::SE3d & tf );
    void transform_origin( const Eigen::Vector3d & orig );
    void transform_pose ( const Sophus::SE3d & tf );
    static Ptr create( Ptr other = nullptr);
};

struct MarsColorPointCloud : MarsPointCloudBase
{
    typedef std::shared_ptr<MarsColorPointCloud> Ptr;
    Eigen::Matrix3Xf m_color;
    void resize( const int & n, const bool & all = true);
    static Ptr create( Ptr other = nullptr);
};


struct MarsSemanticPointCloud : MarsPointCloudBase
{
    typedef std::shared_ptr<MarsSemanticPointCloud> Ptr;
    Eigen::VectorXi m_class;
    Eigen::MatrixXf m_semantic;
    void resize( const int & n, const bool & all = true);
    static Ptr create( Ptr other = nullptr);
};

template<typename T>
constexpr bool has_semantics_v = false;
template<>
constexpr bool has_semantics_v<MarsSemanticPointCloud> = true;

template<typename T>
constexpr bool has_scan_id_v = false;
template<>
constexpr bool has_scan_id_v<MarsPointCloudBase> = true;

template<typename T>
constexpr bool has_intensity_v = false;
template<>
constexpr bool has_intensity_v<MarsPointCloudBase> = true;

template<typename T>
constexpr bool has_color_v = false;
template<>
constexpr bool has_color_v<MarsColorPointCloud> = true;

template<typename T>
constexpr bool is_mars_cloud_v = false;
template<>
constexpr bool is_mars_cloud_v<MarsPointCloudBase> = true;
template<>
constexpr bool is_mars_cloud_v<MarsColorPointCloud> = true;
template<>
constexpr bool is_mars_cloud_v<MarsSemanticPointCloud> = true;

