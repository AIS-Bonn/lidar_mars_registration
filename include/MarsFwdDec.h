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

#include "MarsPointTypes.h"
#include <memory>

typedef MarsPointCloudBase MarsMapPointCloud;
typedef std::shared_ptr<MarsMapPointCloud> MarsMapPointCloudPtr;

class DenseVoxelGrid;
class SparseVoxelGrid;
class BlockSparseVoxelGrid;
template<int P>
class PermutohedralLattice;
typedef PermutohedralLattice<3> Lattice;


template<typename T>
struct MarsMapBase;

typedef MarsMapBase<SparseVoxelGrid> SparseMarsMap;
typedef MarsMapBase<BlockSparseVoxelGrid> BlockSparseMarsMap;
//typedef MarsMapBase<DenseVoxelGrid> DenseMarsMap;
typedef MarsMapBase<Lattice> LatticeMarsMap;
//typedef DenseMarsMap MarsMap;
//typedef BlockSparseMarsMap MarsMap;
//typedef SparseMarsMap MarsMap;
typedef LatticeMarsMap MarsMap;

typedef MarsMap MarsMapType;
typedef std::shared_ptr<MarsMapType> MarsMapTypePtr;

class MarsMapWindow;
typedef std::shared_ptr<MarsMapWindow> MarsMapWindowPtr;

#define OPT_THREAD_NUM 3
//constexpr int MarsSplineDegree = 2;
constexpr int MarsSplineDegree = 3;
//constexpr int MarsSplineDegree = 4;

namespace TimeConversion
{
    constexpr double ns_to_s = 1e-9;  ///< Nanosecond to second conversion
    constexpr double s_to_ns = 1e9;   ///< Second to nanosecond conversion
    inline double to_s ( const int64_t & ns )
    {
        return ns * ns_to_s;
    }
    inline int64_t to_ns ( const double & s )
    {
        return static_cast<int64_t>( s * s_to_ns );
    }

    inline Eigen::VectorXd to_s ( const Eigen::VectorXt & ns )
    {
        return ns.cast<double>() * TimeConversion::ns_to_s;
    }
    inline Eigen::VectorXt to_ns ( const Eigen::VectorXd & s )
    {
        return ( s * TimeConversion::s_to_ns ).cast<int64_t>();
    }
    inline Eigen::VectorXt to_ns ( const Eigen::VectorXt & s )
    {
        return s * int64_t(TimeConversion::s_to_ns);
    }
}

namespace basalt
{
    template<int _N, typename _Scalar>
    struct Se3Spline;
}
typedef basalt::Se3Spline<MarsSplineDegree,double> MarsSplineType;
typedef basalt::Se3Spline<1,double> MarsConstraintSplineType;
typedef std::shared_ptr<MarsSplineType> MarsSplineTypePtr;
typedef std::shared_ptr<MarsConstraintSplineType> MarsConstraintSplineTypePtr;
