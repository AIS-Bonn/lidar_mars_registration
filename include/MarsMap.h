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

#include "MarsPointTypes.h"
#include "MarsMapParameters.h"
#include "MarsFwdDec.h"
#include "MarsSurfel.h"
#include "MarsCell.h"

#include "MarsDenseVoxelGrid.h"
#include "MarsSparseVoxelGrid.h"
#include "MarsBlockSparseVoxelGrid.h"
#include "MarsPermutohedralLattice.h"

class UID{
public :
    static int get();
};

template<typename T>
struct MarsMapBase
{
    typedef Cell CellType;
    typedef SurfelInfo::SurfelPtr SurfelPtr;
    typedef std::shared_ptr<MarsMapBase<T>> Ptr;
    typedef typename MapParameters::IndexType IndexType;
    typedef typename MapParameters::LevelIndexType LevelIndexType;

    T m_storage;

    bool m_adaptor_updated = false;
    std::uint32_t m_id = 0;

    Sophus::SE3d m_pose_init_w;
    SurfelInfoVector m_surfels;

    size_t m_num_points = 0;
    const MapParameters m_map_params;
    Sophus::SE3d m_pose_w;
    Eigen::Vector3d m_scan_origin = Eigen::Vector3d::Zero();


    const Sophus::SE3d & getMapInitPose() const;
    const Sophus::SE3d & getMapPose() const;
    size_t getNumPoints ( ) const;

    static Ptr create ( const MapParameters & params = MapParameters::getMapParameters(), const bool & allocate = true);
    static Ptr create ( const MarsMapBase & otherMap, const bool & allocate = true );

    MarsMapBase ( const MapParameters & params = MapParameters::getMapParameters(), const bool & allocate = true );

    void clear();

    void initParams ( const MarsMapBase & otherMap );
    void transform ( const Sophus::SE3d & pose_w );
    void transform_shift ( const Eigen::Vector3d & t_w );
    void setPose ( const Sophus::SE3d & pose_w );

    void setCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_w = Sophus::SE3d(), const Eigen::Vector3d & scan_origin = Eigen::Vector3d::Zero() );
    void addCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_w2 = Sophus::SE3d() );
    void removeCloudByScanId ( const IndexType & scan_id );

    void setCells ( const SurfelInfoVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void addCells ( const SurfelInfoVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void setCellsOnGrid ( const SurfelInfoVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void addCellsOnGrid ( const SurfelInfoVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void setCells ( const SurfelInfoConstPtrVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void addCells ( const SurfelInfoConstPtrVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void setCellsOnGrid ( const SurfelInfoConstPtrVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void addCellsOnGrid ( const SurfelInfoConstPtrVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate = true, const LevelIndexType & level_modifier = 0 );
    void updateCells ( );

    SurfelInfoVector iterate() ;
    int getSensorCells( const SurfelInfoVector & data ) ;

    void getCells ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center = true ) const;
    void getCellsFused ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center = true ) const;
    void getCellsScanWise ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center = true, const LevelIndexType & level_modifier = 0 ) const;
    int countSensorCell( const SurfelInfoVector & data ) const;

    bool getCell( const Eigen::Vector3f & pt_w, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors = 0) const;
    bool getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors = 0) const;

    void getCellsNoCenter ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w ) const;
    void getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center = true ) const;
    void getCells ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w ) ; //const;
    double getCellSize ( const LevelIndexType & lvl ) const;

    Eigen::Vector3i getOriginCellShift( const Sophus::SE3d & map_pose, Sophus::SE3d & new_map_pose, Sophus::SE3d & map_integration_pose, const int & lvl = 0 );

    void getCellsAdaptive ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w) const;
    void getCellsAdaptive ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w ) ;

    void update( const bool & use_adaptive = false );
};
