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

#include <absl/container/fixed_array.h>
#include <Eigen/Geometry>
#include <array>

#include "MarsPointTypes.h"
#include "MarsMapParameters.h"
#include "MarsFwdDec.h"
#include "MarsSurfel.h"
#include "MarsCell.h"

class DenseVoxelGrid {

public:
    typedef float Scalar;
    typedef std::shared_ptr<DenseVoxelGrid> Ptr;

    typedef Cell CellType;
    typedef SurfelInfo::SurfelPtr SurfelPtr;
    typedef typename MapParameters::IndexType IndexType;
    typedef typename MapParameters::LevelIndexType LevelIndexType;
    typedef typename MapParameters::CellIndexType CellIndexType;
    typedef std::vector<CellType> MapLevelType;
    typedef absl::FixedArray<MapLevelType,MapParameters::maxAllowedMapLevels> MapType;

    typedef std::vector<SurfelInfoConstPtr> MapLevelSurfelType;
    absl::FixedArray<MapLevelSurfelType,MapParameters::maxAllowedMapLevels> m_maps_surfel_info;

    MapType m_maps;
    const MapParameters m_map_params;

    static Ptr create( const MapParameters & params )
    {
        return std::make_shared<DenseVoxelGrid>( params );
    }

    int addCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_s1s2, const Eigen::Vector3d & origin );

    inline int addCellSurfel ( const SurfelInfo & surfel_info, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier, Sophus::SE3d & n_pose_s1s2 );
    inline int addCellSurfelOnGrid( const SurfelInfo & surfel_info, const Eigen::Vector3f & t_s1s2, const LevelIndexType & level_modifier, const Sophus::SE3d & id );
    template <typename SurfelInfoT>
    int addCells ( const std::vector<SurfelInfoT> & surfels, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier );
    template <typename SurfelInfoT>
    int addCellsOnGrid ( const std::vector<SurfelInfoT> & surfels, const Eigen::Vector3d & orig_t_s1s2, const LevelIndexType & level_modifier );

    void clear();

    void allocate();
    SurfelInfoVector iterate() const;

    int update( const SurfelInfoVector & surfels );// updates surfel info
    void updateCells();

    void getOriginCellShift( const Eigen::Vector3f & rel_translation, const LevelIndexType & lvl, Eigen::Vector3f & map_translation, Eigen::Vector3i & num_cells ) const;

    void removeCloudByScanId( const IndexType & scan_id, int & cellsWithRemovedPoints, int & pointsRemoved, int & cellsRemoved);
    void getCellsAdaptive( SurfelInfoVector & surfels, int & normal_level_cnter, int & coarser_level_cnter ) const;
    void getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center = true ) const;
    void getCellsScanWiseOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center = true, const LevelIndexType & level_modifier = 0 ) const;
    void getCellsFusedOnLevel ( const IndexType & lvlIdx, SurfelInfoConstPtrVector & surfels, const bool & omit_center ) const;

    int getSensorCells( const SurfelInfoVector & data, const IndexType & neighbors = 0 ) const;
    bool getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors = 0) const;

public:
    DenseVoxelGrid( const MapParameters & params );
    ~DenseVoxelGrid(){}
};
