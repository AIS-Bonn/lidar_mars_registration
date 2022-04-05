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
#include "MarsBlockSparseVoxelGrid.h"

#include <absl/container/flat_hash_set.h>

#include <Tracy.hpp>

//#define OPT_THREAD_NUM 3
//#define USE_TBB
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

BlockSparseVoxelGrid::BlockSparseVoxelGrid( const MapParameters & params )
    : m_map_params ( params ), m_maps ( params.m_num_levels ),
      m_maps_surfel_info ( params.m_num_levels )
{

}

SurfelInfoVector BlockSparseVoxelGrid::iterate() const
{
    const bool should_omit_center = true;
    SurfelInfoVector surfels;
    surfels.clear();
    surfels.reserve(1e5);
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        const bool omit_center = should_omit_center && lvlIdx < m_map_params.m_num_levels-1;
        for ( const auto & block_volume : m_maps[lvlIdx] )
        {
            const BlockType & block = block_volume.second;
            for ( const auto it : block )
            {
                const CellType & cell = it.second;

                if ( cell.isEmpty() || ! cell.m_surfel ) continue;
                if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
                surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center);
            }
        }
    }
    return surfels;
}

void BlockSparseVoxelGrid::allocate()
{
    ZoneScopedN("BlockSparseVoxelGrid::allocate");
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        if ( map_level.empty() )
            map_level.reserve(std::min(500.,std::pow(m_map_params.getNumCells(lvlIdx)/m_map_params.block_size,3)));
        if ( m_maps_surfel_info[lvlIdx].empty() )
            m_maps_surfel_info[lvlIdx].reserve(std::min(500.,std::pow(m_map_params.getNumCells(lvlIdx),3)));
    }
}


void BlockSparseVoxelGrid::clear()
{
    ZoneScopedN("BlockSparseVoxelGrid::clear");
    for ( LevelIndexType lvl = 0; lvl < m_map_params.m_num_levels; ++lvl )
    {
        {
        ZoneScopedN("BlockSparseVoxelGrid::clear::map_lvl");
        m_maps[lvl].clear();
        }
//        {
//        ZoneScopedN("BlockSparseVoxelGrid::clear::map_lvl_erase");
//        m_maps[lvl].erase(m_maps[lvl].begin(), m_maps[lvl].end());
//        }
        {
        ZoneScopedN("BlockSparseVoxelGrid::clear::surfel_info_lvl");
        m_maps_surfel_info[lvl].clear();
        }
    }
}

int BlockSparseVoxelGrid::update( const SurfelInfoVector & surfels )
{
    int num_points = 0;
    for ( auto & map : m_maps_surfel_info )
    {
        map.clear();
        map.reserve(surfels.size()/m_map_params.block_volume);
    }
    for ( const SurfelInfo & surfel : surfels )
    {
        if ( ! surfel.m_surfel ) continue;

        const Eigen::Vector3i cell_idx_vec = Eigen::Map<const Eigen::Vector3i>(surfel.m_index.data());
        const Eigen::Vector3i block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( cell_idx_vec );

        const Eigen::Vector3i block_idx_inside_vec = m_map_params.fromCellIndexVectorToBlockInsideIdx( cell_idx_vec, block_idx_vec );
        const IndexType blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( block_idx_inside_vec );

        const CellIndexType blockIdx({block_idx_vec[0],block_idx_vec[1],block_idx_vec[2]});

        BlockSurfelType & block = m_maps_surfel_info[surfel.m_level][blockIdx];
        block[blockInsideIdx] = &surfel;
        num_points += surfel.m_surfel->getNumPoints();
    }
    return num_points;
}

void BlockSparseVoxelGrid::removeCloudByScanId( const IndexType & scan_id, int & cellsWithRemovedPoints, int & pointsRemoved, int & cellsRemoved)
{
    cellsRemoved = 0;
    pointsRemoved = 0;
    cellsWithRemovedPoints = 0;
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        int cellsWithRemoved = 0;
        MapLevelType & map_level = m_maps[lvlIdx];
        std::vector<CellIndexType> blockIndicesToBeRemoved;
        for ( auto & block_pair : map_level )
        {
            bool allEmpty = true;
            for ( auto & cellIt : block_pair.second )
            {
                CellType & cell = cellIt.second;
                if ( cell.isEmpty() ) continue;
                const int prevPoints = cell.getNumPoints();
                const bool removedSome = cell.removeByScanId ( scan_id );
                cellsWithRemoved+=removedSome;
                if ( !cell.isEmpty() )
                    allEmpty = false;
                else
                    cell.updateSurfel( true );
                const int afterPoints = cell.getNumPoints();
                pointsRemoved += prevPoints - afterPoints;
            }
            if ( allEmpty )
                blockIndicesToBeRemoved.emplace_back(block_pair.first);
        }
        int num_erased = 0;
        for ( const CellIndexType & blockRemoveIdx : blockIndicesToBeRemoved )
        {
            num_erased += map_level.erase( blockRemoveIdx );
        }
        cellsWithRemovedPoints+=cellsWithRemoved;
        cellsRemoved += num_erased;
    }
}

void BlockSparseVoxelGrid::getOriginCellShift( const Eigen::Vector3f & rel_translation, const LevelIndexType & lvl, Eigen::Vector3f & map_translation, Eigen::Vector3i & num_cells ) const
{
    const double min_cell_size = m_map_params.getCellSizeOnLevel( lvl );
    num_cells = (rel_translation.array() / min_cell_size).cast<int>();
    map_translation = num_cells.cast<float>() * min_cell_size;
}

void BlockSparseVoxelGrid::getCellsAdaptive ( SurfelInfoVector & surfels, int & normal_level_cnter, int & coarser_level_cnter ) const
{
    normal_level_cnter = 0;
    coarser_level_cnter = 0;

    static Eigen::Matrix3Xi neighborAdd = Eigen::Matrix3Xi::Zero(3,8);
    static bool neighborAdded = false;
    if ( ! neighborAdded )
    {
        constexpr int neighbors = 2; // due to subdivision by 2, there are only 2 entries per axis
        int vi = 0;
        for ( IndexType iz = 0; iz < neighbors; ++iz )
        {
            for ( IndexType iy = 0; iy < neighbors; ++iy )
            {
                for ( IndexType ix = 0; ix < neighbors; ++ix )
                {
                    neighborAdd.col(vi) << ix,iy,iz;
                    ++vi;
                }
            }
        }
        neighborAdded = true;
    }

    const auto checkPlanarity = [](const CellType & cell ){
        return ( ( cell.m_surfel->eigen_values_(0) < MapParameters::first_ev_planarity_threshold )
              || ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold )
              || ( cell.m_surfel->eigen_values_(0) < MapParameters::plane_scale_factor * cell.m_surfel->eigen_values_(1)));};
    const auto checkIsValid = [] ( const CellType & cell ) { return cell.m_surfel->valid_ && cell.m_surfel->getNumPoints() >= 10.; };
    const auto checkNotInCenter = [] ( const CellType & cell ) { return ! cell.m_surfel_info.m_in_center; };

    absl::flat_hash_set<CellIndexType> prevCheckingCoarserFailed, prevUseCoarserLevel;
    for ( LevelIndexType lvlIdx = m_map_params.m_num_levels-1; lvlIdx >= 1; --lvlIdx ) // bottom up // top down
    {
        const auto & mapLevel = m_maps[lvlIdx];
        absl::flat_hash_set<CellIndexType> useCoarserLevel, checkingCoarserFailed;
        useCoarserLevel.reserve(mapLevel.size());
        checkingCoarserFailed.reserve(mapLevel.size());
        //SurfelInfoVector coarser_surfels;
        //coarser_surfels.reserve(mapLevel.size());
        for ( const auto & pa : mapLevel )
        {
            const CellIndexType & block_index_vec = pa.first;
            const BlockType & block = pa.second;
            for ( const auto & cellIt : block )
            {
                const CellType & cell = cellIt.second;
                if ( ! cell.m_surfel ) LOG(FATAL) << "no surfel set."; // Should never happen

                const Eigen::Vector3i search_idx_vec_coarse = m_map_params.toCellIndexVector( cell.m_surfel_info.m_center_s, lvlIdx-1 );
                const CellIndexType coarse_index ({search_idx_vec_coarse[0],search_idx_vec_coarse[1], search_idx_vec_coarse[2]});

                // c = use coarser
                // m = merge
                // a = accepted
                // - = not set
                // x = failed
                //     0     |     1     |     2 =>     a     |     a     |     x
                //  0  |  1  |  2  |  3  |  4  | =>  m  |  m  |  m  |  c  |  x  |  a
                // 0 1 | 2 3 | 4 5 | 6 7 | 8 9 | => c c | c c | c - | - - | c x | c c

                if ( prevCheckingCoarserFailed.count(pa.first) > 0 ){ // left case (c x -> x) on the right
                    // this one was rejected on previous level
                    checkingCoarserFailed.insert(coarse_index); // thus the one above also fails.
                    continue;
                }

                //top-down: If this one is invalid, we do not need to check smaller levels since they will also be invalid!
                const bool isValid = checkIsValid ( cell );

                bool coarser_is_valid = true;
                bool checkCoarser = true;
                if ( isValid )
                    checkCoarser = checkPlanarity( cell );

                if ( useCoarserLevel.count(coarse_index) > 0 )  // corresponding coarser was already added and should be used, skip this one
                    continue;

                // finishes last case ( c c -> a ) on the right
                if ( checkingCoarserFailed.count(coarse_index) > 0 ) // corresponding coarser was already added and this one was checked!
                {
                    checkCoarser = false;
                }


                if ( checkCoarser ) // should check upper one
                {
                    // check neighbors:
                    const Eigen::Matrix3Xi v2 = neighborAdd.colwise() + 2 * search_idx_vec_coarse;
                    const Eigen::Vector3i block_idx_vec = Eigen::Map<const Eigen::Vector3i>(block_index_vec.data());
                    const IndexType ownBlockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx( Eigen::Map<const Eigen::Vector3i>(cell.m_surfel_info.m_index.data()), block_idx_vec ) );

                    //LOG(1) << "v2:\n" << v2.row(0) << "\n" << v2.row(1) << "\n" << v2.row(2);
                    Eigen::Vector3f fine_normal = cell.m_surfel->normal_;
                    int numChecked = 0;
                    for ( int nIdx = 0; nIdx < v2.cols(); ++nIdx)
                    {
                        const CellIndexType nI ({v2(0,nIdx),v2(1,nIdx),v2(2,nIdx)});
                        if ( prevCheckingCoarserFailed.count(nI) > 0 ){ // inverse of left case ( a x -> x ) on the right
                            // this one previously failed!
                            checkCoarser = false;
                            break;
                        }
                        const IndexType blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx( v2.col(nIdx), block_idx_vec ) );
                        if ( blockInsideIdx != ownBlockInsideIdx ) // only neighbors, not this one.
                        {
                            const auto bit = block.find(blockInsideIdx);
                            if ( bit == block.end() ) continue;
                            const auto & other_cell = bit->second;
                            if ( ! other_cell.m_surfel ) continue;
                            if ( ! checkIsValid ( other_cell ) ) continue;
                            checkCoarser &= checkPlanarity ( other_cell );
                            fine_normal+= other_cell.m_surfel->normal_;
                            ++numChecked;
                        }
                        if ( ! checkCoarser ) break; // need to add all.
                    }
                    if ( checkCoarser )
                    {
                        // check the actual coarser one
                        const Eigen::Vector3i coarser_block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( search_idx_vec_coarse );
                        const CellIndexType coarser_block_index ( {coarser_block_idx_vec[0],coarser_block_idx_vec[1],coarser_block_idx_vec[2]} );
                        const auto it = m_maps[lvlIdx-1].find(coarser_block_index);
                        if ( it != m_maps[lvlIdx-1].end() )
                        {
                            const auto & coarser_block = it->second;

                            const IndexType coarserBlockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx ( search_idx_vec_coarse, coarser_block_idx_vec ) );
                            const auto cit = coarser_block.find(coarserBlockInsideIdx);
                            if ( cit == coarser_block.end() ) continue;
                            const auto & coarser_cell = cit->second;
                            coarser_is_valid = checkIsValid ( coarser_cell );


                            if ( coarser_cell.m_surfel && coarser_is_valid && checkPlanarity ( coarser_cell )
                                 && ( numChecked == 0 || (numChecked > 0 && std::abs(fine_normal.normalized().dot(coarser_cell.m_surfel->normal_)) > 0.8) )
                                 )
                            {
                                // use coarser one.
                                useCoarserLevel.insert(coarse_index);
                            }
                            else
                                checkCoarser = false; // should add the neighbors too
                        }
                    }

                    if ( ! checkCoarser ) // upper failed, add all neighbors
                    {
                        // adding all available neighbors:
                        for ( int nIdx = 0; nIdx < v2.cols(); ++nIdx)
                        {
                            const CellIndexType nI ({v2(0,nIdx),v2(1,nIdx),v2(2,nIdx)});
                            if ( prevCheckingCoarserFailed.count(nI) > 0 ) continue; // this one failed before, thus should not be added.

                            const IndexType blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx( v2.col(nIdx), block_idx_vec ) );
                            if ( blockInsideIdx != ownBlockInsideIdx ) // only neighbors, not this one.
                            {
                                const auto bit = block.find(blockInsideIdx);
                                if ( bit == block.end() ) continue;
                                const auto & other_cell = bit->second;
                                if ( ! other_cell.m_surfel ) continue;
                                if ( ! checkIsValid ( other_cell ) ) continue;
                                if ( other_cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                                surfels.emplace_back(other_cell.m_surfel, other_cell.m_surfel_info.m_center_s, lvlIdx, other_cell.m_surfel_info.m_index, other_cell.m_cov_per_scan.empty() ? -1 : other_cell.m_cov_per_scan.front().first, other_cell.m_surfel_info.m_in_center, other_cell.m_class);
                                if ( ! other_cell.m_surfel_info.m_in_center )
                                    ++coarser_level_cnter;
                                else
                                    ++normal_level_cnter;
                            }
                        }
                    }
                }
                if ( ! checkCoarser && coarser_is_valid )
                    checkingCoarserFailed.insert(coarse_index);

                if ( !checkCoarser && isValid ) // add current one
                {
                    if ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                    surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
                    if ( ! cell.m_surfel_info.m_in_center )
                        ++coarser_level_cnter;
                    else
                        ++normal_level_cnter;
                }
            }
        }
        prevCheckingCoarserFailed = checkingCoarserFailed;
        prevUseCoarserLevel = useCoarserLevel;
    }

    int notInCenterL1 = 0;
    for( const auto & s : surfels )
        if ( s.m_level == 1 && !s.m_in_center )
            ++notInCenterL1;

    // add outer rim if available
    const auto & mapLevel = m_maps[0];
    for ( const auto & pa : mapLevel )
    {
        const BlockType & block = pa.second;
        for ( const auto & cellIt : block )
        {
            const CellType & cell = cellIt.second;
            if ( ! cell.m_surfel ) LOG(FATAL) << "nullptr on surfel?"; // Should never happen

            //top-down: If this one is invalid, we do not need to check smaller levels since they will also be invalid!
            const bool isInvalid = ! checkIsValid ( cell );
            if ( isInvalid ) continue;

            if ( checkNotInCenter(cell) || prevUseCoarserLevel.count(pa.first) > 0 )
            {
                if ( ! cell.m_surfel_info.m_in_center )
                    ++coarser_level_cnter;
                else
                    ++normal_level_cnter;

                surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, 0, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center,  cell.m_class);
            }
        }
    }

    std::vector<int> surfelsPerLevel(m_maps.size(),0);
    for ( const auto & s : surfels )
        ++surfelsPerLevel[s.m_level];
    LOG(1) << "NewAdaptive: per lvl: " << Eigen::Map<Eigen::VectorXi>(surfelsPerLevel.data(),surfelsPerLevel.size(), 1).transpose() << " prCoarser: " << prevUseCoarserLevel.size() << " notCL1: " << notInCenterL1;
}

void BlockSparseVoxelGrid::getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        for ( const auto & cellIt : pa.second )
        {
            const CellType & cell = cellIt.second;
            if ( cell.isEmpty() || ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center);
        }
    }
}

void BlockSparseVoxelGrid::getCellsFusedOnLevel ( const IndexType & lvlIdx, SurfelInfoConstPtrVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        for ( const auto & cellIt : pa.second )
        {
            const CellType & cell = cellIt.second;
            if ( cell.isEmpty() || ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            surfels.emplace_back(cell.getSurfelInfo());
        }
    }
}

void BlockSparseVoxelGrid::getCellsScanWiseOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center, const LevelIndexType & level_modifier  ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        for ( auto & cellIt : pa.second )
        {
            const CellType & cell = cellIt.second;
            if ( cell.isEmpty() || ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : cell.m_cov_per_scan)
                surfels.emplace_back(surfelInfo.second, cell.m_surfel_info.m_center_s, lvlIdx+level_modifier, cell.m_surfel_info.m_index, surfelInfo.first, cell.m_surfel_info.m_in_center, cell.m_class);
        }
    }
}

void BlockSparseVoxelGrid::updateCells()
{
    ZoneScopedN("BlockSparseVoxelGrid::updateChangedCells");
#ifdef USE_TBB
    tbb::parallel_for(tbb::blocked_range<LevelIndexType>(0,m_map_params.m_num_levels,1), [&](const tbb::blocked_range<LevelIndexType> & range)
    {
        for( LevelIndexType lvlIdx=range.begin(); lvlIdx!=range.end(); ++lvlIdx)
        {
#else
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
#endif
        ZoneScopedN("BlockSparseVoxelGrid::updateChangedCells::PerLevel");
        MapLevelType & map_level = m_maps[lvlIdx];
        for ( auto & block : map_level )
        {
            for ( auto & it : block.second )
            {
                CellType & cell = it.second;
                if ( cell.isUpdated() ) continue;
                cell.updateSurfel(true);
            }
        }
    }
#ifdef USE_TBB
    });
#endif
}

template <typename SurfelInfoT>
int BlockSparseVoxelGrid::addCells ( const std::vector<SurfelInfoT> & surfels, const Sophus::SE3d & orig_pose_s1s2, const LevelIndexType & level_modifier )
{
    ZoneScopedN("BlockSparseVoxelGrid::addCells");
    Sophus::SE3d pose_s1s2 = orig_pose_s1s2;
    if ( pose_s1s2.unit_quaternion().isApprox(Eigen::Quaterniond::Identity()) )
        pose_s1s2.so3().setQuaternion(Eigen::Quaterniond::Identity());
    if ( pose_s1s2.translation().isApprox(Eigen::Vector3d::Zero()))
        pose_s1s2.translation().setZero();

    int num_points = 0;
    Sophus::SE3d n_pose_s1s2 = pose_s1s2;
    const size_t num_surfels = surfels.size();
    for ( size_t idx = 0; idx < num_surfels; ++idx)
    {
        if constexpr ( std::is_pointer<SurfelInfoT>::value )
            num_points += addCellSurfel ( *(surfels[idx]), pose_s1s2, level_modifier, n_pose_s1s2);
        else
            num_points += addCellSurfel ( surfels[idx], pose_s1s2, level_modifier, n_pose_s1s2);
    }
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(1) << "num_points: " << num_points << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose();
    return num_points;
}

inline int BlockSparseVoxelGrid::addCellSurfel ( const SurfelInfo & surfel_info, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier, Sophus::SE3d & n_pose_s1s2 )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;
    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & center_s = surfel_info.m_center_s;
    const Eigen::Vector3f meanPt = (pose_s1s2 * (surfel->mean_ + center_s).template cast<double>()).template cast<float>();
    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector ( meanPt, lvlIdx );
    const Eigen::Vector3i block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( cellIdxVec );
    const IndexType blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx( cellIdxVec, block_idx_vec ) );

    const CellIndexType blockIdx ( {block_idx_vec[0],block_idx_vec[1],block_idx_vec[2]} );
    BlockType & block = map_level[blockIdx];
    CellType & cell = block[blockInsideIdx];

    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > std::sqrt(std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            return 0;
        }
        cell.create( new_cell_center, (pose_s1s2.translation().template cast<float>()-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, CellIndexType({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]}) );
    }

    // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
    // mean_new = meanPt - center_new - surfel mean

    n_pose_s1s2.translation() =  (meanPt - surfel->mean_ - cell.m_surfel_info.m_center_s).template cast<double>();

    cell.addSurfel(surfel, n_pose_s1s2, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}


inline int BlockSparseVoxelGrid::addCellSurfelOnGrid( const SurfelInfo & surfel_info, const Eigen::Vector3f & t_s1s2, const LevelIndexType & level_modifier, const Sophus::SE3d & id )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;
    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & old_center = surfel_info.m_center_s;
    const Eigen::Vector3f surfel_mean = surfel->mean_ + old_center;
    const Eigen::Vector3f meanPt = surfel_mean + t_s1s2;
    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector ( meanPt, lvlIdx );
    const Eigen::Vector3i block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( cellIdxVec );
    const IndexType blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( m_map_params.fromCellIndexVectorToBlockInsideIdx( cellIdxVec, block_idx_vec ) );

    const CellIndexType blockIdx ( {block_idx_vec[0],block_idx_vec[1],block_idx_vec[2]} );
    BlockType & block = map_level[blockIdx];
    CellType & cell = block[blockInsideIdx];

    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > std::sqrt(std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            return 0;
        }
        cell.create( new_cell_center,  (t_s1s2-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, CellIndexType({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]}) );
    }
    // -> offset does not matter, since the shift is on the grid!
    cell.addSurfel(surfel, id, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}

template <typename SurfelInfoT>
int BlockSparseVoxelGrid::addCellsOnGrid ( const std::vector<SurfelInfoT> & surfels, const Eigen::Vector3d & orig_t_s1s2, const LevelIndexType & level_modifier )
{
    // Assumes surfels are on the same grid with a shift of n_i per dimension. -> no rotation or fractions within translation
    ZoneScopedN("BlockSparseVoxelGrid::addCells");

    Eigen::Vector3f t_s1s2 = Eigen::round(orig_t_s1s2.template cast<float>().array());
    const Sophus::SE3d id;
    int num_points = 0;
    const size_t num_surfels = surfels.size();
    for ( size_t idx = 0; idx < num_surfels; ++idx)
    {
        if constexpr ( std::is_pointer<SurfelInfoT>::value )
            num_points += addCellSurfelOnGrid ( *(surfels[idx]), t_s1s2, level_modifier, id);
        else
            num_points += addCellSurfelOnGrid ( surfels[idx], t_s1s2, level_modifier, id);
    }
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(1) << "num_points: " << num_points << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose();
    return num_points;
}

template int BlockSparseVoxelGrid::addCellsOnGrid<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Eigen::Vector3d & , const LevelIndexType &);
template int BlockSparseVoxelGrid::addCellsOnGrid<SurfelInfo>( const SurfelInfoVector & , const Eigen::Vector3d & , const LevelIndexType & );
template int BlockSparseVoxelGrid::addCells<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Sophus::SE3d & , const LevelIndexType & );
template int BlockSparseVoxelGrid::addCells<SurfelInfo>( const SurfelInfoVector & , const Sophus::SE3d & , const LevelIndexType & );

int BlockSparseVoxelGrid::addCloud( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_s1s2, const Eigen::Vector3d & origin )
{
    ZoneScopedN("BlockSparseVoxelGrid::addCloudBlockSparse");
    const Eigen::Vector3f originf = origin.cast<float>();
    const Eigen::Matrix3Xf pts_s = (pose_s1s2.so3().cast<float>().matrix() * cloud->m_points).colwise() + pose_s1s2.translation().cast<float>();
    const int scan_id = cloud->single_scan() ? cloud->m_scan_id(0) : -1; // TODO
    const int num_pts = cloud->size();


    LOG(1) << "scan id: " << scan_id << " np: " << cloud->size() << " pts: " << cloud->m_points.cols() << " " << pts_s.cols();
    absl::FixedArray<absl::flat_hash_set<CellIndexType>,MapParameters::maxAllowedMapLevels> changedBlocks ( m_map_params.m_num_levels );
    int num_points = 0;
#ifdef USE_TBB
    std::vector<int> num_points_v( m_map_params.m_num_levels, 0);
    tbb::parallel_for(tbb::blocked_range<LevelIndexType>(0,m_map_params.m_num_levels,1), [&](const tbb::blocked_range<LevelIndexType> & range)
    {
        for( LevelIndexType lvlIdx=range.begin(); lvlIdx!=range.end(); ++lvlIdx)
        {
            int & num_points = num_points_v[lvlIdx];
#else
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:num_points)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
#endif
        ZoneScopedN("BlockSparseVoxelGrid::addCloudBlockSparse::PerLevel");
        MapLevelType & map_level = m_maps[lvlIdx];
        absl::flat_hash_set<CellIndexType> & changed_level_blocks = changedBlocks[lvlIdx];
        IndexType blockInsideIdx = 0;
        for ( int pointIdx = 0; pointIdx < num_pts; ++pointIdx)
        {
            const Eigen::Vector3f pt_s = pts_s.col(pointIdx);
            const bool isInside = m_map_params.isInBounds( pt_s, lvlIdx );
            if ( ! isInside ) continue;
            // add data.

            const Eigen::Vector3i cell_idx_vec = m_map_params.toCellIndexVector ( pt_s, lvlIdx);
            const CellIndexType cellIdx ({cell_idx_vec[0],cell_idx_vec[1],cell_idx_vec[2]});

            const Eigen::Vector3i block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( cell_idx_vec );
            const Eigen::Vector3i block_idx_inside_vec = m_map_params.fromCellIndexVectorToBlockInsideIdx( cell_idx_vec, block_idx_vec );
            blockInsideIdx = m_map_params.fromBlockIndexInsideVectorToBlockIdx ( block_idx_inside_vec );

            const CellIndexType blockIdx({block_idx_vec[0],block_idx_vec[1],block_idx_vec[2]});

            BlockType & block = map_level[blockIdx];

            CellType & cell = block[blockInsideIdx];
            if ( cell.isEmpty() )
            {
                const Eigen::Vector3f center_s = m_map_params.fromCellIdxToCenter<float> ( cell_idx_vec, lvlIdx);
                cell.create( center_s, (pt_s-originf).normalized(), m_map_params.isInBounds<float>( center_s, lvlIdx+1 ), lvlIdx, cellIdx );
            }

            changed_level_blocks.insert(blockIdx);

            if ( cell.isUpdated() )
            {
                cell.addNewScan( scan_id, (pt_s-originf).normalized() );
            }

            cell.addPoint( pt_s );
            //#pragma omp atomic
            ++num_points;
        }
        LOG(1) << "lvl: " << lvlIdx << " addedPts: " << num_points << " clb: " << changed_level_blocks.size();
    }
#ifdef USE_TBB
    });
    for ( const auto & np : num_points_v )
        num_points+=np;
#endif

    int updated = 0;
    {
        ZoneScopedN("BlockSparseVoxelGrid::addCloud::UpdateSurfel");
#ifdef USE_TBB
        std::vector<int> updated_v(changedBlocks.size(),0);
    tbb::parallel_for(tbb::blocked_range<size_t>(0,changedBlocks.size(),1), [&](const tbb::blocked_range<size_t> & range)
    {
        for( size_t lvlIdx=range.begin(); lvlIdx!=range.end(); ++lvlIdx)
        {
            int & updated = updated_v[lvlIdx];
#else
        #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:updated)
        for ( size_t lvlIdx = 0; lvlIdx < changedBlocks.size(); ++lvlIdx )
        {
#endif
            ZoneScopedN("BlockSparseVoxelGrid::addCloud::UpdateSurfel::PerLevel");
            const absl::flat_hash_set<CellIndexType> & changed_level_blocks = changedBlocks[lvlIdx];
            auto & map_level = m_maps[lvlIdx];
            for ( const auto & blockIdx : changed_level_blocks )
            {
                const auto & blockIt = map_level.find(blockIdx);
                if ( blockIt != map_level.end() )
                {
                    BlockType & block = blockIt->second;
                    for ( auto & cellIt : block )
                    {
                        CellType & cell = cellIt.second;
                        if ( cell.isUpdated() ) continue;
                        cell.updateSurfel( );
                        ++updated;
                    }
                }
            }
        }
#ifdef USE_TBB
        });
        for ( const auto & up : updated_v ) updated+=up;
#endif
    }
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsInclFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        for ( const auto & c : m_maps[lvlIdx] )
            cellsInclFilled[lvlIdx] += c.second.size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(1) << "num_points: " << num_points << " c: " << cloud->size() << " updatedCells: " << updated << " cellsPerLevel: " << cellsFilled.transpose() << " ( " << cellsInclFilled.transpose() << " ) max: " << cellsMaxFilled.transpose();
    if ( updated == 0 ) LOG(FATAL) << "wtf?";
    return num_points;
}

bool BlockSparseVoxelGrid::getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors ) const
{
    //ZoneScopedN("BlockSparseVoxelGrid::getSensorCell");
    cellPtrs.clear();
    if ( search_lvl < 0 || search_lvl >= LevelIndexType(m_maps.size()) ) { LOG(1) << "ooL? " << pt_s.transpose() << " lvl: " << search_lvl; return false; } // check lvl bounds
    if ( ! m_map_params.isInBounds ( pt_s , search_lvl ) ) {  LOG(1) << "oob: " << pt_s.transpose() << " lvl: " << search_lvl; return false; }

    cellPtrs.reserve(27);

    if ( neighbors > 1 ) LOG(FATAL) <<"Neighbors only implemented for 1 in every direction.";
    if ( m_map_params.block_size <= 1 ) LOG(FATAL) << "block size should be larger than: " << m_map_params.block_size;
    if ( m_map_params.block_size <= 1 ) LOG(FATAL) << "block size should be larger than: " << m_map_params.block_size;

    constexpr int min_cell = -(m_map_params.NUM_CELLS>>1);
    constexpr int max_cell = (m_map_params.NUM_CELLS>>1)-1;
    constexpr int num_blocks_per_side_x = 1;
    constexpr int num_blocks_per_side_y = m_map_params.block_size;
    constexpr int num_blocks_per_side_z = m_map_params.block_size * m_map_params.block_size;
    constexpr int block_size_m1 = m_map_params.block_size-1;

    const auto & map_level = m_maps_surfel_info[search_lvl];

    const Eigen::Vector3i orig_search_idx_vec = m_map_params.toCellIndexVector( pt_s, search_lvl );

    const int orig_neg_neighbors_z = orig_search_idx_vec.z()-neighbors  < min_cell ? 0 : -neighbors;
    const int orig_neg_neighbors_y = orig_search_idx_vec.y()-neighbors  < min_cell ? 0 : -neighbors;
    const int orig_neg_neighbors_x = orig_search_idx_vec.x()-neighbors  < min_cell ? 0 : -neighbors;
    const int orig_pos_neighbors_z = orig_search_idx_vec.z()+neighbors >= max_cell ? 0 : +neighbors;
    const int orig_pos_neighbors_y = orig_search_idx_vec.y()+neighbors >= max_cell ? 0 : +neighbors;
    const int orig_pos_neighbors_x = orig_search_idx_vec.x()+neighbors >= max_cell ? 0 : +neighbors;

    const int len_pos_neighbors_z = orig_pos_neighbors_z-orig_neg_neighbors_z; // -neg will become positive
    const int len_pos_neighbors_y = orig_pos_neighbors_y-orig_neg_neighbors_y;
    const int len_pos_neighbors_x = orig_pos_neighbors_x-orig_neg_neighbors_x;

    const Eigen::Vector3i search_idx_vec = orig_search_idx_vec + Eigen::Vector3i( orig_neg_neighbors_x, orig_neg_neighbors_y, orig_neg_neighbors_z ); // shifted to first pos in window
    const Eigen::Vector3i block_idx_vec = m_map_params.fromCellIndexVectorToBlockIdx( search_idx_vec );
    const Eigen::Vector3i block_idx_inside_vec = m_map_params.fromCellIndexVectorToBlockInsideIdx ( search_idx_vec, block_idx_vec );

    const int block_idx_inside_vec_pos_neighbors_z = block_idx_inside_vec.z()+len_pos_neighbors_z;
    const int block_idx_inside_vec_pos_neighbors_y = block_idx_inside_vec.y()+len_pos_neighbors_y;
    const int block_idx_inside_vec_pos_neighbors_x = block_idx_inside_vec.x()+len_pos_neighbors_x;

    const int pos_outside_z = block_idx_inside_vec_pos_neighbors_z >= m_map_params.block_size ? +1 : 0 ;
    const int pos_outside_y = block_idx_inside_vec_pos_neighbors_y >= m_map_params.block_size ? +1 : 0 ;
    const int pos_outside_x = block_idx_inside_vec_pos_neighbors_x >= m_map_params.block_size ? +1 : 0 ;

    int pos_inside_neighbors_z = 0;
    int pos_inside_neighbors_y = 0;
    int pos_inside_neighbors_x = 0;
    int new_block_inside_idx_z = 0;
    int new_block_inside_idx_y = 0;
    int new_block_inside_idx_x = 0;

    for ( IndexType iz = 0; iz <= pos_outside_z; ++iz )
    {
        const int bz = block_idx_vec.z()+iz;

        if ( iz == 0 ) // same block
        {
            new_block_inside_idx_z = block_idx_inside_vec.z();
            pos_inside_neighbors_z = std::min(block_size_m1-new_block_inside_idx_z,len_pos_neighbors_z) ;
        }
        else
        {
            new_block_inside_idx_z = 0;
            pos_inside_neighbors_z = block_idx_inside_vec_pos_neighbors_z - m_map_params.block_size; // left outmost point
        }

        //LOG(1) << "nz: " << pos_outside_z << " iz: " << iz << " nbi: " << new_block_inside_idx_z << " pin: " << pos_inside_neighbors_z;

        for ( IndexType iy = 0; iy <= pos_outside_y; ++iy )
        {
            if ( iy == 0 ) // same block
            {
                new_block_inside_idx_y = block_idx_inside_vec.y();
                pos_inside_neighbors_y = std::min(block_size_m1-new_block_inside_idx_y,len_pos_neighbors_y);
            }
            else
            {
                new_block_inside_idx_y = 0;
                pos_inside_neighbors_y = block_idx_inside_vec_pos_neighbors_y - m_map_params.block_size; // left outmost point
            }

            //LOG(1) << "ny: " << pos_outside_y << " iy: " << iy << " nbi: " << new_block_inside_idx_y  << " pin: " << pos_inside_neighbors_y;

            const int by = block_idx_vec.y()+iy;
            for ( IndexType ix = 0; ix <= pos_outside_x; ++ix )
            {
                const int bx = block_idx_vec.x()+ix;
                const CellIndexType blockIdx ( {bx,by,bz} );

                // work on block:
                const auto it = map_level.find( blockIdx );
                if ( it != map_level.end() )
                {
                    const BlockSurfelType & block = it->second;

                    if ( ix == 0 ) // same block
                    {
                        new_block_inside_idx_x = block_idx_inside_vec.x();
                        pos_inside_neighbors_x = std::min(block_size_m1-new_block_inside_idx_x,len_pos_neighbors_x);
                    }
                    else
                    {
                        new_block_inside_idx_x = 0; // left outmost point
                        pos_inside_neighbors_x = block_idx_inside_vec_pos_neighbors_x - m_map_params.block_size;
                    }

                    //LOG(1) << "nx: " << pos_outside_x << " ix: " << ix << " nbi: " << new_block_inside_idx_x << " pin: " << pos_inside_neighbors_x << " nbiv: " << new_block_idx_vec.transpose() << " nbii: " << Eigen::Vector3i(new_block_inside_idx_x,new_block_inside_idx_y,new_block_inside_idx_z).transpose();

                    for ( IndexType iiz = 0; iiz <= pos_inside_neighbors_z; ++iiz )
                    {
                        const IndexType viiz = (new_block_inside_idx_z+iiz)*num_blocks_per_side_z;
                        for ( IndexType iiy = 0; iiy <= pos_inside_neighbors_y; ++iiy )
                        {
                            const IndexType viiyz = viiz + (new_block_inside_idx_y+iiy) * num_blocks_per_side_y;
                            for ( IndexType iix = 0; iix <= pos_inside_neighbors_x; ++iix )
                            {
                                //LOG(1) << " nbii: " << Eigen::Vector3i(new_block_inside_idx_x+iix,new_block_inside_idx_y+iiy,new_block_inside_idx_z+iiz).transpose();
                                const IndexType blockInsideIdx = viiyz + (new_block_inside_idx_x+iix) * num_blocks_per_side_x;
                                const auto cit = block.find(blockInsideIdx);
                                if ( cit == block.end() ) continue;
                                const auto & c = cit->second;
                                {
                                    cellPtrs.emplace_back(c);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return !cellPtrs.empty();
}


int BlockSparseVoxelGrid::getSensorCells( const SurfelInfoVector & data, const IndexType & neighbors ) const
{
    ZoneScopedN("BlockSparseVoxelGrid::getSensorCells");
    int numNeighbors = 0;
    std::vector<SurfelInfoConstPtrVector> cnPtrs (data.size());

    const size_t num_cells = data.size();
    for ( size_t data_idx = 0; data_idx < num_cells; ++data_idx )
    {
        const LevelIndexType & search_lvl = data[data_idx].m_level;
        const Eigen::Vector3f pt_s = data[data_idx].m_surfel->mean_ + data[data_idx].m_center_s;
        SurfelInfoConstPtrVector & cellPtrs = cnPtrs[data_idx];

        getSensorCell(pt_s, search_lvl, cellPtrs, neighbors );
        numNeighbors+= cellPtrs.size();
    }

    return numNeighbors;
}
