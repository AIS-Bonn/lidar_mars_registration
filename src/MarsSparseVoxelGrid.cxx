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

#include "MarsSparseVoxelGrid.h"
#include <absl/container/flat_hash_set.h>
#include <Tracy.hpp>

//#define OPT_THREAD_NUM 3

SparseVoxelGrid::SparseVoxelGrid( const MapParameters & params )
    : m_map_params ( params ), m_maps ( params.m_num_levels ), m_maps_surfel_info ( params.m_num_levels )
{
    m_last_capacity = Eigen::VectorXi::Zero(params.m_num_levels,1);
    const Eigen::VectorXd d = getDistancesBetweenVertices();
    LOG(1) << "dist between vertices: " << d.transpose();
}


Eigen::VectorXd SparseVoxelGrid::getDistancesBetweenVertices() const
{
    Eigen::VectorXd dists = Eigen::VectorXd::Zero(m_map_params.m_num_levels,1);
    dists.setConstant(std::numeric_limits<double>::max());
    constexpr int neighbors = 1;
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        const Eigen::Vector3f pt_c = m_map_params.fromCellIdxToCenter<float> (Eigen::Vector3i(0,0,0), lvlIdx);
        for ( int z = -neighbors; z <= neighbors; ++z )
            for ( int y = -neighbors; y <= neighbors; ++y )
                for ( int x = -neighbors; x <= neighbors; ++x )
                {
                    if ( (x == y) && (y == z) && x == 0 ) continue;
                    const Eigen::Vector3f pt_n = m_map_params.fromCellIdxToCenter<float>(Eigen::Vector3i(x,y,z), lvlIdx);
                    const double d_min = (pt_n-pt_c).norm();
                    if ( d_min < dists[lvlIdx] )
                        dists[lvlIdx] = d_min;
                }
    }
    return dists;
}

void SparseVoxelGrid::allocate()
{
    ZoneScopedN("SparseVoxelGrid::allocate");
    constexpr int minCellsCoarse = 512;
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        m_last_capacity[lvlIdx] = minCellsCoarse * std::pow(2,lvlIdx) - 1;
        if ( map_level.empty() )
            map_level.reserve(m_last_capacity[lvlIdx]/2);
        if ( m_maps_surfel_info[lvlIdx].empty() )
            m_maps_surfel_info[lvlIdx].reserve(m_last_capacity[lvlIdx]/2);
    }
}

SurfelInfoVector SparseVoxelGrid::iterate() const
{
    const bool should_omit_center = true;
    SurfelInfoVector surfels;
    surfels.clear();
    surfels.reserve(1e5);
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        const bool omit_center = should_omit_center && lvlIdx < m_map_params.m_num_levels-1;
        for ( const auto & pa : m_maps[lvlIdx] )
        {
            const CellType & cell = pa.second;
            if ( ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
        }
    }
    return surfels;
}

void SparseVoxelGrid::clear()
{
    for ( LevelIndexType lvl = 0; lvl < LevelIndexType(m_maps.size()); ++lvl )
    {
        m_last_capacity[lvl] = m_maps[lvl].capacity();
        m_maps[lvl].clear();
        m_maps_surfel_info[lvl].clear();
    }
}

int SparseVoxelGrid::update( const SurfelInfoVector & surfels )
{
    int num_points = 0, num_classes = 0;
    for ( auto & map : m_maps_surfel_info )
    {
        map.clear();
        map.reserve(surfels.size());
    }
    for ( const SurfelInfo & surfel : surfels )
    {
        if ( ! surfel.m_surfel ) continue;

        m_maps_surfel_info[surfel.m_level][surfel.m_index] = &surfel;
        num_points += surfel.m_surfel->getNumPoints();
        num_classes += surfel.m_class != nullptr;
    }
    LOG(1) << "update of Grid - num_clases: " << num_classes;
    return num_points;
}

void SparseVoxelGrid::getOriginCellShift( const Eigen::Vector3f & rel_translation, const LevelIndexType & lvl, Eigen::Vector3f & map_translation, Eigen::Vector3i & num_cells ) const
{
    const double min_cell_size = m_map_params.getCellSizeOnLevel( lvl );
    num_cells = (rel_translation.array() / min_cell_size).cast<int>();
    map_translation = num_cells.cast<float>() * min_cell_size;
}

void SparseVoxelGrid::removeCloudByScanId( const IndexType & scan_id, int & cellsWithRemovedPoints, int & pointsRemoved, int & cellsRemoved)
{
    cellsRemoved = 0;
    pointsRemoved = 0;
    cellsWithRemovedPoints = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:cellsWithRemovedPoints) reduction(+:pointsRemoved) reduction(+:cellsRemoved)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        int cellsWithRemoved = 0;
        MapLevelType & map_level = m_maps[lvlIdx];
        std::vector<CellIndexType> cellIndicesToBeRemoved;
        for ( auto & cell_pair : map_level )
        {
            CellType & cell = cell_pair.second;
            const int prevPoints = cell.getNumPoints();
            const bool removedSome = cell.removeByScanId ( scan_id );
            cellsWithRemoved+=removedSome;
            if ( cell.isEmpty() )
                cellIndicesToBeRemoved.emplace_back(cell_pair.first);
            else
                cell.updateSurfel( true );
            const int afterPoints = cell.getNumPoints();
            pointsRemoved += (prevPoints-afterPoints);
        }
        int num_erased = 0;
        for ( const CellIndexType & cellRemoveIdx : cellIndicesToBeRemoved )
            num_erased += map_level.erase( cellRemoveIdx );

        cellsWithRemovedPoints+=cellsWithRemoved;
        cellsRemoved += num_erased;
    }
}
void SparseVoxelGrid::getCellsAdaptive( SurfelInfoVector & surfels, int & normal_level_cnter, int & coarser_level_cnter ) const
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
        for ( const auto & pa : mapLevel )
        {
            const CellType & cell = pa.second;
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

                    if ( nI != pa.first ) // only neighbors, not this one.
                    {
                        const auto it = mapLevel.find(nI);
                        if ( it == mapLevel.end() ) continue;
                        const auto & other_cell = it->second;
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
                    const auto it = m_maps[lvlIdx-1].find(coarse_index);
                    if ( it != m_maps[lvlIdx-1].end() )
                    {
                        const auto & coarser_cell = it->second;
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
                        if ( nI != pa.first ) // only neighbors, not this one.
                        {
                            const auto it = mapLevel.find(nI);
                            if ( it == mapLevel.end() ) continue;
                            const auto & other_cell = it->second;
                            if ( ! other_cell.m_surfel ) continue;
                            if ( ! checkIsValid ( other_cell ) ) continue;
                            if ( other_cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                            surfels.emplace_back(other_cell.m_surfel, other_cell.m_surfel_info.m_center_s, lvlIdx, pa.first, other_cell.m_cov_per_scan.empty() ? -1 : other_cell.m_cov_per_scan.front().first, other_cell.m_surfel_info.m_in_center, other_cell.m_class);
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
                surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
                if ( ! cell.m_surfel_info.m_in_center )
                    ++coarser_level_cnter;
                else
                    ++normal_level_cnter;
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
        const CellType & cell = pa.second;
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

            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, 0, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center,  cell.m_class);
            //continue; // top down, outside can be added directly
        }
    }

    std::vector<int> surfelsPerLevel(m_maps.size(),0);
    for ( const auto & s : surfels )
        ++surfelsPerLevel[s.m_level];
    LOG(1) << "Adaptive: per lvl: " << Eigen::Map<Eigen::VectorXi>(surfelsPerLevel.data(),surfelsPerLevel.size(), 1).transpose() << " prCoarser: " << prevUseCoarserLevel.size() << " notCL1: " << notInCenterL1;
}

void SparseVoxelGrid::getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        if ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
        surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

void SparseVoxelGrid::getCellsFusedOnLevel ( const IndexType & lvlIdx, SurfelInfoConstPtrVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        surfels.emplace_back(cell.getSurfelInfo());
    }
}

void SparseVoxelGrid::getCellsScanWiseOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center, const LevelIndexType & level_modifier  ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : cell.m_cov_per_scan)
            surfels.emplace_back(surfelInfo.second, cell.m_surfel_info.m_center_s, lvlIdx+level_modifier, pa.first, surfelInfo.first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

void SparseVoxelGrid::updateCells()
{
    ZoneScopedN("SparseVoxelGrid::updateChangedCells");
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        for ( auto & cell : map_level )
        {
            cell.second.updateSurfel(true);
        }
    }
}

template <typename SurfelInfoT>
int SparseVoxelGrid::addCells ( const std::vector<SurfelInfoT> & surfels, const Sophus::SE3d & orig_pose_s1s2, const LevelIndexType & level_modifier )
{
    ZoneScopedN("SparseVoxelGrid::addCells");
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
    LOG(1) << "num_points: " << num_points << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose(); // << " cp: " << cp;
    return num_points;
}

inline int SparseVoxelGrid::addCellSurfel ( const SurfelInfo & surfel_info, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier, Sophus::SE3d & n_pose_s1s2 )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;
    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & center_s = surfel_info.m_center_s;
    const Eigen::Vector3f meanPt = (pose_s1s2 * (surfel->mean_ + center_s).template cast<double>()).template cast<float>();
    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector ( meanPt, lvlIdx );
    const CellIndexType cellIdx ( {cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]} );
    CellType & cell = map_level[cellIdx];
    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            LOG(FATAL) << "ncc: " << new_cell_center.transpose() <<" mp: " << meanPt.transpose() << " ct: " << center_s.transpose() << " cs: " << m_map_params.getCellSizeOnLevel(lvlIdx) << " d2: " << (new_cell_center-meanPt).squaredNorm() << " e: " << (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3);
            //continue;
        }
        // add this one.
        cell.create( new_cell_center, (pose_s1s2.translation().template cast<float>()-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, cellIdx );
    }

    // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
    // mean_new = meanPt - center_new - surfel mean

    n_pose_s1s2.translation() =  (meanPt - surfel->mean_ - cell.m_surfel_info.m_center_s).template cast<double>();

    cell.addSurfel(surfel, n_pose_s1s2, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}


inline int SparseVoxelGrid::addCellSurfelOnGrid( const SurfelInfo & surfel_info, const Eigen::Vector3f & t_s1s2, const LevelIndexType & level_modifier, const Sophus::SE3d & id )
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
    const CellIndexType cellIdx ( {cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]} );
    CellType & cell = map_level[cellIdx];
    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            LOG(FATAL) << "ncc: " << new_cell_center.transpose() <<" mp: " << meanPt.transpose() << " ct: " << old_center.transpose() << " cs: " << m_map_params.getCellSizeOnLevel(lvlIdx) << " d2: " << (new_cell_center-meanPt).squaredNorm() << " e: " << (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3);
            //continue;
        }
        // add this one.
        cell.create( new_cell_center, (t_s1s2-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, cellIdx );
    }

    // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
    // mean_new = meanPt - center_new - surfel mean

    // -> this does not matter, since the shift is on the grid!

    cell.addSurfel(surfel, id, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}

template <typename SurfelInfoT>
int SparseVoxelGrid::addCellsOnGrid ( const std::vector<SurfelInfoT> & surfels, const Eigen::Vector3d & orig_t_s1s2, const LevelIndexType & level_modifier )
{
    // Assumes surfels are on the same grid with a shift of n_i per dimension. -> no rotation or fractions within translation
    ZoneScopedN("SparseVoxelGrid::addCells");

    Eigen::Vector3f t_s1s2 = //orig_t_s1s2.template cast<float>();
                             Eigen::round(orig_t_s1s2.template cast<float>().array());
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


template int SparseVoxelGrid::addCellsOnGrid<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Eigen::Vector3d & , const LevelIndexType &);
template int SparseVoxelGrid::addCellsOnGrid<SurfelInfo>( const SurfelInfoVector & , const Eigen::Vector3d & , const LevelIndexType & );
template int SparseVoxelGrid::addCells<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Sophus::SE3d & , const LevelIndexType & );
template int SparseVoxelGrid::addCells<SurfelInfo>( const SurfelInfoVector & , const Sophus::SE3d & , const LevelIndexType & );

int SparseVoxelGrid::addCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_s1s2, const Eigen::Vector3d & origin )
{
    ZoneScopedN("SparsVoxelGrid::addCloud");
    absl::FixedArray<absl::flat_hash_set<CellIndexType>,MapParameters::maxAllowedMapLevels> changedCells ( m_map_params.m_num_levels );
    const Eigen::Vector3f originf = origin.cast<float>();
    const Eigen::Matrix3Xf pts_s = (pose_s1s2.so3().cast<float>().matrix() * cloud->m_points).colwise() + pose_s1s2.translation().cast<float>();
    const int scan_id = cloud->single_scan() ? cloud->m_scan_id(0) : -1; // TODO
    const int num_pts = cloud->size();

    int num_points = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:num_points)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        if ( map_level.capacity() < m_last_capacity[lvlIdx] )
            map_level.reserve(m_last_capacity[lvlIdx] / 2);
        absl::flat_hash_set<CellIndexType> & changed_level_cells = changedCells[lvlIdx];
        changed_level_cells.reserve(m_last_capacity[lvlIdx] / 2);

        for ( int pointIdx = 0; pointIdx < num_pts; ++pointIdx)
        {
            const Eigen::Vector3f pt_s = pts_s.col(pointIdx);
            const bool isInside = m_map_params.isInBounds( pt_s, lvlIdx );
            if ( ! isInside ) continue;
            // add data.
            //const IndexType cellIdx = m_map_params.toCellIndex ( pt_s, lvlIdx );
            const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector(pt_s,lvlIdx);
            const CellIndexType cellIdx ({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]});
            CellType & cell = map_level[cellIdx];
            if ( cell.isEmpty() )
            {
                // add this one.
                const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
                cell.create( new_cell_center, (pt_s-originf).normalized(), m_map_params.isInBounds<float>( new_cell_center, lvlIdx+1 ), lvlIdx, cellIdx );
            }
            if ( !changed_level_cells.contains(cellIdx) )
            {
                changed_level_cells.insert(cellIdx);
                cell.addNewScan( scan_id, (pt_s-originf).normalized() );
            }
            cell.addPoint( pt_s );
            ++num_points;
        }
    }

    int updated = 0, num_classes = 0;
    {
        ZoneScopedN("SparseVoxelGrid::addCloud::UpdateSurfel");
        for ( size_t lvlIdx = 0; lvlIdx < changedCells.size(); ++lvlIdx )
        {
            const absl::flat_hash_set<CellIndexType> & changed_level_cells = changedCells[lvlIdx];
            auto & map_level = m_maps[lvlIdx];
            for ( const CellIndexType & cellIdx : changed_level_cells )
            {
                const auto & cellIt = map_level.find(cellIdx);
                if ( cellIt != map_level.end() )
                {
                    CellType & cell = cellIt->second;
                    cell.updateSurfel( );
                    ++updated;
                    num_classes += ( cell.m_class != nullptr );
                }
            }
        }
    }
    //LOG(1) << "num_points: " << num_points << " updatedCells: " << updated << " numClasses: " << num_classes;
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(1) << "num_points: " << num_points << " updatedCells: " << updated << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose();
    return num_points;
}


bool SparseVoxelGrid::getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors ) const
{
    ZoneScopedN("SparseVoxelGrid::getSensorCell");
    cellPtrs.clear();
    if ( search_lvl < 0 || search_lvl >= LevelIndexType(m_maps.size()) ) { LOG(1) << "ooL? " << pt_s.transpose() << " lvl: " << search_lvl; return false; } // check lvl bounds
    if ( ! m_map_params.isInBounds ( pt_s , search_lvl ) ) {  LOG(1) << "oob: " << pt_s.transpose() << " lvl: " << search_lvl; return false; }

    cellPtrs.reserve(27);

    constexpr int min_cell = -(m_map_params.NUM_CELLS>>1);
    constexpr int max_cell = (m_map_params.NUM_CELLS>>1)-1;

    const Eigen::Vector3i search_idx_vec = m_map_params.toCellIndexVector( pt_s, search_lvl );
    const auto & map_level = m_maps_surfel_info[search_lvl];

    const int neg_neighbors_z = search_idx_vec.z()-neighbors  < min_cell ? 0 : -neighbors;
    const int pos_neighbors_z = search_idx_vec.z()+neighbors >= max_cell ? 0 : +neighbors;
    const int neg_neighbors_y = search_idx_vec.y()-neighbors  < min_cell ? 0 : -neighbors;
    const int pos_neighbors_y = search_idx_vec.y()+neighbors >= max_cell ? 0 : +neighbors;
    const int neg_neighbors_x = search_idx_vec.x()-neighbors  < min_cell ? 0 : -neighbors;
    const int pos_neighbors_x = search_idx_vec.x()+neighbors >= max_cell ? 0 : +neighbors;

    for ( IndexType iz = neg_neighbors_z; iz <= pos_neighbors_z; ++iz )
    {
        const int vz = search_idx_vec.z()+iz;
        for ( IndexType iy = neg_neighbors_y; iy <= pos_neighbors_y; ++iy )
        {
            const int vy = search_idx_vec.y()+iy;
            for ( IndexType ix = neg_neighbors_x; ix <= pos_neighbors_x; ++ix )
            {
                const int vx = search_idx_vec.x()+ix;

                const CellIndexType v ({vx,vy,vz});
                const auto it = map_level.find( v );
                if ( it != map_level.end() )
                {
                    cellPtrs.emplace_back(it->second);
                }
            }
        }
    }
    return !cellPtrs.empty();
}

int SparseVoxelGrid::getSensorCells( const SurfelInfoVector & data, const IndexType & neighbors ) const
{
    ZoneScopedN("SparseVoxelGrid::getSensorCells");
    int numNeighbors = 0;
    std::vector<SurfelInfoConstPtrVector> cnPtrs (data.size());

    for ( size_t data_idx = 0; data_idx < data.size(); ++data_idx )
    {
        const SparseVoxelGrid::LevelIndexType & search_lvl = data[data_idx].m_level;
        const Eigen::Vector3f pt_s = data[data_idx].m_surfel->mean_ + data[data_idx].m_center_s;
        SurfelInfoConstPtrVector & cellPtrs = cnPtrs[data_idx];

        getSensorCell(pt_s, search_lvl, cellPtrs, neighbors );
        numNeighbors+= cellPtrs.size();
    }
    return numNeighbors;
}
