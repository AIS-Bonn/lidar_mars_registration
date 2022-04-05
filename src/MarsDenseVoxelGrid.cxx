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
#include "MarsDenseVoxelGrid.h"
#include <absl/container/flat_hash_set.h>
#include <Tracy.hpp>

//#define OPT_THREAD_NUM 3

DenseVoxelGrid::DenseVoxelGrid( const MapParameters & params )
    : m_map_params ( params ), m_maps ( params.m_num_levels ), m_maps_surfel_info ( params.m_num_levels )
{

}

void DenseVoxelGrid::allocate()
{
    ZoneScopedN("DenseVoxelGrid::getDenseSpace");
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        if ( map_level.empty() )
            map_level.resize(std::pow(m_map_params.getNumCells(lvlIdx),3));
        if ( m_maps_surfel_info[lvlIdx].empty() )
            m_maps_surfel_info[lvlIdx].reserve(std::pow(m_map_params.getNumCells(lvlIdx),3));
    }
}

SurfelInfoVector DenseVoxelGrid::iterate() const
{
    const bool should_omit_center = true;
    SurfelInfoVector surfels;
    surfels.clear();
    surfels.reserve(1e5);
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        const bool omit_center = should_omit_center && lvlIdx < m_map_params.m_num_levels-1;
        for ( const CellType & cell : m_maps[lvlIdx] )
        {
            if ( cell.isEmpty() || ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center);
        }
    }
    return surfels;
}


void DenseVoxelGrid::clear()
{
    for ( LevelIndexType lvl = 0; lvl < LevelIndexType(m_maps.size()); ++lvl )
    {
        m_maps[lvl].clear();
        m_maps[lvl].resize(std::pow(m_map_params.getNumCells(lvl),3));
        m_maps_surfel_info[lvl].clear();
        m_maps_surfel_info[lvl].resize(std::pow(m_map_params.getNumCells(lvl),3));
    }
}

int DenseVoxelGrid::update( const SurfelInfoVector & surfels )
{
    int num_points = 0;
    for ( size_t lvl = 0; lvl < m_maps.size(); ++lvl )
    {
        auto & map = m_maps_surfel_info[lvl];
        map.clear();
        map.resize( m_maps[lvl].size(), nullptr);
    }
    for ( const SurfelInfo & surfel : surfels )
    {
        if ( ! surfel.m_surfel ) continue;
        m_maps_surfel_info[surfel.m_level][m_map_params.toCellIndexFromIndexVector(Eigen::Map<const Eigen::Vector3i>(surfel.m_index.data()))] = &surfel;
        num_points += surfel.m_surfel->getNumPoints();
    }
    return num_points;
}

int DenseVoxelGrid::addCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_s1s2, const Eigen::Vector3d & origin )
{
    ZoneScopedN("DenseVoxelGrid::addCloudDense");

    absl::FixedArray<absl::flat_hash_set<IndexType>,MapParameters::maxAllowedMapLevels> changedCells ( m_map_params.m_num_levels );
    const Eigen::Vector3f originf = origin.cast<float>();
    const Eigen::Matrix3Xf pts_s = (pose_s1s2.so3().cast<float>().matrix() * cloud->m_points).colwise() + pose_s1s2.translation().cast<float>();
    const int scan_id = cloud->single_scan() ? cloud->m_scan_id(0) : -1; // TODO
    const int num_pts = cloud->size();
    int num_points = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:num_points)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        absl::flat_hash_set<IndexType> & changed_level_cells = changedCells[lvlIdx];
        changed_level_cells.reserve(0.1*num_pts);
        for ( int pointIdx = 0; pointIdx < num_pts; ++pointIdx)
        {
            const Eigen::Vector3f pt_s = pts_s.col(pointIdx);
            const bool isInside = m_map_params.isInBounds( pt_s, lvlIdx );
            if ( ! isInside ) continue;
            // add data.
            const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector ( pt_s, lvlIdx );
            const IndexType cellIdx = m_map_params.toCellIndexFromIndexVector( cellIdxVec );

            if ( cellIdx > int(map_level.size()) )
                LOG(FATAL)  << " cellIdx: " << cellIdx << " lvl: " << lvlIdx << " s: " << map_level.size();

            if ( map_level[cellIdx].isEmpty() )
            {
                // add this one.
                const Eigen::Vector3f center_s = m_map_params.fromCellIndexToCenter <float>(cellIdx,lvlIdx);
                map_level[cellIdx].create( center_s, (pt_s-originf).normalized(), m_map_params.isInBounds<float>( center_s, lvlIdx+1 ), lvlIdx, CellIndexType({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]}));
            }

            if ( !changed_level_cells.contains(cellIdx) )
            {
                changed_level_cells.insert(cellIdx);
                map_level[cellIdx].addNewScan( scan_id, (pt_s-originf).normalized() );
            }
            map_level[cellIdx].addPoint( pt_s );
            ++num_points;
        }
    }

    int updated = 0;
    {
        ZoneScopedN("DenseVoxelGrid::addCloudDense::UpdateSurfel");
        for ( size_t lvlIdx = 0; lvlIdx < changedCells.size(); ++lvlIdx )
        {
            const absl::flat_hash_set<IndexType> & changed_level_cells = changedCells[lvlIdx];
            auto & map_level = m_maps[lvlIdx];
            for ( const IndexType & cellIdx : changed_level_cells )
            {
                if ( cellIdx < int(map_level.size()))
                {
                    CellType & cell = map_level[cellIdx];
                    cell.updateSurfel( );
                    ++updated;
                }
            }
        }
    }
    LOG(1) << "num_points: " << num_points << " updatedCells: " << updated;
    return num_points;
}

void DenseVoxelGrid::updateCells()
{
    ZoneScopedN("DenseVoxelGrid::updateChangedCells");
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        for ( auto & cell : map_level )
        {
            if ( cell.isUpdated() ) continue;
            cell.updateSurfel(true);
        }
    }
}

void DenseVoxelGrid::getOriginCellShift( const Eigen::Vector3f & rel_translation, const LevelIndexType & lvl, Eigen::Vector3f & map_translation, Eigen::Vector3i & num_cells ) const
{
    const float min_cell_size = m_map_params.getCellSizeOnLevel( lvl );
    num_cells = (rel_translation.array() / min_cell_size).cast<int>();
    map_translation = num_cells.cast<float>() * min_cell_size;
}

void DenseVoxelGrid::removeCloudByScanId( const IndexType & scan_id, int & cellsWithRemovedPoints, int & pointsRemoved, int & cellsRemoved)
{
    cellsRemoved = 0;
    pointsRemoved = 0;
    cellsWithRemovedPoints = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)  reduction(+:cellsWithRemovedPoints) reduction(+:pointsRemoved) reduction(+:cellsRemoved)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        int num_erased = 0;
        int cellsWithRemoved = 0;
        MapLevelType & map_level = m_maps[lvlIdx];
        for ( CellType & cell : map_level )
        {
            if ( cell.isEmpty() ) continue;
            const int prevPoints = cell.getNumPoints();
            const bool removedSome = cell.removeByScanId ( scan_id );
            cellsWithRemoved+=removedSome;
            if ( cell.isEmpty() )
                ++num_erased;
            else
                cell.updateSurfel( true );
            const int afterPoints = cell.getNumPoints();
            pointsRemoved += prevPoints - afterPoints;
        }
        cellsWithRemovedPoints+=cellsWithRemoved;
        cellsRemoved += num_erased;
    }
}

void DenseVoxelGrid::getCellsAdaptive ( SurfelInfoVector & surfels, int & normal_level_cnter, int & coarser_level_cnter ) const
{
    normal_level_cnter = 0;
    coarser_level_cnter = 0;
    for ( LevelIndexType lvlIdx = m_map_params.m_num_levels-1; lvlIdx >=0; --lvlIdx )
    {
        const auto & mapLevel = m_maps[lvlIdx];
        absl::flat_hash_set<CellIndexType> usedCoarserLevel;
        usedCoarserLevel.reserve(mapLevel.size());
        SurfelInfoVector coarser_surfels;
        coarser_surfels.reserve(mapLevel.size());
        for ( const CellType & cell : mapLevel )
        {
            if ( cell.isEmpty() || ! cell.m_surfel ) continue;
            if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10. ) continue;
            if ( cell.m_surfel_info.m_in_center ) continue; // bottom up, coarser ones are chosen in the center already:

            bool checkCoarser = true;
            if ( lvlIdx > 0 && (( cell.m_surfel->eigen_values_(0) < MapParameters::first_ev_planarity_threshold )
                                || ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold )
                                || ( cell.m_surfel->eigen_values_(0) < MapParameters::plane_scale_factor * cell.m_surfel->eigen_values_(1))) )
            {
                // check neighbors:
                const Eigen::Vector3i search_idx_vec_coarse = m_map_params.toCellIndexVector( cell.m_surfel_info.m_center_s, lvlIdx-1 );
                const CellIndexType coarse_index ({search_idx_vec_coarse[0],search_idx_vec_coarse[1], search_idx_vec_coarse[2]});

                if ( usedCoarserLevel.count(coarse_index) > 0 ) // corresponding coarser was already added and this one was checked!
                    continue;

                const IndexType idx = m_map_params.toCellIndexFromIndexVector(Eigen::Map<const Eigen::Vector3i>(cell.m_surfel_info.m_index.data()));

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
                const Eigen::Matrix3Xi v2 = neighborAdd.colwise() + search_idx_vec_coarse*2;
                int numChecked = 0;
                Eigen::Vector3f fine_normal = cell.m_surfel->normal_;
                for ( int nIdx = 0; nIdx < v2.cols(); ++nIdx)
                {
                    const IndexType nI = m_map_params.toCellIndexFromIndexVector(v2.col(nIdx));
                    if ( nI != idx ) // only neighbors, not this one.
                    {
                        const CellType & other_cell = mapLevel[nI];
                        if ( other_cell.isEmpty() || ! other_cell.m_surfel ) continue;
                        if ( ! other_cell.m_surfel->valid_ || other_cell.m_surfel->getNumPoints() < 10.0 ) continue;
                        checkCoarser &=   (other_cell.m_surfel->eigen_values_(0) < MapParameters::first_ev_planarity_threshold)
                                          || (other_cell.m_surfel->eigen_values_(0) < MapParameters::plane_scale_factor * other_cell.m_surfel->eigen_values_(1));
                        fine_normal+= other_cell.m_surfel->normal_;
                        ++numChecked;
                    }
                    if ( ! checkCoarser ) break;
                }
                if ( numChecked == 0 )
                    checkCoarser = false;

                if ( checkCoarser )
                {
                    const CellType & coarser_cell = m_maps[lvlIdx-1][m_map_params.toCellIndexFromIndexVector(Eigen::Map<const Eigen::Vector3i>(coarse_index.data()))];
                    if ( coarser_cell.m_surfel && coarser_cell.m_surfel->valid_ && coarser_cell.m_surfel->getNumPoints() >= 10. &&
                         ( coarser_cell.m_surfel->eigen_values_(0) < MapParameters::first_ev_planarity_threshold
                           || (coarser_cell.m_surfel->eigen_values_(0) < MapParameters::plane_scale_factor * coarser_cell.m_surfel->eigen_values_(1))
                           ) && std::abs(fine_normal.normalized().dot(coarser_cell.m_surfel->normal_)) > 0.8
                         && (coarser_cell.m_surfel->eigen_values_(1) > MapParameters::second_ev_degenerate_threshold) )
                    {
                        // use coarser one.
                        coarser_surfels.emplace_back(coarser_cell.m_surfel, coarser_cell.m_surfel_info.m_center_s, lvlIdx-1,  coarse_index, coarser_cell.m_cov_per_scan.empty() ? -1 : coarser_cell.m_cov_per_scan.front().first, coarser_cell.m_surfel_info.m_in_center, coarser_cell.m_class );
                        ++coarser_level_cnter;
                        usedCoarserLevel.insert(coarse_index);
                    }
                    else
                        checkCoarser = false;
                }
            }
            else
                checkCoarser = false;
            if ( !checkCoarser )
            {
                if ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
                ++normal_level_cnter;
            }
        }
        surfels.insert(surfels.end(),std::make_move_iterator(coarser_surfels.begin()),std::make_move_iterator(coarser_surfels.end()));
    }
}

void DenseVoxelGrid::getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center ) const
{
    for ( const CellType & cell : m_maps[lvlIdx] )
    {
        if ( cell.isEmpty() || ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, cell.m_surfel_info.m_index, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

void DenseVoxelGrid::getCellsFusedOnLevel ( const IndexType & lvlIdx, SurfelInfoConstPtrVector & surfels, const bool & omit_center ) const
{
    for ( const CellType & cell : m_maps[lvlIdx] )
    {
        if ( cell.isEmpty() || ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        surfels.emplace_back(cell.getSurfelInfo());
    }
}

void DenseVoxelGrid::getCellsScanWiseOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center, const LevelIndexType & level_modifier  ) const
{
    for ( const CellType & cell : m_maps[lvlIdx] )
    {
        if ( cell.isEmpty() || ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : cell.m_cov_per_scan)
            surfels.emplace_back(surfelInfo.second, cell.m_surfel_info.m_center_s, lvlIdx+level_modifier, cell.m_surfel_info.m_index, surfelInfo.first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

template <typename SurfelInfoT>
int DenseVoxelGrid::addCells ( const std::vector<SurfelInfoT> & surfels, const Sophus::SE3d & orig_pose_s1s2, const LevelIndexType & level_modifier )
{
    ZoneScopedN("DenseVoxelGrid::addCells");
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

inline int DenseVoxelGrid::addCellSurfel ( const SurfelInfo & surfel_info, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier, Sophus::SE3d & n_pose_s1s2 )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;
    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & center_s = surfel_info.m_center_s;
    const Eigen::Vector3f meanPt = (pose_s1s2 * (surfel->mean_ + center_s).template cast<double>()).template cast<float>();
    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    const Eigen::Vector3i cellIdxVec = m_map_params.toCellIndexVector ( meanPt, lvlIdx );
    const IndexType cellIdx = m_map_params.toCellIndexFromIndexVector( cellIdxVec );
    CellType & cell = map_level[cellIdx];
    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > std::sqrt(std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            return 0;
        }
        cell.create( new_cell_center, (pose_s1s2.translation().template cast<float>()-meanPt).normalized(), m_map_params.isInBounds<float>( meanPt, lvlIdx+1), lvlIdx, CellIndexType({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]}) );
    }

    // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
    // mean_new = meanPt - center_new - surfel mean

    n_pose_s1s2.translation() =  (meanPt - surfel->mean_ - cell.m_surfel_info.m_center_s).cast<double>();

    cell.addSurfel(surfel, n_pose_s1s2, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}

inline int DenseVoxelGrid::addCellSurfelOnGrid( const SurfelInfo & surfel_info, const Eigen::Vector3f & t_s1s2, const LevelIndexType & level_modifier, const Sophus::SE3d & id )
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
    const IndexType cellIdx = m_map_params.toCellIndexFromIndexVector( cellIdxVec );
    CellType & cell = map_level[cellIdx];
    if ( cell.isEmpty() )
    {
        const Eigen::Vector3f new_cell_center = m_map_params.fromCellIdxToCenter<float>(cellIdxVec,lvlIdx);
        if ( (new_cell_center-meanPt).squaredNorm() > (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3) )
        {
            LOG(FATAL) << "ncc: " << new_cell_center.transpose() <<" mp: " << meanPt.transpose() << " ct: " << old_center.transpose() << " cs: " << m_map_params.getCellSizeOnLevel(lvlIdx) << " d2: " << (new_cell_center-meanPt).squaredNorm() << " e: " << (std::pow(m_map_params.getCellSizeOnLevel(lvlIdx),2)/4*3);
        }
        // add this one.
        cell.create( new_cell_center, (t_s1s2-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, CellIndexType({cellIdxVec[0],cellIdxVec[1],cellIdxVec[2]}) );
    }

    // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
    // mean_new = meanPt - center_new - surfel mean

    // -> this does not matter, since the shift is on the grid!
    cell.addSurfel(surfel, id, surfel_info.m_scan_id, surfel_info.m_class);
    return surfel->getNumPoints();
}

template <typename SurfelInfoT>
int DenseVoxelGrid::addCellsOnGrid ( const std::vector<SurfelInfoT> & surfels, const Eigen::Vector3d & orig_t_s1s2, const LevelIndexType & level_modifier )
{
    // Assumes surfels are on the same grid with a shift of n_i per dimension. -> no rotation or fractions within translation
    ZoneScopedN("DenseVoxelGrid::addCells");
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

template int DenseVoxelGrid::addCellsOnGrid<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Eigen::Vector3d & , const LevelIndexType &);
template int DenseVoxelGrid::addCellsOnGrid<SurfelInfo>( const SurfelInfoVector & , const Eigen::Vector3d & , const LevelIndexType & );
template int DenseVoxelGrid::addCells<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Sophus::SE3d & , const LevelIndexType & );
template int DenseVoxelGrid::addCells<SurfelInfo>( const SurfelInfoVector & , const Sophus::SE3d & , const LevelIndexType & );

bool DenseVoxelGrid::getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors ) const
{
    ZoneScopedN("DenseVoxelGrid::getSensorCell");
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
                const auto v = m_map_params.toCellIndexFromIndexVector( Eigen::Vector3i(vx,vy,vz) );
                const auto & c = map_level[v];
                if ( c != nullptr )
                    cellPtrs.emplace_back(c);
            }
        }
    }
    return !cellPtrs.empty();
}

int DenseVoxelGrid::getSensorCells( const SurfelInfoVector & data, const IndexType & neighbors ) const
{
    ZoneScopedN("DenseVoxelGrid::getSensorCells");
    int numNeighbors = 0;
    std::vector<SurfelInfoConstPtrVector> cnPtrs (data.size());
    for ( size_t data_idx = 0; data_idx < data.size(); ++data_idx )
    {
        const DenseVoxelGrid::LevelIndexType & search_lvl = data[data_idx].m_level;
        const Eigen::Vector3f pt_s = data[data_idx].m_surfel->mean_ + data[data_idx].m_center_s;
        SurfelInfoConstPtrVector & cellPtrs = cnPtrs[data_idx];

        getSensorCell(pt_s, search_lvl, cellPtrs, neighbors );
        numNeighbors+= cellPtrs.size();
    }
    return numNeighbors;
}
