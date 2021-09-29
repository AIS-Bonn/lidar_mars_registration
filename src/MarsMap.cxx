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
#include "MarsMap.h"

#include <loguru.hpp>

#include <Tracy.hpp>

int UID::get() {
    static std::atomic<std::uint32_t> uid { 0 };  // <<== initialised
    //    uid = 0;    <<== removed
    return ++uid;
}

template <typename T>
typename MarsMapBase<T>::Ptr MarsMapBase<T>::create ( const MarsMapBase<T> & otherMap, const bool & allocate )
{
    MarsMapBase<T>::Ptr map = std::make_shared<MarsMapBase<T>>(otherMap.m_map_params, allocate );
    map->initParams(otherMap);
    return map;
}

template <typename T>
typename MarsMapBase<T>::Ptr MarsMapBase<T>::create ( const MapParameters & params, const bool & allocate )
{
    return std::make_shared<MarsMapBase<T>>(params, allocate);
}

template <typename T>
MarsMapBase<T>::MarsMapBase ( const MapParameters & params, const bool & allocate )
    : m_map_params ( params ), m_storage ( params )
{
    m_id = UID::get();
//    LOG(INFO) << "creating map["<<m_id<<"]. lvls: " << m_map_params.m_num_levels << " " << m_storage.m_maps.size()  <<" mapSideLength: " << m_map_params.m_size << " cls: " << m_map_params.m_num_cells
//              << " min=[" << m_map_params.m_min_origin_bounds.transpose() << "] to max=[" << m_map_params.m_max_origin_bounds.transpose() << "]";
    if ( allocate )
        m_storage.allocate();
}
template <typename T>
void MarsMapBase<T>::update( const bool & use_adaptive )
{
    m_surfels.clear();
    bool adaptive_update = false;
    Sophus::SE3d pose_w;
    {
        ZoneScopedN("MarsMap::updateAdapt");
        if ( use_adaptive && m_map_params.m_omit_center ) // adaptive does not make sense without omit center.
        {
            MarsMapBase<T>::getCellsAdaptive(m_surfels, pose_w);
            //LOG(INFO) << "Map. updating adaptive: " << m_surfels.size();
            adaptive_update = true;
        }
        else
        {
            MarsMapBase<T>::getCells(m_surfels, pose_w, m_map_params.m_omit_center);
            //LOG(INFO) << "Map. updating normal: " << m_surfels.size();
        }
    }
    m_num_points = m_storage.update(m_surfels);
    m_adaptor_updated = true;

    LOG(INFO) << "Map. updated: "<< m_num_points << " adaptive: " << adaptive_update << " surfels: " << m_surfels.size();
}

template <typename T>
void MarsMapBase<T>::initParams ( const MarsMapBase<T> & otherMap )
{
    m_pose_init_w = otherMap.m_pose_init_w;
    m_pose_w = otherMap.m_pose_w;
    m_num_points = 0;
    m_surfels.clear();
    m_storage.clear();
}

template <typename T>
void MarsMapBase<T>::clear()
{
    ZoneScopedN("MarsMap::clear");
    m_pose_init_w = Sophus::SE3d();
    m_pose_w = Sophus::SE3d();
    m_adaptor_updated = false;
    m_num_points = 0;
    m_storage.clear();
}

template <typename T>
void MarsMapBase<T>::setPose ( const Sophus::SE3d & pose_w )
{
    m_pose_w = pose_w;
}

template <typename T>
void MarsMapBase<T>::transform ( const Sophus::SE3d & pose_w )
{
    m_pose_w = pose_w * m_pose_w;
}
template <typename T>
void MarsMapBase<T>::transform_shift ( const Eigen::Vector3d & t_w )
{
    m_pose_w.translation() += t_w;
}

template <typename T>
void MarsMapBase<T>::setCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_w, const Eigen::Vector3d & scan_origin  )
{
    m_adaptor_updated = false;
    m_pose_init_w = pose_w;
    m_scan_origin = scan_origin;
    addCloud ( cloud ); // at identity
    transform ( pose_w );
}

template <typename T>
void MarsMapBase<T>::addCloud ( MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_w2 )
{
    if ( ! cloud ) return;
    m_adaptor_updated = false;

    ZoneScopedN("MarsMap::addCloud");
    const Sophus::SE3d pose_s1s2 = m_pose_w.inverse() * pose_w2;
    const Eigen::Vector3d origin = m_scan_origin + pose_s1s2.translation();

    m_num_points += m_storage.addCloud ( cloud, pose_s1s2, origin );
}

template <typename T>
void MarsMapBase<T>::removeCloudByScanId ( const IndexType & scan_id )
{
    int cellsWithRemovedPoints = 0, pointsRemoved = 0, cellsRemoved = 0;
    m_storage.removeCloudByScanId( scan_id, cellsWithRemovedPoints, pointsRemoved, cellsRemoved );
    m_num_points -= pointsRemoved;
    if ( cellsRemoved > 0 ) m_adaptor_updated = false;
    LOG(INFO) << "removed cloud with id: " << scan_id << " removed some: " <<cellsWithRemovedPoints;
}

template <typename T>
void MarsMapBase<T>::getCellsAdaptive ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w) const
{
    int normal_level_cnter = 0, coarser_level_cnter = 0;
    surfels.clear();
    surfels.reserve(1e4);
    m_storage.getCellsAdaptive( surfels, normal_level_cnter, coarser_level_cnter );
    pose_w = m_pose_w;
    LOG(INFO) << "used from coarser level: " << coarser_level_cnter << " used from normal level: " << normal_level_cnter;
}

template <typename T>
void MarsMapBase<T>::getCells ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center) const
{
    surfels.clear();
    surfels.reserve(1e4);
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        m_storage.getCellsOnLevel(lvlIdx, surfels, omit_center && lvlIdx < m_map_params.m_num_levels-1);
    }
    pose_w = m_pose_w;
}

template <typename T>
void MarsMapBase<T>::getCellsScanWise ( SurfelInfoVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center, const LevelIndexType & level_modifier ) const
{
    surfels.clear();
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        m_storage.getCellsScanWiseOnLevel(lvlIdx, surfels, omit_center && lvlIdx < m_map_params.m_num_levels-1, level_modifier);
    }
    pose_w = m_pose_w;
}

template <typename T>
void MarsMapBase<T>::getCellsFused ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w, const bool & omit_center ) const
{
    surfels.clear();
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        m_storage.getCellsFusedOnLevel(lvlIdx, surfels, omit_center && lvlIdx < m_map_params.m_num_levels-1);
    }
    pose_w = m_pose_w;
}

template <typename T>
void MarsMapBase<T>::getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center ) const
{
    m_storage.getCellsOnLevel(lvlIdx, surfels, omit_center);
}

template <typename T>
double MarsMapBase<T>::getCellSize ( const MarsMapBase<T>::LevelIndexType & lvl ) const
{
    return m_map_params.getCellSizeOnLevel(lvl);
}

template <typename T>
void MarsMapBase<T>::getCells ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w ) //const
{
    if ( ! m_adaptor_updated )
        update();
    pose_w = m_pose_w;
    surfels.clear();
    surfels.reserve(m_surfels.size());
    //int num_classes = 0;
    for  ( const SurfelInfo & surfel : m_surfels )
    {
        surfels.emplace_back(&surfel);
        //num_classes += surfel.m_class != nullptr;
    }
    //LOG(INFO) << "getCells: num_clases " << num_classes << " of " << m_surfels.size();
}

template <typename T>
void MarsMapBase<T>::getCellsAdaptive ( SurfelInfoConstPtrVector & surfels, Sophus::SE3d & pose_w )
{
    if ( ! m_adaptor_updated )
        update(true);
    pose_w = m_pose_w;
    surfels.clear();
    surfels.reserve(m_surfels.size());
    //int num_classes = 0;
    for  ( const SurfelInfo & surfel : m_surfels )
    {
        surfels.emplace_back(&surfel);
        //num_classes += surfel.m_class != nullptr;
    }
    //LOG(INFO) << "getCells: num_clases " << num_classes;
}

template <typename T>
void MarsMapBase<T>::addCells ( const SurfelInfoVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    const Sophus::SE3d pose_s1s2 = m_pose_w.inverse() * pose_w2;

    LOG(INFO) << "pose_w: " << m_pose_w.params().transpose() << " pose_w2: " << pose_w2.params().transpose() << " pose_s1s2: " << pose_s1s2.params().transpose();
    m_num_points += m_storage.addCells( surfels, pose_s1s2, level_modifier );
    if ( shouldUpdate )
    {
        ZoneScopedN("MarsMap::updateChangedCells");
        updateCells();
    }
}
template <typename T>
void MarsMapBase<T>::addCells ( const SurfelInfoConstPtrVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    const Sophus::SE3d pose_s1s2 = m_pose_w.inverse() * pose_w2;

    LOG(INFO) << "pose_w: " << m_pose_w.params().transpose() << " pose_w2: " << pose_w2.params().transpose() << " pose_s1s2: " << pose_s1s2.params().transpose();
    m_num_points += m_storage.addCells( surfels, pose_s1s2, level_modifier );
    if ( shouldUpdate )
    {
        ZoneScopedN("MarsMap::updateChangedCells");
        updateCells();
    }
}

template <typename T>
void MarsMapBase<T>::addCellsOnGrid ( const SurfelInfoVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    const Eigen::Vector3d t_s1s2 = -m_pose_w.translation() + t_w2;

    LOG(INFO) << "pose_w: " << m_pose_w.params().transpose() << " t_w2: " << t_w2.transpose() << " t_s1s2: " << t_s1s2.transpose();
    m_num_points += m_storage.addCellsOnGrid( surfels, t_s1s2, level_modifier );
    if ( shouldUpdate )
    {
        ZoneScopedN("MarsMap::updateChangedCells");
        updateCells();
    }
}

template <typename T>
void MarsMapBase<T>::addCellsOnGrid ( const SurfelInfoConstPtrVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    const Eigen::Vector3d t_s1s2 = -m_pose_w.translation() + t_w2;

    LOG(INFO) << "pose_w: " << m_pose_w.params().transpose() << " t_w2: " << t_w2.transpose() << " t_s1s2: " << t_s1s2.transpose();
    m_num_points += m_storage.addCellsOnGrid( surfels, t_s1s2, level_modifier );
    if ( shouldUpdate )
    {
        ZoneScopedN("MarsMap::updateChangedCells");
        updateCells();
    }
}

template <typename T>
void MarsMapBase<T>::setCells ( const SurfelInfoVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    addCells( surfels, Sophus::SE3d(), shouldUpdate, level_modifier );
    transform ( pose_w2 );
}

template <typename T>
void MarsMapBase<T>::setCells ( const SurfelInfoConstPtrVector & surfels, const Sophus::SE3d & pose_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    addCells( surfels, Sophus::SE3d(), shouldUpdate, level_modifier );
    transform ( pose_w2 );
}

template <typename T>
void MarsMapBase<T>::setCellsOnGrid ( const SurfelInfoVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    addCellsOnGrid( surfels, Eigen::Vector3d::Zero(), shouldUpdate, level_modifier );
    transform_shift ( t_w2 );
}
template <typename T>
void MarsMapBase<T>::setCellsOnGrid ( const SurfelInfoConstPtrVector & surfels, const Eigen::Vector3d & t_w2, const bool & shouldUpdate, const LevelIndexType & level_modifier )
{
    m_adaptor_updated = false;
    addCellsOnGrid( surfels, Eigen::Vector3d::Zero(), shouldUpdate, level_modifier );
    transform_shift ( t_w2 );
}

template <typename T>
void MarsMapBase<T>::updateCells ( )
{
    m_storage.updateCells();
    LOG(INFO) << "updated cells.";
}

template <typename T>
bool MarsMapBase<T>::getSensorCell ( const Eigen::Vector3f & pt_s, const MarsMapBase<T>::LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors) const
{
    cellPtrs.clear();
    //LOG(INFO) << "search_level: " << search_lvl << " pt: " << pt_w.transpose();
    if ( search_lvl < 0 || search_lvl >= LevelIndexType(m_storage.m_maps.size()) ) { /*LOG(INFO) << "ooL? " << pt_s.transpose() << " lvl: " << search_lvl; */ return false; } // check lvl bounds
    //LOG(INFO) << "pt_s: " << pt_s.transpose();
    if ( ! m_map_params.isInBounds<float> ( pt_s , search_lvl ) ) { /* LOG(INFO) << "oob: " << pt_s.transpose() << " lvl: " << search_lvl; */ return false; }
    //LOG(INFO) << "inside.";
    cellPtrs.reserve(27);

    return m_storage.getSensorCell(pt_s, search_lvl, cellPtrs, neighbors);
}

template <typename T>
const Sophus::SE3d & MarsMapBase<T>::getMapPose() const
{
    return m_pose_w;
}
template <typename T>
const Sophus::SE3d & MarsMapBase<T>::getMapInitPose() const
{
    return m_pose_init_w;
}

template <typename T>
size_t MarsMapBase<T>::getNumPoints ( ) const
{
    return m_num_points;
}

template <typename T>
SurfelInfoVector MarsMapBase<T>::iterate()
{
    if ( ! m_adaptor_updated )
        update ( m_map_params.m_use_adaptive );
    return m_storage.iterate();
}

template <typename T>
int MarsMapBase<T>::getSensorCells( const SurfelInfoVector & data )
{
    if ( ! m_adaptor_updated )
        update ( m_map_params.m_use_adaptive );
    return m_storage.getSensorCells(data, 1);
}

template <typename T>
Eigen::Vector3i MarsMapBase<T>::getOriginCellShift( const Sophus::SE3d & map_pose, Sophus::SE3d & new_map_pose, Sophus::SE3d & map_integration_pose, const int & lvl )
{
    new_map_pose = map_pose;
    map_integration_pose = Sophus::SE3d();
    const Sophus::SE3d rel_map_pose = getMapPose().inverse() * map_pose;
    map_integration_pose.so3() = rel_map_pose.so3();
    const Eigen::Vector3f rel_translation = rel_map_pose.translation().template cast<float>();

    Eigen::Vector3i num_cells = Eigen::Vector3i::Zero();
    Eigen::Vector3f map_translation = Eigen::Vector3f::Zero();
    m_storage.getOriginCellShift( rel_translation, lvl, map_translation, num_cells );

    // move map only on non full numbers of grid cells -> aligns grids
    map_integration_pose.translation() = (rel_translation - map_translation).template cast<double>();
    new_map_pose.so3() = getMapPose().so3();
    new_map_pose.translation() = getMapPose() * map_translation.template cast<double>();
    return num_cells.array().abs();
}


//template class MarsMapBase<DenseVoxelGrid>;
template class MarsMapBase<SparseVoxelGrid>;
template class MarsMapBase<BlockSparseVoxelGrid>;
template class MarsMapBase<Lattice>;
