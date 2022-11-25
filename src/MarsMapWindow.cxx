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
#include "MarsMapWindow.h"
#include "MarsMap.h"

#ifdef USE_EASY_PBR
#include "VisUtils.h"
#endif

#include <absl/container/flat_hash_set.h>

#define LOGURU_REPLACE_GLOG 1
#include <loguru.hpp>
#include <Tracy.hpp>

//#define OPT_THREAD_NUM 3
#define USE_TBB
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

MarsMapWindow::MarsMapWindow ( const size_t & cur_queue_size, const bool & map_per_cloud, const bool & use_adaptive)
    : queue_size ( cur_queue_size ),
      store_map_per_cloud ( map_per_cloud ),
      m_use_adaptive ( use_adaptive ),
      oldMaps( map_per_cloud ? cur_queue_size : ( MapParameters::getMapParameters().m_num_levels ), nullptr )
{
    const MapParameters standard_params = MapParameters::getMapParameters();
    const float map_side_length = standard_params.getMapSideLength(0); // coarsest side length
    const float coarsest_lvl_cell_size = standard_params.getCellSizeOnLevel(0);
    const int num_levels = standard_params.m_num_levels;
    const int num_cells = standard_params.m_num_cells;
    // scene map
    if ( map_per_cloud )
    {
        const MapParameters params = MapParameters::getMapParameters(map_side_length, num_cells, num_levels, use_adaptive, true);
        std::generate( oldMaps.begin(), oldMaps.end(), [params](){ return MarsMapType::create(params); } ); // this one omits the centers
    }
    else
    {
            // local map
            for ( size_t i = 0; i < oldMaps.size(); ++i )
            {
                constexpr bool increase_bound = false; // should just add the boundary towards the nearest larger size.
                const int lvl = (i%num_levels);
                const float cur_lvl_cell_size = standard_params.getCellSizeOnLevel(lvl);
                const int cur_lvl_cells_in_coarsest = std::round<int>(coarsest_lvl_cell_size / cur_lvl_cell_size);
                const int num_local_multi_resolution_cells = increase_bound && lvl > 0 ? m_min_num_cells_for_moving_map * cur_lvl_cells_in_coarsest + num_cells : num_cells;
                const int cur_lvl_num_cells = num_local_multi_resolution_cells;
                const MapParameters oldParams = MapParameters::getMapParameters( standard_params.getMapSideLength(lvl), cur_lvl_num_cells, 1 );
                oldMaps[oldMaps.size()-1-i] = MarsMapType::create( oldParams );
            }
            asOneMap = MarsMapType::create( MapParameters::getMapParameters( map_side_length, num_cells, num_levels, use_adaptive, false ) ); // sets omit_center to false! which is important to have all surfels available
            asOneTmpMap = MarsMapType::create( asOneMap->m_map_params );
    }
    LOG(1) << "MarsMapWindow: numMaps: "<< oldMaps.size();
}

MarsMapWindow::Ptr MarsMapWindow::create( const size_t & cur_queue_size, const bool & store_map_per_cloud, const bool & use_adaptive )
{
    return std::make_shared<MarsMapWindow>( cur_queue_size, store_map_per_cloud, use_adaptive );
}

std::vector<MarsMapType*> MarsMapWindow::getMapPtrVec()
{
    if ( ! store_map_per_cloud ) return std::vector<MarsMapType *>();
    std::vector<MarsMapType*> m;
    for ( auto it = maps.begin(); it != maps.end(); ++it)
        if ( *it )
            m.emplace_back(it->get());
    return m;
}

std::vector<MarsMapType*> MarsMapWindow::getTransformedMapPtrVec( const Sophus::SE3d & pose )
{
    std::vector<MarsMapType*> m ( maps.size(), nullptr );
    auto adapIt = m.begin();
    auto mapsIt = maps.begin();
    for ( ; mapsIt != maps.end() && adapIt != m.end(); ++mapsIt, ++adapIt)
    {
        if ( ! *mapsIt ) LOG(FATAL) << "why is there an empty cloud ptr?";
        (*adapIt) = mapsIt->get();
        (*adapIt)->setPose( pose );
        //LOG(1) << "Map has: " << (*adapIt)->m_id << " Map has: " << (*mapsIt)->m_id;
    }
    return m;
}
std::vector<MarsMapPointCloud*> MarsMapWindow::getCloudPtrVec()
{
    std::vector<MarsMapPointCloud*> c;
    for ( auto it = clouds.begin(); it != clouds.end(); ++it)
        if ( *it )
            c.emplace_back(it->get());
    return c;
}
std::vector<MarsMapPointCloud::Ptr> MarsMapWindow::getCloudSharedPtrVec()
{
    std::vector<MarsMapPointCloud::Ptr> c;
    for ( auto it = clouds.begin(); it != clouds.end(); ++it)
        if ( *it )
            c.emplace_back(*it);
    return c;
}
std::vector<Sophus::SE3d> MarsMapWindow::getCloudPosesVec() const
{
    std::vector<Sophus::SE3d> p;
    for ( auto it = clouds_pose.begin(); it != clouds_pose.end(); ++it)
        p.emplace_back(*it);
    return p;
}

Eigen::VectorXd MarsMapWindow::getTimeFactors() const
{
    if ( clouds.empty() ) return Eigen::VectorXd();
    std::vector<double> timeFactors;
    for ( auto it = clouds.begin(); it != clouds.end(); ++it )
    {
        if ( *it)
        {
            timeFactors.emplace_back((*it)->m_stamp);
        }
    }
    Eigen::VectorXd timeScalingFactor = Eigen::VectorXd::Zero(timeFactors.size(),1);
    if ( timeFactors.size() == 1 )
    {
        timeFactors.front() = 1.;
        timeScalingFactor(0) = 1;
    }
    else
    {
        const double last = last_ts;
        const double latest = timeFactors.back();
        const double scaleVal = std::max( latest - last, 1e-5 );
        for ( size_t i = 0; i < timeFactors.size(); ++i )
        {
            timeScalingFactor(i) = (timeFactors[i]-last) / scaleVal;
            if ( timeScalingFactor(i) < 0 || timeScalingFactor(i) > 1 )
                LOG(FATAL) << "Factors are out of reach! tF[i]= " << timeFactors[i]<< " " << timeScalingFactor(i) << " s: " << scaleVal << " f: " << latest << " l: "<< last;
        }
    }
    if ( !timeScalingFactor.allFinite() ) LOG(FATAL) << "Factors are not finite: " << timeScalingFactor.transpose();
    return timeScalingFactor;
}
double MarsMapWindow::getLastTime() const
{
    double stamp = 0;
    for ( auto it = clouds.begin(); it != clouds.end(); ++it )
    {
        if ( *it && (*it)->m_stamp > stamp )
            stamp = (*it)->m_stamp;
    }
    return stamp;
}

Eigen::VectorXt MarsMapWindow::getTimesSince ( const uint64_t & last_ts ) const
{
    std::vector<uint64_t> times = getTimes();
    Eigen::VectorXt timesSinceLast = Eigen::VectorXt::Zero(times.size(),1);
    for ( size_t i = 0; i < times.size(); ++i )
        timesSinceLast[i] = (times[i] - last_ts);
    //LOG(1) << "getTimesSinceLast: last: " << last_ts << " tsl: " << timesSinceLast.transpose();
    return timesSinceLast;
}

std::vector<uint64_t> MarsMapWindow::getTimes( ) const
{
    if ( clouds.empty() ) return std::vector<uint64_t>();
    std::vector<uint64_t> stamps;
    for ( auto it = clouds.begin(); it != clouds.end(); ++it )
    {
        if ( *it)
        {
            stamps.emplace_back((*it)->m_stamp);
        }
    }
    if ( stamps.empty() ) LOG(FATAL) << "No timestamps in the map window? This should not happen!";
    return stamps;
}

bool MarsMapWindow::hasOutOfWindow() const
{
    return m_out_of_window_cloud != nullptr;
}

MarsMapPointCloud::Ptr MarsMapWindow::getOutOfWindowCloud() const
{
    return m_out_of_window_cloud;
}
Sophus::SE3d MarsMapWindow::getOutOfWindowCloudPose() const
{
    return m_out_of_window_cloud_pose;
}

void MarsMapWindow::transform ( const Sophus::SE3d & tf)
{
    for ( auto & map : maps )
        map->transform ( tf );
}

MarsMapPointCloud::Ptr MarsMapWindow::getCurrentOriginCloud() const
{
    if ( clouds.empty() ) return nullptr;
    return clouds.front();
}

bool MarsMapWindow::hasFullWindow() const
{
    return clouds.size() >= queue_size;
}

bool MarsMapWindow::isEmpty() const
{
    return clouds.empty();
}

bool MarsMapWindow::moveWindowToNext()
{
    if ( !clouds.front() ) LOG(FATAL) << "why? is this not true?";

    // this only removes clouds if queue really full, but ensures the maps / oldMaps are popped/emplaced
    const bool shouldRemove = hasFullWindow();

    if ( shouldRemove )
    {
        m_out_of_window_cloud = clouds.front();
        m_out_of_window_cloud_pose = clouds_pose.front();
    }

    const int iters = store_map_per_cloud ? 1 : asOneMap->m_map_params.m_num_levels;
    for ( int i  = 0; i < iters; ++i )
    {
        //LOG(1) << "maps.size():" << maps.size() << " oldMaps.size(): " << oldMaps.size() << " param: " << maps.front()->m_map_params.m_num_cells << " " << maps.front()->m_map_params.m_size;
        oldMaps.emplace_front(maps.front());
        maps.pop_front();
    }

    if ( shouldRemove ) // this only removes clouds if queue really full, but ensures the maps / oldMaps are popped/emplaced
    {
        last_ts = clouds.front()->m_stamp;
        clouds_pose.pop_front();
        clouds.pop_front();
    }
    return shouldRemove;
}

Sophus::SE3d MarsMapWindow::getClosestKeyFramePose( const Sophus::SE3d & scene_pose ) const
{
    Sophus::SE3d closest_pose;
    if ( clouds.empty() )
        return closest_pose;
    double min_map_distance = std::numeric_limits<double>::max();
    for ( const Sophus::SE3d & map_pose : clouds_pose )
    {
        const double cur_map_distance = (map_pose.translation() - scene_pose.translation()).squaredNorm();
        if ( cur_map_distance < min_map_distance )
        {
            min_map_distance = cur_map_distance;
            closest_pose = map_pose;
        }
    }
    return closest_pose;
}

void MarsMapWindow::addCloud ( MarsMapPointCloud::Ptr newCloud, const Sophus::SE3d & map_pose )
{
    if ( !newCloud ) LOG(FATAL) << "empty cloud.";
    //LOG(1) << "cloud size: " << clouds.size() << " " << maps.size() << " qs: " << queue_size;
    if ( store_map_per_cloud && maps.size() != clouds.size() ) LOG(FATAL) << " sizes are inconsistent.";

    local_map_moved = false;
    const bool with_forgetting = true;
    m_out_of_window_cloud = nullptr;
    m_out_of_window_cloud_pose = Sophus::SE3d();
    if ( hasFullWindow() || ( !maps.empty() && ! store_map_per_cloud ) )
    {
        //LOG(1) << "Moving Local Map.";
        const bool thereIsSomethingToForget = moveWindowToNext();
        if ( with_forgetting && thereIsSomethingToForget )
        {
            if ( ! store_map_per_cloud && m_out_of_window_cloud )
            {
                //LOG(1) << "finding the right ones to forget.";
                ZoneScopedN("MarsMapWindow::remove_old_cloud");
                absl::flat_hash_set<MarsMapType::IndexType> scan_ids;
                if ( m_out_of_window_cloud->single_scan() )
                {
                    asOneMap->removeCloudByScanId ( m_out_of_window_cloud->m_scan_id(0) );
                }
                else
                {
                    const int num_elems = m_out_of_window_cloud->size();
                    //                for ( const MarsMapPoint & pt : m_out_of_window_cloud->points )
                    for ( int idx = 0; idx < num_elems; ++idx )
                        if ( scan_ids.find(m_out_of_window_cloud->m_scan_id(idx)) == scan_ids.end() )
                            scan_ids.insert(m_out_of_window_cloud->m_scan_id(idx));
                    for ( const MarsMapType::IndexType & scan_id : scan_ids )
                        asOneMap->removeCloudByScanId ( scan_id );
                }
            }
        }
        //LOG(1) << "Moved and removed from Local Map.";
    }
    clouds.emplace_back(newCloud);
    clouds_pose.emplace_back(map_pose);

    if ( store_map_per_cloud )
    {
        // scene map
        ZoneScopedN("MarsMapWindow::AddSingleCloud");
        maps.emplace_back(oldMaps.back());
        oldMaps.pop_back();
        {
            ZoneScopedN("MarsMapWindow::AddSingleCloud::clearPrev");
            maps.back()->clear();
            //maps.back() = MarsMapType::create(maps.back()->m_map_params);
        }
        {
            ZoneScopedN("MarsMapWindow::AddSingleCloud::setCloud");
            maps.back()->setCloud( newCloud, map_pose );
        }

#ifdef USE_EASY_PBR
        constexpr bool check_map_cells = false;
        if ( check_map_cells )
        {
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector cells;
            maps.back()->getCells( cells, pose );
            int num_invalid = 0, num_less_pts = 0, num_valid = 0;
            for ( const SurfelInfoConstPtr & cell : cells )
            {
                if ( ! cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                num_invalid += !cell->m_surfel->valid_;
                num_less_pts += (cell->m_surfel->getNumPoints() < 10.0);
    //            if ( cell->m_surfel->getNumPoints() < 10.0 )
    //                LOG(1) <<"MarsMapWindow:: Cell: lvl: " << cell->m_level << " idx: " << Eigen::Map<const Eigen::Vector3i>(cell->m_index.data()).transpose() << " center: " << cell->m_center_s.transpose() << " pts: "<<  cell->m_surfel->getNumPoints();
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                ++num_valid;
            }
            LOG(1) << "MarsMapWindow:: Map-Cells: invalid: " << num_invalid << " less: " << num_less_pts << " valid: " << num_valid;
            maps.back()->update(true);
            maps.back()->getCells( cells, pose );
            num_invalid = 0, num_less_pts = 0, num_valid = 0;
            for ( const SurfelInfoConstPtr & cell : cells )
            {
                if ( ! cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                num_invalid += !cell->m_surfel->valid_;
                num_less_pts += (cell->m_surfel->getNumPoints() < 10.0);
    //            if ( cell->m_surfel->getNumPoints() < 10.0 )
    //                LOG(1) <<"MarsMapWindow:: Cell: lvl: " << cell->m_level << " idx: " << Eigen::Map<const Eigen::Vector3i>(cell->m_index.data()).transpose() << " center: " << cell->m_center_s.transpose() << " pts: "<<  cell->m_surfel->getNumPoints();
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                ++num_valid;
            }
            LOG(1) << "MarsMapWindow:: Map-Cells: invalid: " << num_invalid << " less: " << num_less_pts << " valid: " << num_valid;
        }
#endif
    }
    else
    {
        // local map
        std::vector<MarsMapTypePtr> mapPtr;
        for ( int lvl = 0; lvl < asOneMap->m_map_params.m_num_levels; ++lvl )
        {
            //LOG(1) << "maps.size():" << maps.size() << " oldMaps.size(): " << oldMaps.size() << " lvl " << lvl << " param: " << oldMaps.back()->m_map_params.m_num_cells << " " << oldMaps.back()->m_map_params.m_size;
            maps.emplace_back(oldMaps.back());
            mapPtr.emplace_back(maps.back());
            oldMaps.pop_back();
        }
        //LOG(1) << "maps.size():" << maps.size() << " oldMaps.size(): " << oldMaps.size();

        constexpr bool id_map_orientation = true;
        bool shouldMoveMap = false;

        Eigen::Array3i minNumCells = Eigen::Array3i::Constant(std::numeric_limits<int>::max());
        {
            {
                ZoneScopedN("MarsMapWindow::AddCellsOneCloud");
                Eigen::Array3Xi numCells = Eigen::Array3Xi::Constant(3,asOneMap->m_map_params.m_num_levels,std::numeric_limits<int>::max());
#ifndef USE_TBB
                #pragma omp parallel for num_threads(OPT_THREAD_NUM)
                for ( int lvl = 0; lvl < asOneMap->m_map_params.m_num_levels; ++lvl )
                {
#else
                tbb::parallel_for(tbb::blocked_range<int>(0,asOneMap->m_map_params.m_num_levels,1), [&](const tbb::blocked_range<int> & range)
                {
                    for( int lvl=range.begin(); lvl!=range.end(); ++lvl)
                    {
#endif

                    Sophus::SE3d new_map_pose = map_pose;
                    Sophus::SE3d map_integration_pose;
                    MarsMapPointCloud::Ptr newCloudMap = MarsMapPointCloud::create(newCloud);

                    if ( m_empty )
                    {
                        map_integration_pose = Sophus::SE3d();
                        new_map_pose = map_pose;
                        if ( ! new_map_pose.unit_quaternion().isApprox(Eigen::Quaterniond::Identity()))
                            LOG(FATAL) << "I have an initial rotation? " << new_map_pose.params().transpose();
                    }
                    else
                        numCells.col(lvl) = asOneMap->getOriginCellShift( map_pose, new_map_pose, map_integration_pose, lvl );
                    if ( !map_integration_pose.params().isApprox(Sophus::SE3d().params()) )
                        newCloudMap->transform_points(map_integration_pose);

                    {
                        // clear down below after processing instead
                        //ZoneScopedN("MarsMapWindow::AddCellsOneCloud::clear_map_lvl");
                        //mapPtr[lvl]->clear();
                        //mapPtr[lvl] = MarsMapType::create(mapPtr[lvl]->m_map_params);
                    }
                    {
                        ZoneScopedN("MarsMapWindow::AddCellsOneCloud::set_cloud_lvl");
                        mapPtr[lvl]->setCloud( newCloudMap, new_map_pose, map_integration_pose.inverse().translation() );
                    }
                    LOG(1) << "lvl=["<<lvl<<"] nc: " << numCells.col(lvl).transpose() << " mapPose: " <<  map_pose.params().transpose() << " newMapPose: " << new_map_pose.params().transpose() << " mapIntPose: " << map_integration_pose.params().transpose();
                }
#ifdef USE_TBB
            });
#endif
                minNumCells = numCells.rowwise().minCoeff();
            }
            shouldMoveMap = (minNumCells>=m_min_num_cells_for_moving_map).any() && minNumCells[0] != std::numeric_limits<int>::max();
            LOG(1) << "minNumCells: " << minNumCells.transpose() <<  " th: " << m_min_num_cells_for_moving_map << " shouldMove: " << shouldMoveMap;
        }
        if ( shouldMoveMap )
        {
            ZoneScopedN("MarsMapWindow::MovingAsOneMap");
            int max_translation_lvl = 0;
            Sophus::SE3d max_translation_lvl_pose_inv;
            {
                ZoneScopedN("MarsMapWindow::MovingAsOneMap::Clear");
                asOneTmpMap->clear();
                local_map_moved = true;
                if ( id_map_orientation )
                {
                    local_map_last_moving_pose = mapPtr[max_translation_lvl]->m_pose_w;
                    local_map_pose.translation() += mapPtr[max_translation_lvl]->m_pose_w.translation();
                    max_translation_lvl_pose_inv.translation() = -mapPtr[max_translation_lvl]->m_pose_w.translation();
                    for ( int mapIdx = 0; mapIdx < int(maps.size()); ++mapIdx )
                        mapPtr[mapIdx]->transform_shift(max_translation_lvl_pose_inv.translation());
                    asOneMap->transform_shift(max_translation_lvl_pose_inv.translation());
                }
                else
                {
                    local_map_last_moving_pose = mapPtr[max_translation_lvl]->m_pose_w;
                    local_map_pose = local_map_pose * mapPtr[max_translation_lvl]->m_pose_w;
                    max_translation_lvl_pose_inv = mapPtr[max_translation_lvl]->m_pose_w.inverse();
                    for ( int mapIdx = 0; mapIdx < int(maps.size()); ++mapIdx )
                        mapPtr[mapIdx]->transform(max_translation_lvl_pose_inv);
                    asOneMap->transform(max_translation_lvl_pose_inv);
                }

                LOG(1) << "lm_last_moving_pose: " << local_map_last_moving_pose.params().transpose() << " lmp: " << local_map_pose.params().transpose()
                          << " mtlpi: " << max_translation_lvl_pose_inv.params().transpose() << " aOM: " << asOneMap->m_pose_w.params().transpose();
                for ( int mapIdx = 0; mapIdx < int(maps.size()); ++mapIdx )
                    LOG(1) << " M"<<mapIdx <<": " <<  maps[mapIdx]->m_pose_w.params().transpose();
                for ( int mapIdx = 0; mapIdx < int(mapPtr.size()); ++mapIdx )
                    LOG(1) << " M"<<mapIdx <<": " <<  mapPtr[mapIdx]->m_pose_w.params().transpose();
            }
            {
                ZoneScopedN("MarsMapWindow::MovingAsOneMap::addOneMap");
                SurfelInfoVector surfels;
                Sophus::SE3d pose_w;
                {
                    ZoneScopedN("MarsMapWindow::MovingAsOneMap::addOneMap::getCellsScanWise");
                    asOneMap->getCellsScanWise ( surfels, pose_w, false );
                }
                {
                    ZoneScopedN("MarsMapWindow::MovingAsOneMap::addOneMap::addCells");
                    if ( id_map_orientation )
                        asOneTmpMap->addCellsOnGrid ( surfels, pose_w.translation(), false ); // will be updated when adding last of the current ones.
                    else
                        asOneTmpMap->addCells ( surfels, pose_w, false ); // will be updated when adding last of the current ones.
                }
                asOneMap.swap(asOneTmpMap);
            }
        }
#ifdef USE_EASY_PBR
        constexpr bool showSurfelsScanWise = false;
        if ( showSurfelsScanWise )
        {
            SurfelInfoConstPtrVector surfels;
            Sophus::SE3d pose_w;
            for ( int mapIdx = 0; mapIdx < int(maps.size()); ++mapIdx )
            {
                const int lvl = mapIdx % asOneMap->m_map_params.m_num_levels;
                if ( lvl < 2 ) continue;
                //maps[mapIdx]->getCellsScanWise ( surfels, pose_w, false, lvl );
                maps[mapIdx]->getCellsFused ( surfels, pose_w, false );
                //LOG(1) << "Adding Cells: " << surfels.size() << " lvl :" << lvl << " id: " << maps[mapIdx]->m_id;
                //asOneMap->addCells ( surfels, pose_w, false, lvl ); // will be updated when adding last of the current ones.

                Eigen::Vector3f rc = VisMesh::getPlasma( float(mapIdx) / maps.size() );
                VisMesh v;
                for ( const SurfelInfoConstPtr & cell : surfels )
                {
                    if ( cell == nullptr || !cell->m_surfel ) continue;
                    if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                    v.addSurfelTriangleEllipsoid(cell->m_surfel->eigen_vectors_.cast<float>(), cell->m_surfel->eigen_values_.cast<float>(), (pose_w * (cell->m_center_s + cell->m_surfel->mean_).cast<double>()).template cast<float>(), rc);
                }
                v.showSurfelTriangleEllipsoidMesh("marsSurfels"+std::to_string(mapIdx), true, local_map_pose);
            }
        }
#endif
        LOG(1) << "showed status in between! #pts: " << asOneMap->m_num_points;

        {
            ZoneScopedN("MarsMapWindow::AddedCells::AddSurfelToAsOneMap");
            SurfelInfoVector surfels;
            Sophus::SE3d pose_w;
            for ( int lvl = 0; lvl < asOneMap->m_map_params.m_num_levels; ++lvl )
            {
                mapPtr[lvl]->getCellsScanWise ( surfels, pose_w, false );

                LOG(1) << "Adding Cells: " << surfels.size() << " lvl :" << lvl << " param: " << mapPtr[lvl]->m_map_params.m_num_cells << " " << mapPtr[lvl]->m_map_params.m_size;
                if ( shouldMoveMap ) LOG(1) << "pose_w: " << pose_w.params().transpose() << " oneMap: "<< asOneMap->m_pose_w.params().transpose();

                if ( id_map_orientation )
                {
                    if ( m_empty && lvl==0 )
                        asOneMap->setCellsOnGrid ( surfels, pose_w.translation(), lvl == (asOneMap->m_map_params.m_num_levels-1), lvl );
                    else
                        asOneMap->addCellsOnGrid ( surfels, pose_w.translation(), lvl == (asOneMap->m_map_params.m_num_levels-1), lvl );
                }
                else
                {
                    if ( m_empty && lvl==0 )
                        asOneMap->setCells ( surfels, pose_w, lvl == (asOneMap->m_map_params.m_num_levels-1), lvl );
                    else
                        asOneMap->addCells ( surfels, pose_w, lvl == (asOneMap->m_map_params.m_num_levels-1), lvl );
                }
                m_empty = false;
            }
        }
        asOneMap->update(m_use_adaptive);
        LOG(1) << "showed status after adding! #pts: " << asOneMap->m_num_points;
#ifdef USE_EASY_PBR
        {
            ZoneScopedN("MarsMapWindow::AddedCells::getCells");
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector cells;
            asOneMap->getCells( cells, pose );
            int num_invalid = 0, num_less_pts = 0, num_valid = 0;
            for ( const SurfelInfoConstPtr & cell : cells )
            {
                if ( !cell->m_surfel ) continue;
                if ( !cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend! " << cell->m_level << " idx: " << Eigen::Map<const Eigen::Vector3i>(cell->m_index.data()).transpose() << " center: " << cell->m_center_s.transpose() << " pts: "<<  cell->m_surfel->getNumPoints();
                num_invalid += !cell->m_surfel->valid_;
                num_less_pts += (cell->m_surfel->getNumPoints() < 10.0);
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                ++num_valid;
            }
            LOG(1) << "MarsMapWindow:: AsOneMap-Cells: invalid: " << num_invalid << " less: " << num_less_pts << " valid: " << num_valid;
            num_invalid = 0, num_less_pts = 0, num_valid = 0;
            for ( int lvl = 0; lvl < asOneMap->m_map_params.m_num_levels; ++lvl )
            {
                mapPtr[lvl]->getCells( cells, pose );
                for ( const SurfelInfoConstPtr & cell : cells )
                {
                    if ( !cell->m_surfel ) continue;
                    if ( !cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    num_invalid += !cell->m_surfel->valid_;
                    num_less_pts += (cell->m_surfel->getNumPoints() < 10.0);
                    if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                    ++num_valid;
                }
            }
            LOG(1) << "oneMapPose: " <<  map_pose.params().transpose() << " asOneMapPose: " << asOneMap->getMapPose().params().transpose();
            LOG(1) << "MarsMapWindow::LvlCells: invalid: " << num_invalid << " less: " << num_less_pts << " valid: " << num_valid;
        }
#endif
        for ( int lvl = 0; lvl < asOneMap->m_map_params.m_num_levels; ++lvl )
        {
            mapPtr[lvl]->clear();
        }
    }
    //LOG(1) << "cloud size: " << clouds.size() << " " << maps.size() << " qs: " << queue_size;
}

MarsMapTypePtr MarsMapWindow::getAsOneCenteredMap()
{
    return asOneMap;
}
