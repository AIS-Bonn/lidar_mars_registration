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

#include "MarsSplineRegistrator.h"
#include <Eigen/Geometry>

#define LOGURU_REPLACE_GLOG 1
#include "loguru.hpp"
//configuru
#define CONFIGURU_WITH_EIGEN 1
#define CONFIGURU_IMPLICIT_CONVERSIONS 1
#include <configuru.hpp>
using namespace configuru;

//boost
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "MarsPointTypes.h"
#include "MarsMap.h"
#include "MarsMapWindow.h"

#include "basalt/spline/se3_spline.h"
#include "MarsSplineRegistration.h"

#include "StopWatch.h"

#include <Tracy.hpp>

#ifdef USE_EASY_PBR
#include "easy_pbr/Mesh.h"
#include "easy_pbr/Scene.h"
#include "easy_pbr/LabelMngr.h"
#include "VisUtils.h"
#else
#include "ros/ros.h"
#endif

MarsSplineRegistrator::MarsSplineRegistrator(const std::string & config_file)
{
    init_params(config_file);
}

void MarsSplineRegistrator::init_params(const std::string & config_file){

    //read all the parameters
    std::string config_file_abs;
    if (fs::path(config_file).is_relative()){
        config_file_abs=(fs::path(PROJECT_SOURCE_DIR) / config_file).string();
    }else{
        config_file_abs=config_file;
    }
    Config cfg = configuru::parse_file(config_file_abs, CFG);
    Config dyn_config=cfg["registration_adaptor"];
    m_scene_map_window = dyn_config["scene_map_window"];
    m_local_map_window = dyn_config["local_map_window"];

    m_min_range = dyn_config["min_range"];
    m_min_initial_scans = dyn_config["min_initial_scans"];
    m_max_iterations = dyn_config["max_iterations"];

    m_keyframe_use_rotation = dyn_config["keyframe_use_rotation"];
    m_min_keyframe_rotation = double(dyn_config["min_keyframe_rotation"]) * M_PI/180.;
    m_min_keyframe_distance = dyn_config["min_keyframe_distance"];
    m_use_closest_kf = dyn_config["use_closest_keyframe"];
    m_use_adaptive = dyn_config["use_adaptive"];
    bool sca_use_constraint = dyn_config["sca_use_constraint"];
    double sca_constraint_factor = dyn_config["sca_constraint_factor"];

    MapParameters::plane_scale_factor = dyn_config["plane_scale_factor"];
    MapParameters::first_ev_planarity_threshold = dyn_config["first_ev_planarity_threshold"];
    MapParameters::second_ev_degenerate_threshold = dyn_config["second_ev_degenerate_threshold"];

    m_interpolationTimeOffset = dyn_config["interpolation_time_offset"];

    m_init_spline_window = dyn_config["init_spline_window"];

    m_surfel_type = SurfelType((int)dyn_config["surfel_type"]);

#ifdef USE_EASY_PBR
    m_show_mars_map = dyn_config["show_mars_map"];
    m_show_scene_map = dyn_config["show_scene_map"];
    m_show_large_errors = dyn_config["show_large_errors"];
    m_show_scene_surfels = dyn_config["show_scene_surfels"];
    m_show_mars_surfels = dyn_config["show_mars_surfels"];
    m_show_initial_surfels = dyn_config["show_initial_surfels"];
    m_show_subdiv_surfels = dyn_config["show_subdiv_surfels"];
    m_show_knots = dyn_config["show_knots"];
    m_show_scene_knots = dyn_config["show_scene_knots"];
    m_show_opt_scene = dyn_config["show_opt_scene"];
#endif

    if ( m_interpolationTimeOffset < 0.) LOG(FATAL) << "Shall not interpolate backwards in time.";

    m_scan_id = 0;
    m_firstMarsMesh = true;
    m_localMarsMap = MarsMapWindow::create( m_local_map_window, false, true );
    m_sceneMarsMap = MarsMapWindow::create( m_scene_map_window, true, m_use_adaptive );
    m_marsRegistrator = MarsSplineRegistration::create( );
    m_marsRegistrator->setUseAdaptive ( m_use_adaptive );
    m_marsRegistrator->setScaConstraint ( sca_use_constraint, sca_constraint_factor );

    const double dt = 0.1;
    m_trajectory_spline = std::make_shared<MarsSplineType>( TimeConversion::to_ns(dt), 0 );
    m_trajectory_spline->setKnots(Sophus::SE3d(), MarsSplineType::N );
    m_mars_map_mesh = MarsMapPointCloud::create();
    m_mars_local_map_mesh = MarsMapPointCloud::create();
    m_mars_scene_map_mesh = MarsMapPointCloud::create();

#ifdef USE_EASY_PBR
    if constexpr ( has_semantics_v<MarsMapPointCloud> )
    {
    Config mngr_config = cfg["label_mngr"];
    m_label_manager = std::make_shared<easy_pbr::LabelMngr>(mngr_config);
    }
#endif

}

MarsSplineRegistrator::~MarsSplineRegistrator()
{
}

#ifdef USE_EASY_PBR
template< typename T>
typename VisMesh::Ptr showSemanticColored( typename T::Ptr scene_cloud, const Sophus::SE3d & interp_scene_pose )
{
    VisMesh::Ptr vm = VisMesh::create();
    if constexpr ( has_semantics_v<T> )
    {
        const Eigen::Matrix3Xf & pts = scene_cloud->toMatrix3Xf();
        const Eigen::VectorXf & pt_class = scene_cloud->m_class;
        const Eigen::MatrixXf & pt_prob = scene_cloud->m_class;
        const Eigen::VectorXf & intensity = scene_cloud->intensity();
        const Eigen::VectorXu & reflectivity = scene_cloud->reflectivity();
        const int num_pts = scene_cloud->size();
        for ( int i = 0; i < num_pts; ++i )
        {
            if ( pt_prob.col(i).isApproxToConstant(pt_prob.col(i)(0),1e-10) ) continue;
            int idx = 0;
            pt_prob.col(i).maxCoeff(&idx);
            //LOG(1) << "s["<<i<<"]: " << pt_prob.col(i).transpose() << " " << idx;
            vm->addPoint( (interp_scene_pose * pts.col(i).cast<double>()).cast<float>(), Eigen::Vector3f::Constant(idx), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr );
        }
    }
    return vm;
}
#endif

// cur_pose is in field frame
Sophus::SE3d MarsSplineRegistrator::register_cloud ( MarsMapPointCloudPtr sceneCloud )
{
    ++m_scan_id;
    m_was_keyframe = false;
    constexpr bool semantics_only = false && has_semantics_v<MarsMapPointCloud>;

    static Sophus::SE3d new_scene_origin_to_newest_local_scan = Sophus::SE3d();
    static Sophus::SE3d last_interp_pose;
    static Sophus::SE3d prev_scene_pose;
    static Sophus::SE3d last_pose = m_cur_pose;
    static Sophus::SE3d last_pos = m_cur_pos;

    static Eigen::VectorXt oldTimes;

    //LOG(1) << "sceneCloud: " << sceneCloud->size() << " " << sceneCloud->m_points.cols() << " id:" << sceneCloud->m_scan_id;

    if ( ! m_sceneMarsMap ) LOG(FATAL) << "no scene map.";
    if ( ! m_localMarsMap ) LOG(FATAL) << "no local map.";

    static uint64_t firstTime = 0;
    static uint64_t last_time = 0;
    if ( m_firstMarsMesh )
    {
        if ( m_scan_id == 1 )
        {
            m_localMarsMap->local_map_pose = Sophus::SE3d();
            firstTime = sceneCloud->m_stamp;
        }
        m_firstMarsMesh = m_scan_id < (m_min_initial_scans-1);
        static MarsMapPointCloudPtr initial_cloud = nullptr;
        static Sophus::SE3d initial_pose;
        Sophus::SE3d pose;
        if ( m_min_initial_scans == 1 )
        {
            initial_cloud = sceneCloud;
            initial_pose = Sophus::SE3d();
            pose = Sophus::SE3d();
        }
        else
        {
            if ( initial_cloud == nullptr)
            {
                initial_cloud = MarsMapPointCloud::create(sceneCloud);
                initial_pose = Sophus::SE3d();
                pose = Sophus::SE3d();
            }
            else
            {
                // fuse initial clouds.
                initial_cloud->add( sceneCloud, m_localMarsMap->local_map_pose.inverse()*m_cur_pose );
                pose = m_localMarsMap->local_map_pose.inverse()*m_cur_pose;
            }
        }
        if ( ! m_firstMarsMesh )
        {
#ifdef USE_EASY_PBR
            {
                VisMesh v;
                Eigen::Vector3f rc = VisMesh::randColor();
                const Eigen::Matrix3Xf & pts = initial_cloud->toMatrix3Xf();
                const Eigen::VectorXf & intensity = initial_cloud->intensity();
                const Eigen::VectorXu & reflectivity = initial_cloud->reflectivity();
                const int num_pts = initial_cloud->size();
                for ( int idx = 0; idx < num_pts; ++idx)
                {
                    //if ( pts.col(idx).allFinite() )
                    v.addPoint(pts.col(idx),rc, intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                }
                v.showPointMesh("initial_cloud",true,m_localMarsMap->local_map_pose);
            }
#endif
            {
                std::unique_lock<std::mutex> lock(m_localMapMutex);
                m_localMarsMap->addCloud(initial_cloud, initial_pose);
            }
#ifdef USE_EASY_PBR
            if ( m_show_initial_surfels )
            {
                MarsMapType::Ptr asOneCenteredLocalMap = m_localMarsMap->getAsOneCenteredMap();
                for ( int lvl = 0; lvl < asOneCenteredLocalMap->m_map_params.m_num_levels; ++lvl)
                {
                    VisMesh vl;
                    SurfelInfoVector cells;
                    asOneCenteredLocalMap->getCellsOnLevel ( lvl, cells, false);
                    vl.reserveSurfelTriangleEllipsoid(cells.size());
                    for ( const SurfelInfo & cell : cells )
                    {
                        if ( ! cell.m_surfel ) continue;
                        if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                        if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 ) continue;
                        vl.addSurfel(cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, cell.m_center_s + cell.m_surfel->mean_, VisMesh::green, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(lvl)/2, m_surfel_type);
                    }
                    vl.showSurfelMesh("marsSurfels"+std::to_string(lvl), true, m_localMarsMap->local_map_pose, m_surfel_type);
                }
                std::vector<VisMesh> vls(asOneCenteredLocalMap->m_map_params.m_num_levels);
                std::vector<VisMesh> vln(asOneCenteredLocalMap->m_map_params.m_num_levels);
                Sophus::SE3d pose;
                SurfelInfoConstPtrVector cells;
                asOneCenteredLocalMap->getCells (cells, pose );
                for ( int lvl = 0; lvl < asOneCenteredLocalMap->m_map_params.m_num_levels; ++lvl)
                {
                    vls[lvl].reserveSurfelTriangleEllipsoid(cells.size());
                    vln[lvl].reserveSurfelTriangleEllipsoid(cells.size());
                }
                for ( const SurfelInfoConstPtr & cell : cells )
                {
                    if ( ! cell->m_surfel ) continue;
                    if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    if ( ! cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                    vls[cell->m_level].addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, VisMesh::green, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                    vln[cell->m_level].addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, VisMesh::colorFromNormal(cell->m_surfel->normal_.normalized()), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                }
                for ( int lvl = 0; lvl < asOneCenteredLocalMap->m_map_params.m_num_levels; ++lvl)
                {
                    vls[lvl].showSurfelMesh("marsSurfelCells"+std::to_string(lvl), true, m_localMarsMap->local_map_pose,m_surfel_type);
                    vln[lvl].showSurfelMesh("marsSurfelNormalCells"+std::to_string(lvl), true, m_localMarsMap->local_map_pose,m_surfel_type);
                }
            }
#endif
        }

        last_time = sceneCloud->m_stamp;
        m_scene_poses.clear();
        m_scene_poses.emplace_back(pose);
        m_interpolated_scene_pose = pose;
        m_last_key_frame_stamp = sceneCloud->m_stamp;
        return pose;
    }

    // update estimate with diff to last:
    {
        ZoneScopedN("MarsSplineRegistrator::addSceneCloud");
        const Sophus::SE3d new_cur_pose = m_localMarsMap->local_map_pose.inverse() * m_cur_pose;
        std::unique_lock<std::mutex> lock(m_sceneMapMutex);
        m_sceneMarsMap->addCloud(sceneCloud, new_cur_pose);
        //LOG(1) << "lmpi: " <<  m_localMarsMap->local_map_pose.inverse().params().transpose() << " cp: " <<  cur_pose.params().transpose() << " ncp: " << new_cur_pose.params().transpose();
        std::vector<MarsMap*> sceneMapPtrVec = m_sceneMarsMap->getMapPtrVec( ); // sceneMaps now in frame of local map
        while ( m_scene_surfels.size() >= sceneMapPtrVec.size() )
            m_scene_surfels.pop_front();

        m_scene_surfels.emplace_back();
        Sophus::SE3d pose;
        sceneMapPtrVec.back()->getCells(m_scene_surfels.back(), pose);
    }

    Eigen::VectorXt curTimes = m_sceneMarsMap->getTimesSince( last_time );
    if ( oldTimes.size() > 0 )
    {
        const int64_t maxTime = curTimes.tail<1>()[0]+1+TimeConversion::to_ns(m_interpolationTimeOffset);

        MarsSplineTypePtr interpolated_spline = std::make_shared<MarsSplineType>( maxTime, 0 );
        interpolated_spline->setKnots(Sophus::SE3d(), MarsSplineType::N );

        std::vector<Sophus::SE3d> poses;
        for ( int cloudIdx = 0; cloudIdx < oldTimes.rows(); ++cloudIdx)
        {
            const Sophus::SE3d pre_interpolated_pose_se3 = m_trajectory_spline->pose(oldTimes(cloudIdx) );
            poses.emplace_back(pre_interpolated_pose_se3);
        }

        for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
        {
            interpolated_spline->setKnot(poses.front(),knot_idx);
        }

        //orig_diff_to_last
        Eigen::VectorXt poseTimes = Eigen::VectorXt::Zero(oldTimes.rows(),1);
        {
            if ( oldTimes.size() < int(m_sceneMarsMap->queue_size) )
                poseTimes = curTimes.head(curTimes.rows()-1);
            else
            {
                constexpr bool m_use_const_vel_for_init = true;
                if ( m_use_const_vel_for_init )
                {
                    poseTimes = Eigen::VectorXt::Zero(oldTimes.rows()+1,1);
                    poseTimes.tail(curTimes.rows()) = curTimes;

                    Eigen::Matrix3Xd vs_w = Eigen::Matrix3Xd::Zero(3,oldTimes.rows());
                    for ( int idx = 0; idx < oldTimes.size(); ++idx )
                    {
                        const int64_t ts ( oldTimes(idx) );
                        vs_w.col(idx) = m_trajectory_spline->transVelWorld(ts);
                    }

                    const int64_t ts ( oldTimes.tail<1>()(0) );
                    const int64_t last_ts ( curTimes.tail<1>()(0) );
                    const int64_t second_ts ( curTimes.tail<2>()(0) );
                    const double dt = 1e-9 * ( last_ts - second_ts );
                    // get velocity at last pose:
                    const Sophus::SE3d P_wi = m_trajectory_spline->pose(ts);
                    const Eigen::Vector3d r_i = m_trajectory_spline->rotVelBody(ts);
                    //const Eigen::Vector3d a_w = m_trajectory_spline->transAccelWorld(ts);
                    //const Eigen::Vector3d v_w = m_trajectory_spline->transVelWorld(ts);
                    //const Eigen::Vector3d v_w = vs_w.cols() > 1 ? vs_w.col(vs_w.cols()-2) : vs_w.col(0);
                    const Eigen::Vector3d v_w = vs_w.cols() > 1 ? Eigen::Vector3d(vs_w.rowwise().mean()) : vs_w.col(0);

                    //const Sophus::SO3d R_wi = P_wi.so3() * Sophus::SO3d::exp( .5 * dt * r_i );
                    //const Eigen::Vector3d a_wi = R_wi * a_w;

                    Sophus::SE3d new_pose_wi;
                    new_pose_wi.translation() = P_wi.translation() + v_w * dt; // + 0.5 * a_wi * dt * dt;
                    if ( m_rotation_prior_valid )
                        new_pose_wi.so3() = P_wi.so3() * m_rotation_prior;
                    else
                        new_pose_wi.so3() = P_wi.so3() * Sophus::SO3d::exp( dt * r_i );

                    //LOG(1) << "r: " << r_i.transpose() //<< " a: " << a_w.transpose()
                    //          << " v: " << v_w.transpose() << " dt: " << dt << " t: " << P_wi.translation().transpose();
                    //LOG(1) << "|vs_w|: "<< vs_w.colwise().norm();
                    //LOG(1) << "PreInterp=["<<poseTimes.rows()<<"]=["<< new_pose_wi.translation().transpose()<<"]";

                    poses.emplace_back( new_pose_wi );
                }
                else
                    poseTimes.tail(curTimes.rows()-1) = curTimes.head(curTimes.rows()-1);
            }
        }

        //LOG(1) << "numPoses: " << poses.size() << " poseTimes: " << poseTimes.transpose() << " oldTimes: " << oldTimes.transpose() << " curTimes: " << curTimes.transpose();
        m_marsRegistrator->optimizeSplineKnotsFromPoses ( interpolated_spline, poseTimes, poses );
        for ( int cloudIdx = 0; cloudIdx < poseTimes.rows(); ++cloudIdx)
        {
            const int64_t t = poseTimes(cloudIdx);
            const Sophus::SE3d ref_pose = poses[cloudIdx];
            const Sophus::SE3d interpolated_spline_pose = interpolated_spline->pose( t );
            const Sophus::Vector3d trans_error = interpolated_spline->positionResidual(t,ref_pose.translation());
            const Sophus::Vector3d rot_error = interpolated_spline->orientationResidual(t,ref_pose.so3());
            if ( rot_error.norm() > 0.1 || trans_error.norm() > 0.05 )
              LOG(WARNING) << "PreInterp["<<cloudIdx<<"]=["<<interpolated_spline_pose.params().transpose()<<"] diff: " << rot_error.transpose() << " (| "<< rot_error.norm()<< " |)" << trans_error.transpose() << " (|"<< trans_error.norm()<<" |)";
        }
        m_trajectory_spline = interpolated_spline;
    }
    else
    {
        //LOG(1) << "Not enough entries yet.";
        if ( curTimes.rows() == 1 && oldTimes.rows() == 0 ) // just initialized, now set to estimated diff
        {
            const int64_t maxTime = curTimes.tail<1>()[0] + 1 + TimeConversion::to_ns(m_interpolationTimeOffset);
            m_trajectory_spline = std::make_shared<MarsSplineType>( maxTime, 0 );
            m_trajectory_spline->setKnots( Sophus::SE3d(), MarsSplineType::N );
            for (size_t i = 0; i < m_trajectory_spline->numKnots(); ++i) {
                m_trajectory_spline->setKnot(new_scene_origin_to_newest_local_scan,i);
            }
        }
    }

    oldTimes = curTimes;

    if ( m_sceneMarsMap->hasOutOfWindow() )
    {
        LOG(FATAL) << "This should never happen!";
    }

    // Assumptions:
    // -> each cloud is originally centered at Id in its own grid
    // -> center_s is in sensors coord system, and pose_w is needed for transformation: center_w = pose_w * center_s
    // -> youngest cloud in localMarsMap is kept at Id. others are not.

    MarsMapType::Ptr asOneCenteredLocalMap = m_localMarsMap->getAsOneCenteredMap();

    Eigen::VectorXt sceneTimes;
    std::vector<MarsMap*> sceneMapPtrVec;
    {
        sceneTimes = curTimes;
        sceneMapPtrVec = m_sceneMarsMap->getTransformedMapPtrVec( prev_scene_pose ); // sceneMaps now in frame of local map
    }
#ifdef USE_EASY_PBR
    if ( m_show_knots )
    {
        VisMesh::Ptr tm = VisMesh::create();
        VisMesh::Ptr bm = VisMesh::create();
        VisMesh::Ptr km = VisMesh::create();
        VisMesh::Ptr sm = VisMesh::create();
        Eigen::Matrix3Xd scene_translation = Eigen::Matrix3Xd::Zero(3,m_scene_poses.size());
        Eigen::Matrix3Xd interp_translation = Eigen::Matrix3Xd::Zero(3,sceneTimes.rows());
        Eigen::Matrix3Xd knot_translation = Eigen::Matrix3Xd::Zero(3,m_trajectory_spline->numKnots());
        for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
        {
            const Sophus::SE3d knot = m_trajectory_spline->getKnot(knot_idx);
            km->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation()).cast<float>(), VisMesh::blue);
            knot_translation.col(knot_idx) = (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation());
        }
        for ( size_t i = 0; i < m_scene_poses.size(); ++i )
        {
            const Sophus::SE3d & scan_pose = m_scene_poses[i];
            bm->addPoint((m_localMarsMap->local_map_pose * scan_pose).translation().cast<float>(), VisMesh::green);
            scene_translation.col(i) = m_localMarsMap->local_map_pose * scan_pose.translation();
        }
        for ( int i = 0; i < sceneTimes.rows(); ++i)
        {
            const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( sceneTimes[i] ); // curTimes still contains second
            sm->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation()).cast<float>(), VisMesh::red);
            interp_translation.col(i) = (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation());
        }
        constexpr uint64_t num_samples = 1000;
        const uint64_t dt = m_trajectory_spline->getDtNs();
        uint64_t dt_inc = dt / num_samples;
        Sophus::SE3d prev_interp_scene_pose;
        bool prev_is_set = false;
        Eigen::Matrix3Xd points = Eigen::Matrix3Xd::Zero(3,num_samples/50);
        for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
        {
            Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( m_trajectory_spline->minTimeNs() + dt_inc * t_inc ); // curTimes still contains second
            Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            if ( prev_is_set )
                tm->addEdge( prev_interp_scene_pose.translation().cast<float>(), interp_scene_pose.translation().cast<float>(), VisMesh::getPlasma( float(t_inc) / num_samples ) );
            prev_interp_scene_pose = interp_scene_pose;
            prev_is_set = true;
            if ( t_inc % 50 == 0 )
                points.col(t_inc/50) = interp_scene_pose.translation();
        }
        LOG(1) << "lm_t: " << m_localMarsMap->local_map_pose.translation().transpose();
        LOG(1) << "ps_t: " << prev_scene_pose.translation().transpose();
        LOG(1) << "sceneTimes: " << sceneTimes.transpose();
        LOG(1) << "knot_t:\n" << knot_translation;
        LOG(1) << "scene_t:\n" << scene_translation;
        LOG(1) << "interp_t:\n" << interp_translation;
        LOG(1) << "\n" << points;
        km->showPointMesh("InitKnotPoses"+std::to_string(m_scan_id),true);
        bm->showPointMesh("InitScanPose"+std::to_string(m_scan_id),true);
        sm->showPointMesh("InitInterpPose"+std::to_string(m_scan_id),true);
        tm->showEdgeMesh("InitSpline"+std::to_string(m_scan_id),true);
    }
#endif
#ifdef USE_EASY_PBR
    if ( m_show_opt_scene)
    {
        VisMesh::Ptr vm = VisMesh::create();
        std::vector<MarsMapPointCloud::Ptr> scene_clouds = m_sceneMarsMap->getCloudSharedPtrVec();
        for ( size_t cloud_idx = 0; cloud_idx < scene_clouds.size(); ++cloud_idx)
        {
            VisMesh::Ptr cm = VisMesh::create();
            Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( curTimes(cloud_idx) ); // curTimes still contains second
            Sophus::SE3d interp_scene_pose  = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            MarsMapPointCloud::Ptr scene_cloud = scene_clouds[cloud_idx];
            const Eigen::Matrix3Xf & pts = scene_cloud->toMatrix3Xf();
            const Eigen::VectorXf & intensity = scene_cloud->intensity();
            const Eigen::VectorXu & reflectivity = scene_cloud->reflectivity();
            const Eigen::VectorXu & scan_line = scene_cloud->m_scan_line_id;
            const int num_pts = scene_cloud->size();
            for ( int idx = 0; idx < num_pts; ++idx )
            {
                const Eigen::Vector3f pt = (interp_scene_pose * pts.col(idx).cast<double>()).template cast<float>();
                vm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                cm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
            }
            cm->showPointMesh("InitSceneCloud"+std::to_string(m_scan_id)+"_"+std::to_string(cloud_idx), true);
        }
        vm->showPointMesh("InitSceneMap"+std::to_string(m_scan_id), true);
    }
#endif

    bool successfulRegistration = false;
    {
        ZoneScopedN("MarsSplineRegistrator::registration");
        successfulRegistration = m_marsRegistrator->estimate ( asOneCenteredLocalMap.get(), sceneMapPtrVec, sceneTimes, m_trajectory_spline, m_max_iterations );
    }

    //LOG(1) << "Reg was succesful: " << successfulRegistration;

    Sophus::SE3d newest_cloud_pose = prev_scene_pose;
    {
        const Sophus::SE3d interp_scene_pose = m_trajectory_spline->pose( curTimes.tail<1>()[0] + TimeConversion::to_ns(m_interpolationTimeOffset) );
        m_interpolated_scene_pose = prev_scene_pose * interp_scene_pose;
        m_scene_poses.clear();
        m_scene_poses.reserve(curTimes.rows());
        for ( int row = 0; row < curTimes.rows(); ++row)
        {
            const Sophus::SE3d interp_pose = m_trajectory_spline->pose( curTimes[row] );
            //const Sophus::Vector6d interp_velo ( m_trajectory_spline->getVelocity( curTimes.tail<1>()[0] ).matrix() );
            const Sophus::SE3d cloud_pose = prev_scene_pose * interp_pose;
            m_scene_poses.emplace_back(cloud_pose);
            //m_scene_velocity.emplace_back( interp_velo );
        }
        newest_cloud_pose = m_scene_poses.back();

        Eigen::Matrix6d covar = Eigen::Matrix6d::Zero();
        const bool valid_covar = m_marsRegistrator->getCovariance( sceneMapPtrVec.size()-1, newest_cloud_pose, covar );
        if ( ! valid_covar ) covar.setIdentity();
        m_current_pose_cov = covar;
    }

//    for ( int cloudIdx = 0; cloudIdx < curTimes.rows(); ++cloudIdx)
//    {
//        Sophus::SE3d pre_interpolated_pose_se3 = m_trajectory_spline->pose(curTimes(cloudIdx) );
//        LOG(1) << "RegInterp["<<cloudIdx<<"]=["<<pre_interpolated_pose_se3.params().transpose()<<"]";
//    }


#ifdef USE_EASY_PBR
    if ( m_show_knots )
    {
        VisMesh::Ptr tm = VisMesh::create();
        VisMesh::Ptr bm = VisMesh::create();
        VisMesh::Ptr km = VisMesh::create();
        VisMesh::Ptr sm = VisMesh::create();
        Eigen::Matrix3Xd scene_translation = Eigen::Matrix3Xd::Zero(3,m_scene_poses.size());
        Eigen::Matrix3Xd interp_translation = Eigen::Matrix3Xd::Zero(3,sceneTimes.rows());
        Eigen::Matrix3Xd knot_translation = Eigen::Matrix3Xd::Zero(3,m_trajectory_spline->numKnots());
        for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
        {
            const Sophus::SE3d knot = m_trajectory_spline->getKnot(knot_idx);
            km->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation()).cast<float>(), VisMesh::blue);
            knot_translation.col(knot_idx) = (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation());
        }
        for ( size_t i = 0; i < m_scene_poses.size(); ++i )
        {
            const Sophus::SE3d & scan_pose = m_scene_poses[i];
            bm->addPoint((m_localMarsMap->local_map_pose * scan_pose.translation()).cast<float>(), VisMesh::green);
            scene_translation.col(i) = m_localMarsMap->local_map_pose * scan_pose.translation();
        }
        for ( int i = 0; i < sceneTimes.rows(); ++i)
        {
            const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( sceneTimes[i] ); // curTimes still contains second
            sm->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation()).cast<float>(), VisMesh::red);
            interp_translation.col(i) = (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation());
        }
        constexpr uint64_t num_samples = 1000;
        const uint64_t dt = m_trajectory_spline->getDtNs();
        uint64_t dt_inc = dt / num_samples;
        Sophus::SE3d prev_interp_scene_pose;
        bool prev_is_set = false;
        Eigen::Matrix3Xd points = Eigen::Matrix3Xd::Zero(3,num_samples/50);
        for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
        {
            Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( m_trajectory_spline->minTimeNs() + dt_inc * t_inc ); // curTimes still contains second
            Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            if ( prev_is_set )
                tm->addEdge( prev_interp_scene_pose.translation().cast<float>(), interp_scene_pose.translation().cast<float>(), VisMesh::getPlasma( float(t_inc) / num_samples ) );
            prev_interp_scene_pose = interp_scene_pose;
            prev_is_set = true;
            if ( t_inc % 50 == 0 )
                points.col(t_inc/50) = interp_scene_pose.translation();
        }
        //LOG(1) << "sceneTimes: " << sceneTimes.transpose();
        //LOG(1) << "knot_t:\n" << knot_translation;
        //LOG(1) << "scene_t:\n" << scene_translation;
        //LOG(1) << "interp_t:\n" << interp_translation;
        //LOG(1) << "\n" << points;
        km->showPointMesh("KnotPoses"+std::to_string(m_scan_id),true);
        bm->showPointMesh("ScanPose"+std::to_string(m_scan_id),true);
        sm->showPointMesh("InterpPose"+std::to_string(m_scan_id),true);
        tm->showEdgeMesh("Spline"+std::to_string(m_scan_id),true);
    }
#endif
#ifdef USE_EASY_PBR
    if ( m_show_scene_knots )
    {
        VisMesh::Ptr tm = VisMesh::create();
        VisMesh::Ptr bm = VisMesh::create();
        VisMesh::Ptr km = VisMesh::create();
        VisMesh::Ptr sm = VisMesh::create();
        Eigen::Matrix3Xd scene_translation = Eigen::Matrix3Xd::Zero(3,m_scene_poses.size());
        Eigen::Matrix3Xd interp_translation = Eigen::Matrix3Xd::Zero(3,sceneTimes.rows());
        Eigen::Matrix3Xd knot_translation = Eigen::Matrix3Xd::Zero(3,m_trajectory_spline->numKnots());
        for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
        {
            const Sophus::SE3d knot = m_trajectory_spline->getKnot(knot_idx);
            km->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation()).cast<float>(), VisMesh::blue);
            knot_translation.col(knot_idx) = (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation());
        }
        for ( size_t i = 0; i < m_scene_poses.size(); ++i )
        {
            const Sophus::SE3d & scan_pose = m_scene_poses[i];
            bm->addPoint((m_localMarsMap->local_map_pose * scan_pose.translation()).cast<float>(), VisMesh::green);
            scene_translation.col(i) = m_localMarsMap->local_map_pose * scan_pose.translation();
        }
        for ( int i = 0; i < sceneTimes.rows(); ++i)
        {
            const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( sceneTimes[i] ); // curTimes still contains second
            sm->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation()).cast<float>(), VisMesh::red);
            interp_translation.col(i) = (m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose.translation());
        }
        constexpr uint64_t num_samples = 1000;
        const uint64_t dt = m_trajectory_spline->getDtNs();
        uint64_t dt_inc = dt / num_samples;
        Sophus::SE3d prev_interp_scene_pose;
        bool prev_is_set = false;
        Eigen::Matrix3Xd points = Eigen::Matrix3Xd::Zero(3,num_samples/50);
        for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
        {
            Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( m_trajectory_spline->minTimeNs() + dt_inc * t_inc ); // curTimes still contains second
            Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            if ( prev_is_set )
                tm->addEdge( prev_interp_scene_pose.translation().cast<float>(), interp_scene_pose.translation().cast<float>(), VisMesh::getPlasma( float(t_inc) / num_samples ) );
            prev_interp_scene_pose = interp_scene_pose;
            prev_is_set = true;
            if ( t_inc % 50 == 0 )
                points.col(t_inc/50) = interp_scene_pose.translation();
        }
        km->showPointMesh("SceneKnotPoses",true);
        bm->showPointMesh("SceneScanPose",true);
        sm->showPointMesh("SceneInterpPose",true);
        tm->showEdgeMesh("SceneSpline",true);
    }
#endif

#ifdef USE_EASY_PBR
    constexpr bool computeKappa = false;
    if constexpr( computeKappa )
    {
        Eigen::VectorXt newCurTimes = sceneTimes;
        for ( size_t idx = 0; idx < sceneMapPtrVec.size(); ++idx)
        {
            MarsMap * scene = sceneMapPtrVec[idx];
            const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( newCurTimes(idx) ); // curTimes still contains second
            const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            Eigen::Matrix3f a = Eigen::Matrix3f::Zero(), an = Eigen::Matrix3f::Zero();
            int n = 0; double wn = 0;
            const Eigen::Vector3f rc = VisMesh::getPlasma(float(idx+1) / sceneMapPtrVec.size());
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            scene->getCells( scene_cells, pose );
            for ( const SurfelInfoConstPtr & cell : scene_cells )
            {
                if ( cell == nullptr || !cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                //if ( cell->m_surfel->getNumPoints() < 10.0 ) p.addPoint( cell->m_center_s, rc );
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);

                const Eigen::Vector3f h = -cell->m_surfel->normal_.normalized();
                a += h * h.transpose();
                an += h * cell->m_surfel->getNumPoints() * h.transpose();
                ++n;
                wn+=cell->m_surfel->getNumPoints();
            }
            if ( !a.isZero() )
            {
                const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigs(a / n);
                const Eigen::Matrix3f sevs = eigs.eigenvectors();
                const Eigen::Vector3f sev = eigs.eigenvalues().normalized() * 10;

                VisMesh p;
                Eigen::Vector3f v1x = (sevs.col(0) * sev(0));
                if ( v1x.dot(Eigen::Vector3f::UnitX()) < 0 ) v1x *= -1;
                Eigen::Vector3f v1y = (sevs.col(1) * sev(1));
                if ( v1y.dot(Eigen::Vector3f::UnitY()) < 0 ) v1y *= -1;
                Eigen::Vector3f v1z = (sevs.col(2) * sev(2));
                if ( v1z.dot(Eigen::Vector3f::UnitZ()) < 0 ) v1z *= -1;

                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1x, rc );
                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1y, rc);
                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1z, rc );
                p.showEdgeMesh("Eigs_"+std::to_string(idx),true);
                LOG(1) << "Kappa: " << (std::abs(eigs.eigenvalues()(2))/(std::abs(eigs.eigenvalues()(0))));
                LOG(1) << "Att=[" << a.row(0) <<"; "<< a.row(1) << "; " << a.row(2) << "]";
                LOG(1) << "Evs=[" << sevs.row(0) <<"; "<< sevs.row(1) << "; " << sevs.row(2) << "]";
            }
            if ( !an.isZero() )
            {
                const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigs(an / wn);
                const Eigen::Matrix3f sevs = eigs.eigenvectors();
                const Eigen::Vector3f sev = eigs.eigenvalues().normalized() * 10;

                VisMesh p;
                Eigen::Vector3f v1x = (sevs.col(0) * sev(0));
                if ( v1x.dot(Eigen::Vector3f::UnitX()) < 0 ) v1x *= -1;
                Eigen::Vector3f v1y = (sevs.col(1) * sev(1));
                if ( v1y.dot(Eigen::Vector3f::UnitY()) < 0 ) v1y *= -1;
                Eigen::Vector3f v1z = (sevs.col(2) * sev(2));
                if ( v1z.dot(Eigen::Vector3f::UnitZ()) < 0 ) v1z *= -1;

                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1x, rc );
                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1y, rc);
                p.addEdge(interp_scene_pose.translation().cast<float>(), interp_scene_pose.cast<float>()*v1z, rc );
                p.showEdgeMesh("NEigs_"+std::to_string(idx),true);
                LOG(1) << "NKappa: " << (std::abs(eigs.eigenvalues()(2))/(std::abs(eigs.eigenvalues()(0))));
                LOG(1) << "NAtt=[" << an.row(0) <<"; "<< an.row(1) << "; " << an.row(2) << "]";
                LOG(1) << "NEvs=[" << sevs.row(0) <<"; "<< sevs.row(1) << "; " << sevs.row(2) << "]";
            }
        }
    }
#endif

#ifdef USE_EASY_PBR
        const bool showSceneSurfels = m_show_scene_surfels;
        if ( showSceneSurfels )
        {

            StopWatch total, s1,c1;
            //std::vector<MarsMap *> sceneMapPtrVec = m_sceneMarsMap->getTransformedMapPtrVec(prev_scene_pose);
            Eigen::VectorXt newCurTimes = sceneTimes;
            for ( size_t idx = 0; idx < sceneMapPtrVec.size(); ++idx)
            {
                MarsMap * scene = sceneMapPtrVec[idx];
                const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( newCurTimes(idx) ); // curTimes still contains second
                const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                Eigen::Matrix3f a = Eigen::Matrix3f::Zero();
                Eigen::Vector6f h = Eigen::Vector6f::Zero();
                Eigen::Matrix6f A = Eigen::Matrix6f::Zero();
                VisMesh v, n, f, e, p;
                Eigen::Vector3f rc = VisMesh::getPlasma(float(idx+1) / sceneMapPtrVec.size()); //VisMesh::randColor();
                Sophus::SE3d pose;
                SurfelInfoConstPtrVector scene_cells;
                scene->getCells( scene_cells, pose );
                c1.reset();
                v.reserveSurfelTriangleEllipsoid(scene_cells.size());
                for ( const SurfelInfoConstPtr & cell : scene_cells )
                {
                    if ( cell == nullptr || !cell->m_surfel ) continue;
                    if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                    //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                    //if ( cell->m_surfel->getNumPoints() < 10.0 ) p.addPoint( cell->m_center_s, rc );
                    if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                    //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                    //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);

                    const Eigen::Vector3f m = cell->m_center_s+cell->m_surfel->mean_;
                    h.head<3>() = -m.cross(cell->m_surfel->normal_);
                    h.tail<3>() = -cell->m_surfel->normal_;
                    A += h * h.transpose();
                    a += h.tail<3>() * h.tail<3>().transpose();

                    v.addSurfel( cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, (cell->m_center_s + cell->m_surfel->mean_), rc,scene->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                }

                const Eigen::Matrix3f Att = A.block<3,3>(3,3);
                if ( !Att.isZero() )
                {
                    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigs(Att);
                    const Eigen::Matrix3f sevs = eigs.eigenvectors();
                    const Eigen::Vector3f sev = eigs.eigenvalues().normalized() * 10;
                    VisMesh p;
                    for ( int c = 0; c < sevs.cols(); ++c )
                    {
                        const Eigen::Vector3d v1 = interp_scene_pose * Eigen::Vector3d::Zero();
                        const Eigen::Vector3d v2 = interp_scene_pose * (sevs.col(c) * sev(c)).cast<double>();
                        p.addEdge(v1.cast<float>(), v2.cast<float>(), rc);
                    }
                    p.showEdgeMesh("eigs_"+std::to_string(idx),true);
                    LOG(1) << "kappa: " << (std::abs(eigs.eigenvalues()(2))/(std::abs(eigs.eigenvalues()(0))));
                    LOG(1) << "Att=[" << Att.row(0) <<"; "<< Att.row(1) << "; " << Att.row(2) << "]";
                    LOG(1) << "att=[" << a.row(0) <<"; "<< a.row(1) << "; " << a.row(2) << "]";
                    LOG(1) << "evs=[" << sevs.row(0) <<"; "<< sevs.row(1) << "; " << sevs.row(2) << "]";
                }

                LOG(1) << "\n\nc1: " << c1.getTime();
                s1.reset();
                v.showSurfelMesh("sceneSurfelsPreInt_"+std::to_string(idx), true, interp_scene_pose,m_surfel_type);
                LOG(1) << "s1: " << s1.getTime();
                //n.showEdgeMesh("sceneSurfelPreIntNormals_"+std::to_string(idx), true, interp_scene_pose);
                //f.showEdgeMesh("sceneSurfelPreIntFirstViewDir_"+std::to_string(idx), true, interp_scene_pose);
                //e.showPointMesh("sceneSurfelPreIntNotValidMean_"+std::to_string(idx), true, interp_scene_pose);
                //p.showPointMesh("sceneSurfelPreIntTooLittlePoints_"+std::to_string(idx), true, interp_scene_pose);
            }
            LOG(1) << "total: " << total.getTime();
        }
#endif
#ifdef USE_EASY_PBR
    const bool showOptSceneMap = m_show_opt_scene;
    if ( showOptSceneMap)
    {
        VisMesh::Ptr vm = VisMesh::create();
        std::vector<MarsMapPointCloud::Ptr> scene_clouds = m_sceneMarsMap->getCloudSharedPtrVec();
        for ( size_t cloud_idx = 0; cloud_idx < scene_clouds.size(); ++cloud_idx)
        {
            VisMesh::Ptr cm = VisMesh::create();
            const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( curTimes(cloud_idx) ); // curTimes still contains second
            const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
            MarsMapPointCloud::Ptr scene_cloud = scene_clouds[cloud_idx];
            const Eigen::Matrix3Xf & pts = scene_cloud->toMatrix3Xf();
            const Eigen::VectorXf & intensity = scene_cloud->intensity();
            const Eigen::VectorXu & reflectivity = scene_cloud->reflectivity();
            const Eigen::VectorXu & scan_line = scene_cloud->m_scan_line_id;
            const int num_pts = scene_cloud->size();
            for ( int idx = 0; idx < num_pts; ++idx )
            {
                const Eigen::Vector3f pt = (interp_scene_pose * pts.col(idx).cast<double>()).template cast<float>();
                vm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                cm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
            }
            cm->showPointMesh("OptSceneCloud"+std::to_string(m_scan_id)+"_"+std::to_string(cloud_idx), true);
        }
        vm->showPointMesh("OptSceneMap"+std::to_string(m_scan_id), true);
    }
#endif

    if ( m_sceneMarsMap->hasFullWindow() )
    {
        //LOG(1) << "Moving Scene!";
        {
            std::unique_lock<std::mutex> lock(m_sceneMapMutex);
            m_sceneMarsMap->moveWindowToNext();
        }
        //LOG(1) << "Moved Scene!";
        // get last cloud from scene map
        MarsMapPointCloud::Ptr old_last_scene_cloud = m_sceneMarsMap->getOutOfWindowCloud();
        assert ( old_last_scene_cloud != nullptr );

        Eigen::VectorXt newCurTimes = m_sceneMarsMap->getTimesSince( last_time );
        Sophus::SE3d interp_pose = m_trajectory_spline->pose( curTimes(0) ); // old time

        Sophus::SE3d new_scene_pose = prev_scene_pose * interp_pose;
        const Sophus::SE3d new_scene_map_pose = m_localMarsMap->local_map_pose * prev_scene_pose * interp_pose;
        m_old_poses.emplace_back(std::pair<int64_t, Sophus::SE3d>(old_last_scene_cloud->m_stamp,new_scene_map_pose));
        while ( int(m_old_poses.size()) > m_init_spline_window )
            m_old_poses.pop_front();

#ifdef USE_EASY_PBR
        const bool showPrevMarsSurfels = false;
        if ( showPrevMarsSurfels )
        {
            VisMesh v, n, f, e, p;
            Eigen::Vector3f rc = VisMesh::randColor();
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            asOneCenteredLocalMap->getCells( scene_cells, pose );
            v.reserveSurfelTriangleEllipsoid(scene_cells.size());
            for ( const SurfelInfoConstPtr & cell : scene_cells )
            {
                if ( cell == nullptr || !cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                //if ( ! cell->m_surfel->valid_ && cell->m_surfel->getNumPoints() >= 20.0 )
                //{
                //    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eig(cell->m_surfel->cov_);
                //    Eigen::Matrix3f eigen_vectors = eig.eigenvectors();
                //    Eigen::Vector3f eigen_values = eig.eigenvalues();
                //    e.addSurfelTriangleEllipsoid(eigen_vectors, eigen_values, cell->m_center_s + cell->m_surfel->mean_, rc, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2);
                //    //e.addPoint( cell->m_center_s, rc );
                //}
                //if ( ! cell->m_surfel->valid_ && cell->m_surfel->getNumPoints() >= 10.0 ) LOG(1) << "PrevMarsSurfel: invalid. mcs: " << cell->m_center_s.transpose() << " #:" << cell->m_surfel->getNumPoints() << " cov=[" << cell->m_surfel->cov_.row(0)<<"; " << cell->m_surfel->cov_.row(1) << "; " <<cell->m_surfel->cov_.row(2)<<"]; dc: " << cell->m_surfel->cov_.determinant() << " ddc: " << cell->m_surfel->cov_.determinant()<< " eps: " << std::numeric_limits<Surfel::Scalar>::epsilon();
                //if ( cell->m_surfel->getNumPoints() < 10.0 ) p.addPoint( cell->m_center_s, rc );
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);
                v.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, rc, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
            }
            v.showSurfelMesh("marsPrevSurfels", true, m_localMarsMap->local_map_pose,m_surfel_type);
            //n.showEdgeMesh("marsPrevSurfelsNormals", true, asOneCenteredLocalMap->getMapPose());
            //f.showEdgeMesh("marsPrevSurfelsFirstViewDir", true, asOneCenteredLocalMap->getMapPose());
            //e.showSurfelTriangleEllipsoidMesh("marsPrevSurfelNotValidMean", true, asOneCenteredLocalMap->getMapPose());
            //p.showPointMesh("marsPrevSurfelTooLittlePoints", true, asOneCenteredLocalMap->getMapPose());
        }
#endif

        static Sophus::SE3d last_kf_pose = new_scene_pose;
        if ( m_use_closest_kf ) last_kf_pose = m_localMarsMap->getClosestKeyFramePose( new_scene_pose );

        const Sophus::SO3d diff_rot = (last_kf_pose.so3().inverse() * new_scene_pose.so3());
        const bool kf_from_rotation = Eigen::Vector2d(diff_rot.angleX(),diff_rot.angleY()).norm() > m_min_keyframe_rotation;
        const bool shouldCreateKeyFrame = (last_kf_pose.translation()-new_scene_pose.translation()).norm() > m_min_keyframe_distance || (m_keyframe_use_rotation && kf_from_rotation);
        //LOG(1) << "KF? " << shouldCreateKeyFrame << " t: " << (last_kf_pose.translation()-new_scene_pose.translation()).transpose() << " nt: " << (last_kf_pose.translation()-new_scene_pose.translation()).squaredNorm() << " th: " << m_min_keyframe_distance;
        if ( shouldCreateKeyFrame )
        {
            m_was_keyframe = true;
            m_last_key_frame_stamp = old_last_scene_cloud->m_stamp;
            LOG(1) << "last_kf_t: " << last_kf_pose.translation().transpose() << " new_scene_t: "<<  new_scene_pose.translation().transpose();
            //LOG(1) << "last_kf_t: " << (m_localMarsMap->local_map_pose*prev_scene_pose*last_kf_pose.translation()).transpose() << " new_scene_t: "<<  (m_localMarsMap->local_map_pose*prev_scene_pose*new_scene_pose.translation()).transpose() << " diff: " << (m_localMarsMap->local_map_pose*prev_scene_pose*last_kf_pose.translation()-m_localMarsMap->local_map_pose*prev_scene_pose*new_scene_pose.translation()).transpose() << " e: " << (m_localMarsMap->local_map_pose*prev_scene_pose*last_kf_pose.translation()-m_localMarsMap->local_map_pose*prev_scene_pose*new_scene_pose.translation()).norm();

            Eigen::Matrix6d covar = Eigen::Matrix6d::Zero();
            m_marsRegistrator->getCovariance( 0, new_scene_pose, covar );

#ifdef USE_EASY_PBR
            {
                const int num_points = old_last_scene_cloud->size();
                VisMesh v;
                v.reservePoints(num_points);
                for ( int idx = 0; idx < num_points; ++idx )
                    v.addPoint(old_last_scene_cloud->m_points.col(idx), VisMesh::red);

                v.showPointMesh("last_kf",true,m_localMarsMap->local_map_pose * new_scene_pose);
            }
#endif

            {
            ZoneScopedN("MarsSplineRegistrator::addCloudToMap");
            std::unique_lock<std::mutex> lock(m_localMapMutex);
            m_localMarsMap->addCloud( old_last_scene_cloud, new_scene_pose ); // "new_scene_pose" stored in the map might be outdated from now on.
            }
            if ( m_localMarsMap->local_map_moved )
            {
                asOneCenteredLocalMap = m_localMarsMap->getAsOneCenteredMap();
                LOG(1) << "local map moved... pose: " << m_localMarsMap->local_map_last_moving_pose.params().transpose() << "\nnew_scene: " << new_scene_pose.params().transpose() << "\ninterp_pose: " << interp_pose.params().transpose();

                const Sophus::SE3d move_pose = m_localMarsMap->local_map_last_moving_pose.inverse();
                m_interpolated_scene_pose = move_pose * m_interpolated_scene_pose;
                new_scene_pose = move_pose * new_scene_pose;
                interp_pose = move_pose * interp_pose;
                for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
                {
                    m_trajectory_spline->setKnot(move_pose*m_trajectory_spline->getKnot(knot_idx),knot_idx);
                }
                for ( auto & scene_pose : m_scene_poses )
                    scene_pose = move_pose * scene_pose;
                newest_cloud_pose = m_scene_poses.back();

#ifdef USE_EASY_PBR
                const bool showMovedKnots = m_show_knots;
                if ( showMovedKnots )
                {
                    VisMesh::Ptr tm = VisMesh::create();
                    VisMesh::Ptr bm = VisMesh::create();
                    VisMesh::Ptr km = VisMesh::create();
                    for ( size_t knot_idx = 0; knot_idx < m_trajectory_spline->numKnots(); ++knot_idx )
                    {
                        const Sophus::SE3d knot = m_trajectory_spline->getKnot(knot_idx);
                        km->addPoint( (m_localMarsMap->local_map_pose * prev_scene_pose * knot.translation()).cast<float>(), VisMesh::blue);
                    }
                    for ( const Sophus::SE3d & scan_pose : m_scene_poses )
                    {
                        bm->addPoint((m_localMarsMap->local_map_pose * scan_pose.translation()).cast<float>(), VisMesh::green);
                    }
                    constexpr uint64_t num_samples = 1000;
                    const uint64_t dt = m_trajectory_spline->getDtNs();
                    uint64_t dt_inc = dt / num_samples;
                    Sophus::SE3d prev_interp_scene_pose;
                    bool prev_is_set = false;
                    for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
                    {
                        const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( m_trajectory_spline->minTimeNs() + dt_inc * t_inc ); // curTimes still contains second
                        const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                        if ( prev_is_set )
                            tm->addEdge( prev_interp_scene_pose.translation().cast<float>(), interp_scene_pose.translation().cast<float>(), VisMesh::getPlasma( float(t_inc) / num_samples ) );
                        prev_interp_scene_pose = interp_scene_pose;
                        prev_is_set = true;
                    }
                    km->showPointMesh("MovedKnotPoses",true);
                    bm->showPointMesh("MovedScanPose",true);
                    tm->showEdgeMesh("MovedSpline",true);
                }
#endif
            }
            last_kf_pose = new_scene_pose;
        }
#ifdef USE_EASY_PBR
        const bool showMarsMap = m_show_mars_map; // && m_scan_id == 2400;
        if ( showMarsMap )
        {
            Eigen::Vector3f rc = VisMesh::randColor();
            static VisMesh::Ptr vm = nullptr;
            static VisMesh::Ptr am = nullptr;
            static VisMesh::Ptr tm = nullptr;
            static VisMesh::Ptr om = nullptr;
            std::vector<MarsMapPointCloud::Ptr> clouds = m_localMarsMap->getCloudSharedPtrVec();
            std::vector<Sophus::SE3d> cloud_poses = m_localMarsMap->getCloudPosesVec();
            if ( ! am  )
            {
                vm = VisMesh::create();
                am = VisMesh::create();
                tm = VisMesh::create();
                om = VisMesh::create();

                //if ( clouds.size() != 2 ) LOG(FATAL) << "why does it contain more than two clouds?";
                for ( size_t cloudIdx = 0; cloudIdx < clouds.size(); ++cloudIdx)
                {
                    MarsMapPointCloud::Ptr cloud = clouds[cloudIdx];
                    Sophus::SE3d relTf = m_localMarsMap->local_map_pose * cloud_poses[cloudIdx];


                    const Eigen::Matrix3Xf & pts = cloud->toMatrix3Xf();
                    const Eigen::VectorXf & intensity = cloud->intensity();
                    const Eigen::VectorXu & reflectivity = cloud->reflectivity();
                    const Eigen::VectorXu & scan_line = cloud->m_scan_line_id;
                    const int num_pts = cloud->size();
                    for ( int idx = 0; idx < num_pts; ++idx )
                    {
                        const Eigen::Vector3f pt = (relTf * pts.col(idx).cast<double>()).template cast<float>();
                        //am->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                        vm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                    }
                }
                tm->addPoint ( m_localMarsMap->local_map_pose.translation().cast<float>(), rc );
                //am->showAggregatedPointMesh(Sophus::SE3d(),"allMap");
                vm->showAggregatedPointMesh(Sophus::SE3d(),"marsMap", 0, m_scan_id >= 2350);
                tm->showAggregatedPointMesh(Sophus::SE3d(),"marsTraj");

                om->addEdge ( m_localMarsMap->local_map_pose.translation().cast<float>(), (m_localMarsMap->local_map_pose * (Eigen::Vector3d::UnitX()*0.1)).cast<float>(), VisMesh::red );
                om->addEdge ( m_localMarsMap->local_map_pose.translation().cast<float>(), (m_localMarsMap->local_map_pose * (Eigen::Vector3d::UnitY()*0.1)).cast<float>(), VisMesh::blue );
                om->addEdge ( m_localMarsMap->local_map_pose.translation().cast<float>(), (m_localMarsMap->local_map_pose * (Eigen::Vector3d::UnitZ()*0.1)).cast<float>(), VisMesh::green );
                om->showAggregatedEdgeMesh(Sophus::SE3d(),"marsOriTraj");

                rc = VisMesh::randColor();
            }
            else
            {
                if (shouldCreateKeyFrame )
                {
                    constexpr bool showEachKeyFrameSolo = false;
                    if ( !showEachKeyFrameSolo )
                        vm->resetForAggregation();
                    else
                        vm = VisMesh::create();
                    MarsMapPointCloud::Ptr cloud = old_last_scene_cloud;
                    Sophus::SE3d relTf = m_localMarsMap->local_map_pose * new_scene_pose;

                    const Eigen::Matrix3Xf & pts = cloud->toMatrix3Xf();
                    const Eigen::VectorXf  & intensity = cloud->intensity();
                    const Eigen::VectorXu& reflectivity = cloud->reflectivity();
                    const Eigen::VectorXu & scan_line = cloud->m_scan_line_id;
                    const int num_pts = cloud->size();
                    for ( int idx = 0; idx < num_pts; ++idx )
                    {
                        const Eigen::Vector3f pt = (relTf * pts.col(idx).cast<double>()).template cast<float>();
                        vm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                    }
                    if ( !showEachKeyFrameSolo )
                        vm->showAggregatedPointMesh(Sophus::SE3d(),"marsMap", 0, true );//m_scan_id >= 2350);
                    else
                        vm->showAggregatedPointMesh(Sophus::SE3d(),"marsMap_"+std::to_string(m_scan_id), 0, true );//m_scan_id >= 2350);
                }
            }

            //am->resetForAggregation();
            //for ( size_t i = 0; i < old_last_scene_cloud->points.size(); ++i )
            //{
            //    am->addPoint( m_localMarsMap->local_map_pose * new_scene_pose * old_last_scene_cloud->points[i].pt, Eigen::Vector3f::Constant(old_last_scene_cloud->points[i].scan_line_id/128.), &old_last_scene_cloud->points[i].intensity,&old_last_scene_cloud->points[i].reflectivity );
            //}
            tm->resetForAggregation();
            om->resetForAggregation();
            if ( curTimes.rows() > 1 )
            {
                constexpr uint64_t num_samples = 10;
                const uint64_t dt = curTimes(1) - curTimes(0);
                uint64_t dt_inc = dt / num_samples;
                for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
                {
                    const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( curTimes(0) + dt_inc * t_inc ); // curTimes still contains second
                    const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                    tm->addPoint( interp_scene_pose.translation().cast<float>(), rc );
                    if ( t_inc == 4 )
                    {
                        om->addEdge ( interp_scene_pose.translation().cast<float>(), (interp_scene_pose * (Eigen::Vector3d::UnitX()*0.1)).cast<float>(), VisMesh::red );
                        om->addEdge ( interp_scene_pose.translation().cast<float>(), (interp_scene_pose * (Eigen::Vector3d::UnitY()*0.1)).cast<float>(), VisMesh::blue );
                        om->addEdge ( interp_scene_pose.translation().cast<float>(), (interp_scene_pose * (Eigen::Vector3d::UnitZ()*0.1)).cast<float>(), VisMesh::green );
                    }
                }
            }
            tm->addPoint( (m_localMarsMap->local_map_pose*new_scene_pose).translation().cast<float>(), rc );
            tm->showAggregatedPointMesh(Sophus::SE3d(),"marsTraj");
            om->showAggregatedEdgeMesh(Sophus::SE3d(),"marsOriTraj");
            //am->showAggregatedPointMesh(Sophus::SE3d(),"allMap" ); //, 0.8, m_scan_id > 980 );
        }
#endif

#ifdef USE_EASY_PBR
        const bool showMarsSurfels = m_show_mars_surfels;
        if ( showMarsSurfels ) //&& shouldCreateKeyFrame )
        {
            constexpr bool showSemanticColored = has_semantics_v<MarsMapPointCloud>;
            VisMesh v, nc, n, f, e, p,s;
            Eigen::Vector3f rc = VisMesh::getPlasma(0);// VisMesh::randColor();
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            asOneCenteredLocalMap->getCells( scene_cells, pose );
            v.reserveSurfelTriangleEllipsoid(scene_cells.size());
            nc.reserveSurfelTriangleEllipsoid(scene_cells.size());
            if constexpr ( showSemanticColored )
                s.reserveSurfelTriangleEllipsoid(scene_cells.size());
            LOG(1) << "ShowMarsSurfels: " << scene_cells.size() << " pts: " << asOneCenteredLocalMap->m_num_points << " p:" << pose.params().transpose();
            for ( const SurfelInfoConstPtr & cell : scene_cells )
            {
                if ( cell == nullptr || !cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "Not evaluated ? this should not have happend!";
                //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                //if ( cell->m_surfel->getNumPoints() < 10.0 ) p.addPoint( cell->m_center_s, rc );
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);
                //if ( (cell->m_surfel->eigen_values_.array().sqrt().maxCoeff()*3) > 5 )
                //{
                //    LOG(1) << "Upsi Daisy: ct: "<< cell->m_center_s.transpose() << " m: " << cell->m_surfel->mean_.transpose()<< " cov=[" << cell->m_surfel->cov_.row(0) << "; "<< cell->m_surfel->cov_.row(1) << "; " << cell->m_surfel->cov_.row(2)<<"] ";
                //}
                if constexpr ( ! semantics_only )
                {
                    v.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, rc, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                    nc.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, VisMesh::colorFromNormal(cell->m_surfel->normal_.normalized()), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                }

                if constexpr ( showSemanticColored )
                    if ( cell->m_class )
                    s.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, m_label_manager->color_for_label(cell->m_class->getArgMaxClass()).cast<float>(), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
            }
            v.showSurfelMesh("marsSurfels", true, m_localMarsMap->local_map_pose,m_surfel_type);
            nc.showSurfelMesh("marsSurfelsN", true, m_localMarsMap->local_map_pose,m_surfel_type);
            if constexpr ( showSemanticColored )
                s.showSurfelMesh("marsSurfelsSem", true, m_localMarsMap->local_map_pose,m_surfel_type);

            for ( int lvl = 0; lvl < asOneCenteredLocalMap->m_map_params.m_num_levels; ++lvl)
            {
                VisMesh vl, vn, vs;
                SurfelInfoVector cells;
                asOneCenteredLocalMap->getCellsOnLevel ( lvl, cells, false);
                vl.reserveSurfelTriangleEllipsoid(cells.size());
                vn.reserveSurfelTriangleEllipsoid(cells.size());
                if constexpr ( showSemanticColored )
                    vs.reserveSurfelTriangleEllipsoid(cells.size());
                for ( const SurfelInfo & cell : cells )
                {
                    if ( ! cell.m_surfel ) continue;
                    if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 ) continue;
                    if constexpr ( ! semantics_only )
                    {
                        vl.addSurfel(cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, cell.m_center_s + cell.m_surfel->mean_, rc, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell.m_level)/2, m_surfel_type);
                        vn.addSurfel(cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, cell.m_center_s + cell.m_surfel->mean_, VisMesh::colorFromNormal(cell.m_surfel->normal_.normalized()), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell.m_level)/2, m_surfel_type);
                    }
                    if constexpr ( showSemanticColored )
                        if ( cell.m_class )
                        vs.addSurfel(cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, cell.m_center_s + cell.m_surfel->mean_, m_label_manager->color_for_label(cell.m_class->getArgMaxClass()).cast<float>(), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell.m_level)/2, m_surfel_type);
                }
                if constexpr ( ! semantics_only )
                vl.showSurfelMesh("marsSurfels"+std::to_string(lvl), true, m_localMarsMap->local_map_pose, m_surfel_type);
                if constexpr ( ! semantics_only )
                vn.showSurfelMesh("marsSurfelsN"+std::to_string(lvl), true, m_localMarsMap->local_map_pose, m_surfel_type);
                if constexpr ( showSemanticColored )
                    vs.showSurfelMesh("marsSurfelsSem"+std::to_string(lvl), true, m_localMarsMap->local_map_pose, m_surfel_type);
            }

            //n.showEdgeMesh("marsSurfelsNormals", true, asOneCenteredLocalMap->getMapPose());
            //f.showEdgeMesh("marsSurfelsFirstViewDir", true, asOneCenteredLocalMap->getMapPose());
            //e.showPointMesh("marsSurfelNotValidMean", true, asOneCenteredLocalMap->getMapPose());
            //p.showPointMesh("marsSurfelTooLittlePoints", true, asOneCenteredLocalMap->getMapPose());
        }
#endif
#ifdef USE_EASY_PBR
        const bool showSceneMap = m_show_scene_map; // && m_scan_id % 100 == 0;
        if ( showSceneMap )
        {
            constexpr bool showCloudSemanticColored = has_semantics_v<MarsMapPointCloud>;
            constexpr bool showCloudHeightColored = false;
            constexpr bool showSceneViewDir = false;
            VisMesh::Ptr vm = VisMesh::create(), tm = VisMesh::create(), em = VisMesh::create();
            std::vector<MarsMapPointCloud::Ptr> scene_clouds = m_sceneMarsMap->getCloudSharedPtrVec();
            for ( size_t cloud_idx = 0; cloud_idx < scene_clouds.size(); ++cloud_idx)
            {
                //VisMesh::Ptr cm = VisMesh::create();
                const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( newCurTimes(cloud_idx) ); // curTimes still contains second
                const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;

                MarsMapPointCloud::Ptr scene_cloud = scene_clouds[cloud_idx];
                const Eigen::Vector3f rc = VisMesh::randColor();
                vm = VisMesh::create();
                const Eigen::Matrix3Xf & pts = scene_cloud->toMatrix3Xf();
                const Eigen::VectorXf  & intensity = scene_cloud->intensity();
                const Eigen::VectorXu& reflectivity = scene_cloud->reflectivity();
                //const Eigen::VectorXu & scan_line = scene_cloud->m_scan_line_id;
                const int num_pts = scene_cloud->size();
                for ( int idx = 0; idx < num_pts; ++idx )
                {
                    const Eigen::Vector3f pt = (interp_scene_pose * pts.col(idx).cast<double>()).template cast<float>();
                    vm->addPoint( pt,  VisMesh::getPlasma(float(cloud_idx+1)/(scene_clouds.size()+1)), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                    if constexpr (showSceneViewDir )
                    {
                        const Eigen::Vector3f view_dir = (pt - interp_scene_pose.translation().cast<float>()).normalized();
                        em->addEdge( pt, pt - view_dir * 0.1, VisMesh::colorFromNormal( view_dir ) );
                    }
                    //cm->addPoint( pt, scan_line.rows()>0 ? Eigen::Vector3f::Constant(scan_line(idx)/128) : Eigen::Vector3f::Zero(), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
                }

                if constexpr ( showCloudHeightColored )
                {
                    CloudPtr vmm = vm->getPointMesh();
                    vmm->worldROS2worldGL();
                    vmm->m_vis.m_color_type = easy_pbr::MeshColorType::Height;
                    vmm->m_vis.m_point_size = 5;
                    vmm->m_min_max_y_for_plotting = Eigen::Vector2f(-3,1);
                    easy_pbr::Scene::show(vmm, "sceneCloud_"+std::to_string(cloud_idx+1));
                }
                else
                    vm->showPointMesh("sceneCloud_"+std::to_string(cloud_idx+1), true);


                if constexpr ( showCloudSemanticColored )
                {
                    VisMesh::Ptr sm = showSemanticColored<MarsMapPointCloud>( scene_cloud, interp_scene_pose );
                    CloudPtr vmm = sm->getPointMesh();
                    vmm->worldROS2worldGL();

                    if ( vmm->C.rows() > 0 )
                    {
                        vmm->L_pred = vmm->C.col(0).cast<int>();
                        vmm->L_gt = vmm->L_pred;
                    }

                    //vmm->m_vis.m_color_type = easy_pbr::MeshColorType::Height;
                    vmm->m_vis.m_point_size = 5;
                    vmm->m_label_mngr=m_label_manager->shared_from_this();
                    vmm->m_vis.set_color_semanticpred();
                    easy_pbr::Scene::show(vmm, "sceneSemantics"+std::to_string(cloud_idx+1));

                    //sm->showPointMesh("sceneSemantics"+std::to_string(cloud_idx+1),true);
                }

                LOG(1) << "SceneCloud["<<cloud_idx<<"] has " << scene_cloud->size();
                if ( int(cloud_idx)+1 <= curTimes.rows() )
                {
                    constexpr uint64_t num_samples = 10;
                    const uint64_t dt = curTimes(cloud_idx+1) - curTimes(cloud_idx);
                    uint64_t dt_inc = dt / num_samples;
                    for ( size_t t_inc = 0; t_inc < num_samples; ++t_inc )
                    {
                        const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( curTimes(cloud_idx) + dt_inc * t_inc ); // curTimes still contains second
                        const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                        tm->addPoint( interp_scene_pose.translation().cast<float>(), rc );
                    }
                }
                tm->addPoint( interp_scene_pose.translation().cast<float>(), rc );
                //cm->showPointMesh("sceneCloud"+std::to_string(cloud_idx), true);
            }
            //vm->showPointMesh("sceneMap"+std::to_string(m_scan_id), true);

            if constexpr ( showSceneViewDir ) em->showEdgeMesh("sceneViewDir", true);
            tm->showEdgeMeshFromPoints("sceneTraj", true);

            VisMesh::Ptr lm = VisMesh::create();

            const Eigen::Matrix3Xf & pts = old_last_scene_cloud->toMatrix3Xf();
            const Eigen::VectorXf & intensity = old_last_scene_cloud->intensity();
            const Eigen::VectorXu& reflectivity = old_last_scene_cloud->reflectivity();
            //const Eigen::VectorXu & scan_line = old_last_scene_cloud->m_scan_line_id;
            const int num_pts = old_last_scene_cloud->size();
            for ( int idx = 0; idx < num_pts; ++idx )
            {
                lm->addPoint( (new_scene_map_pose*pts.col(idx).cast<double>()).template cast<float>(), VisMesh::getPlasma(0), intensity.rows()>0 ? &intensity(idx) : nullptr, reflectivity.rows()>0 ? &reflectivity(idx) : nullptr);
            }
            if constexpr ( showCloudHeightColored )
            {
                CloudPtr vmm = lm->getPointMesh();
                vmm->worldROS2worldGL();
                vmm->m_vis.m_point_size = 5;
                vmm->m_vis.m_color_type = easy_pbr::MeshColorType::Height;
                vmm->m_min_max_y_for_plotting = Eigen::Vector2f(-3,1);
                easy_pbr::Scene::show(vmm, "sceneCloud_0");
            }
            else
                lm->showPointMesh("sceneCloud_0", true);
            if constexpr ( showCloudSemanticColored )
            {
                VisMesh::Ptr sm = showSemanticColored<MarsMapPointCloud>( old_last_scene_cloud, new_scene_map_pose );
                CloudPtr vmm = sm->getPointMesh();
                vmm->worldROS2worldGL();

                if ( vmm->C.rows() > 0 )
                {
                    vmm->L_pred = vmm->C.col(0).cast<int>();
                    vmm->L_gt = vmm->L_pred;
                }
                //vmm->m_vis.m_color_type = easy_pbr::MeshColorType::Height;
                vmm->m_vis.m_point_size = 5;
                vmm->m_label_mngr=m_label_manager->shared_from_this();
                vmm->m_vis.set_color_semanticpred();
                easy_pbr::Scene::show(vmm, "sceneSemantics0");
            }
        }
#endif

#ifdef USE_EASY_PBR
        const bool showSceneSurfels = m_show_scene_surfels;
        if ( showSceneSurfels )
        {
            StopWatch total;
            constexpr bool showSemanticColored = has_semantics_v<MarsMapPointCloud>;
            MarsMap * last_scene = sceneMapPtrVec[0];
            VisMesh v, n, s;
            const Eigen::Vector3f rc = VisMesh::randColor();
            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            last_scene->getCells( scene_cells, pose );
            v.reserveSurfelTriangleEllipsoid(scene_cells.size());
            n.reserveSurfelTriangleEllipsoid(scene_cells.size());
            if constexpr ( showSemanticColored )
                s.reserveSurfelTriangleEllipsoid(scene_cells.size());
            for ( const SurfelInfoConstPtr & cell : scene_cells )
            {
                if ( cell == nullptr || !cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);
                if ( ! semantics_only )
                {
                    v.addSurfel( cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, (cell->m_center_s + cell->m_surfel->mean_), rc,last_scene->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                    n.addSurfel( cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, (cell->m_center_s + cell->m_surfel->mean_), VisMesh::colorFromNormal(cell->m_surfel->normal_),last_scene->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                }
                if constexpr ( showSemanticColored )
                    if ( cell->m_class )
                    s.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, m_label_manager->color_for_label(cell->m_class->getArgMaxClass()).cast<float>(), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
            }
            v.showSurfelMesh("lastSceneSurfels", true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
            n.showSurfelMesh("lastSceneSurfelsN", true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
            if constexpr ( showSemanticColored )
                s.showSurfelMesh("lastSceneSurfelsSem", true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
            //n.showEdgeMesh("lastSceneSurfelNormals", true, new_scene_pose);
            //f.showEdgeMesh("lastSceneSurfelFirstViewDir", true, new_scene_pose);

            std::vector<MarsMap *> sceneMapPtrVec = m_sceneMarsMap->getTransformedMapPtrVec(prev_scene_pose);
            for ( size_t idx = 0; idx < sceneMapPtrVec.size(); ++idx)
            {
                MarsMap * scene = sceneMapPtrVec[idx];
                const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( newCurTimes(idx) ); // curTimes still contains second
                const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                VisMesh v, n, f, e, p, s;
                const Eigen::Vector3f rc = VisMesh::randColor();
                Sophus::SE3d pose;
                SurfelInfoConstPtrVector scene_cells;
                scene->getCells( scene_cells, pose );
                v.reserveSurfelTriangleEllipsoid(scene_cells.size());
                n.reserveSurfelTriangleEllipsoid(scene_cells.size());
                if constexpr ( showSemanticColored )
                    s.reserveSurfelTriangleEllipsoid(scene_cells.size());
                for ( const SurfelInfoConstPtr & cell : scene_cells )
                {
                    if ( cell == nullptr || !cell->m_surfel ) continue;
                    if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                    //if ( ! cell->m_surfel->valid_ ) e.addPoint( cell->m_center_s, rc );
                    //if ( cell->m_surfel->getNumPoints() < 10.0 ) p.addPoint( cell->m_center_s, rc );
                    if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                    //n.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->normal_, rc);
                    //f.addNormal( (cell->m_center_s + cell->m_surfel->mean_), cell->m_surfel->first_view_dir_, rc);
                    if ( ! semantics_only )
                    {
                        v.addSurfel( cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, (cell->m_center_s + cell->m_surfel->mean_), rc, scene->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                        n.addSurfel( cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, (cell->m_center_s + cell->m_surfel->mean_), VisMesh::colorFromNormal(cell->m_surfel->normal_.normalized()), scene->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                    }
                    if constexpr ( showSemanticColored )
                        if ( cell->m_class )
                        s.addSurfel(cell->m_surfel->eigen_vectors_, cell->m_surfel->eigen_values_, cell->m_center_s + cell->m_surfel->mean_, m_label_manager->color_for_label(cell->m_class->getArgMaxClass()).cast<float>(), asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell->m_level)/2,m_surfel_type);
                }
                v.showSurfelMesh("sceneSurfels_"+std::to_string(idx), true, interp_scene_pose,m_surfel_type);
                n.showSurfelMesh("sceneSurfelsN_"+std::to_string(idx), true, interp_scene_pose,m_surfel_type);
                if constexpr ( showSemanticColored )
                    s.showSurfelMesh("sceneSurfelsSem_"+std::to_string(idx), true, interp_scene_pose,m_surfel_type);
                //n.showEdgeMesh("sceneSurfelNormals_"+std::to_string(idx), true, interp_scene_pose);
                //f.showEdgeMesh("sceneSurfelFirstViewDir_"+std::to_string(idx), true, interp_scene_pose);
                //e.showPointMesh("sceneSurfelNotValidMean_"+std::to_string(idx), true, interp_scene_pose);
                //p.showPointMesh("sceneSurfelTooLittlePoints_"+std::to_string(idx), true, interp_scene_pose);
            }
            LOG(1) << "\n\ntotal_scene_surfel_show: " << total.getTime();
        }
#endif

#ifdef USE_EASY_PBR
        const bool showSubdivSurfels = m_show_subdiv_surfels;
        if ( showSubdivSurfels )
        {
            MarsMap * last_scene = sceneMapPtrVec[0];

            const double plane_scale_factor = MapParameters::plane_scale_factor;
            const double first_ev_planarity_threshold = MapParameters::first_ev_planarity_threshold;
            const double second_ev_degenerate_threshold = MapParameters::second_ev_degenerate_threshold;

            std::vector<VisMesh> levelV(last_scene->m_map_params.m_num_levels);
            std::vector<VisMesh> levelC(last_scene->m_map_params.m_num_levels);
            std::vector<VisMesh> levelN(last_scene->m_map_params.m_num_levels);
            Sophus::SE3d pose;
            SurfelInfoVector scene_cells;
            last_scene->getCells( scene_cells, pose, false );
            for ( size_t idx = 0; idx < levelV.size(); ++idx )
            {
                levelV[idx].reserveSurfelTriangleEllipsoid(scene_cells.size());
                levelC[idx].reserveSurfelTriangleEllipsoid(scene_cells.size());
                levelN[idx].reserveSurfelTriangleEllipsoid(scene_cells.size());
            }
            for ( const SurfelInfo & cell : scene_cells )
            {
                if ( ! cell.m_surfel ) continue;
                if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 )
                {
//                    if ( cell.m_surfel->getNumPoints() > 10 )
//                    {
//                        levelV[cell.m_level].addSurfelTriangleEllipsoid( Eigen::Matrix3d::Identity(), Eigen::Vector3f::Constant(0.01*(4-cell.m_level)), cell.m_center_s, VisMesh::red,last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2);
//                        //if ( cell.m_center_s.norm() < 10 )
//                        LOG(1) << "lvl: " << cell.m_level << " c: " << cell.m_center_s.transpose() << " #pt: " << cell.m_surfel->getNumPoints() << " valid: " << cell.m_surfel->valid_ << " det: "  << cell.m_surfel->cov_.determinant() << " cov=[" << cell.m_surfel->cov_.row(0) << "; " << cell.m_surfel->cov_.row(1) << "; " << cell.m_surfel->cov_.row(2)<<"] ev: " << Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f>(cell.m_surfel->cov_).eigenvalues().transpose()<< " s: " << (Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f>(cell.m_surfel->cov_).eigenvalues().array()<= std::numeric_limits<float>::epsilon()).transpose();
//                    }
                    continue;
                }
                Eigen::Vector3f c = (cell.m_surfel->eigen_values_(0) < first_ev_planarity_threshold || cell.m_surfel->eigen_values_(0) / cell.m_surfel->eigen_values_(1) < plane_scale_factor) ? VisMesh::blue : VisMesh::yellow;
                Eigen::Vector3f r = (cell.m_surfel->eigen_values_(1) < second_ev_degenerate_threshold) ? VisMesh::red : VisMesh::yellow;
                levelV[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), c, last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
                levelC[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), r, last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
                levelN[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), VisMesh::colorFromNormal(cell.m_surfel->normal_.normalized()), last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
            }
            for ( size_t idx = 0; idx < levelV.size(); ++idx )
            {
                levelV[idx].showSurfelMesh("sdSceneL"+std::to_string(idx), true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
                levelC[idx].showSurfelMesh("sdSceneC"+std::to_string(idx), true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
                levelN[idx].showSurfelMesh("sdSceneN"+std::to_string(idx), true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
            }

            last_scene->getCells( scene_cells, pose );

            std::vector<VisMesh> levelOV(last_scene->m_map_params.m_num_levels);
            std::vector<VisMesh> levelOC(last_scene->m_map_params.m_num_levels);
            SurfelInfoVector ov;
            last_scene->getCellsAdaptive ( ov, pose);
            for ( size_t idx = 0; idx < levelOV.size(); ++idx )
            {
                levelOV[idx].reserveSurfelTriangleEllipsoid(ov.size());
                levelOC[idx].reserveSurfelTriangleEllipsoid(ov.size());
            }
            for ( const SurfelInfo & cell : ov )
            {
                if ( ! cell.m_surfel ) continue;
                if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 ) continue;
                const Eigen::Vector3f c = (cell.m_surfel->eigen_values_(0) < first_ev_planarity_threshold || cell.m_surfel->eigen_values_(0) / cell.m_surfel->eigen_values_(1) < plane_scale_factor) ? VisMesh::blue : VisMesh::yellow;
                const Eigen::Vector3f r = (cell.m_surfel->eigen_values_(1) < second_ev_degenerate_threshold) ? VisMesh::red : VisMesh::yellow;
                levelOV[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), c, last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
                levelOC[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), r, last_scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
            }
            for ( size_t idx = 0; idx < levelOV.size(); ++idx )
            {
                levelOV[idx].showSurfelMesh("osdSceneL"+std::to_string(idx), true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
                levelOC[idx].showSurfelMesh("osdSceneL"+std::to_string(idx), true, m_localMarsMap->local_map_pose * new_scene_pose,m_surfel_type);
            }
            LOG(1) << "Normal: " << scene_cells.size() << " Adapted: " << ov.size();

            std::vector<MarsMap *> sceneMapPtrVec = m_sceneMarsMap->getTransformedMapPtrVec(prev_scene_pose);
            for ( size_t idx = 0; idx < sceneMapPtrVec.size(); ++idx)
            {
                continue;
                MarsMap * scene = sceneMapPtrVec[idx];
                const Sophus::SE3d cloud_interp_pose = m_trajectory_spline->pose( newCurTimes(idx) ); // curTimes still contains second
                const Sophus::SE3d interp_scene_pose = m_localMarsMap->local_map_pose * prev_scene_pose * cloud_interp_pose;
                std::vector<VisMesh> levelV(scene->m_map_params.m_num_levels);
                std::vector<VisMesh> levelR(scene->m_map_params.m_num_levels);
                Sophus::SE3d pose;
                SurfelInfoVector scene_cells;
                scene->getCells( scene_cells, pose, false );
                for ( size_t idx = 0; idx < levelV.size(); ++idx )
                {
                    levelV[idx].reserveSurfelTriangleEllipsoid(scene_cells.size());
                    levelR[idx].reserveSurfelTriangleEllipsoid(scene_cells.size());
                }
                for ( const SurfelInfo & cell : scene_cells )
                {
                    if ( ! cell.m_surfel ) continue;
                    if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                    if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 ) continue;
                    const Eigen::Vector3f c = (cell.m_surfel->eigen_values_(0) < first_ev_planarity_threshold || cell.m_surfel->eigen_values_(0) / cell.m_surfel->eigen_values_(1) < plane_scale_factor) ? VisMesh::blue : VisMesh::yellow;
                    const Eigen::Vector3f r = (cell.m_surfel->eigen_values_(1) < second_ev_degenerate_threshold) ? VisMesh::red : VisMesh::yellow;
                    levelV[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), c, scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
                    levelR[cell.m_level].addSurfel( cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, (cell.m_center_s + cell.m_surfel->mean_), r, scene->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
                }
                for ( size_t vIdx = 0; vIdx < levelV.size(); ++vIdx)
                {
                    levelV[vIdx].showSurfelMesh("sdSceneS"+std::to_string(idx)+std::to_string(vIdx), true, interp_scene_pose,m_surfel_type);
                    levelR[vIdx].showSurfelMesh("sdSceneS"+std::to_string(idx)+std::to_string(vIdx), true, interp_scene_pose,m_surfel_type);
                }
            }
            if ( false )
            {
            std::vector<VisMesh> levelV(asOneCenteredLocalMap->m_map_params.m_num_levels);
            Sophus::SE3d pose;
            SurfelInfoVector scene_cells;
            asOneCenteredLocalMap->getCells( scene_cells, pose, false );
            for ( size_t vIdx = 0; vIdx < levelV.size(); ++vIdx)
                levelV[vIdx].reserveSurfelTriangleEllipsoid(scene_cells.size());
            for ( const SurfelInfo & cell : scene_cells )
            {
                if ( ! cell.m_surfel ) continue;
                if ( ! cell.m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                if ( ! cell.m_surfel->valid_ || cell.m_surfel->getNumPoints() < 10.0 ) continue;
                const Eigen::Vector3f c = (cell.m_surfel->eigen_values_(0) < first_ev_planarity_threshold || cell.m_surfel->eigen_values_(0) / cell.m_surfel->eigen_values_(1) < plane_scale_factor) ? VisMesh::blue : VisMesh::yellow;
                levelV[cell.m_level].addSurfel(cell.m_surfel->eigen_vectors_, cell.m_surfel->eigen_values_, cell.m_center_s + cell.m_surfel->mean_, c, asOneCenteredLocalMap->m_map_params.getCellSizeOnLevel(cell.m_level)/2,m_surfel_type);
            }
            for ( size_t vIdx = 0; vIdx < levelV.size(); ++vIdx)
                levelV[vIdx].showSurfelMesh("sdMars"+std::to_string(vIdx), true, m_localMarsMap->local_map_pose,m_surfel_type);
            }
        }
#endif

        prev_scene_pose = Sophus::SE3d(); // for now, need to change this once we move the whole map.

        last_time = old_last_scene_cloud->m_stamp;
        last_interp_pose = interp_pose;
    }
    last_pose = m_cur_pose;
    last_pos = m_cur_pos;

    return m_localMarsMap->local_map_pose * newest_cloud_pose;
}

MarsMapPointCloud::Ptr MarsSplineRegistrator::getLocalMapSurfelCenters()
{
    MarsMapPointCloud::Ptr local_map_centers = MarsMapPointCloud::create();
    SurfelInfoConstPtrVector surfels;
    Sophus::SE3d pose;
    {
        std::unique_lock<std::mutex> lock ( m_localMapMutex );
        MarsMapType::Ptr local_map = m_localMarsMap->getAsOneCenteredMap();
        local_map->getCells(surfels, pose);
        pose = getMapPose();
        local_map_centers->resize(surfels.size());
        for ( const SurfelInfoConstPtr & surfel : surfels )
        {
            local_map_centers->addPoint(surfel->m_center_s + surfel->m_surfel->mean_);
        }
    }
    local_map_centers->m_points = (pose.so3().template cast<float>().matrix() * local_map_centers->m_points).colwise()+pose.translation().template cast<float>();
    return local_map_centers;
}
MarsMapPointCloud::Ptr MarsSplineRegistrator::getSceneMapSurfelCenters()
{
    MarsMapPointCloud::Ptr scene_map_centers = MarsMapPointCloud::create();
    std::vector<Sophus::SE3d> scene_poses;
    Sophus::SE3d map_pose;
    std::vector<MarsMapPointCloud::Ptr> one_scene_map_centers;
    {
        std::unique_lock<std::mutex> lock ( m_sceneMapMutex );
        scene_poses = getSplinePoses();
        const size_t num_scenes = m_scene_surfels.size();
        one_scene_map_centers.resize(num_scenes,nullptr);
        for ( size_t scene_idx = 0; scene_idx < num_scenes; ++scene_idx )
        {
            one_scene_map_centers[scene_idx] = MarsMapPointCloud::create();
            one_scene_map_centers[scene_idx]->resize(m_scene_surfels[scene_idx].size());
            for ( const SurfelInfo & surfel : m_scene_surfels[scene_idx] )
            {
                if ( surfel.m_surfel )
                    one_scene_map_centers[scene_idx]->addPoint(surfel.m_center_s + surfel.m_surfel->mean_);
            }
        }
    }
    if ( one_scene_map_centers.empty() || scene_poses.empty() ) return scene_map_centers; // happens on the first scan
    {
        std::unique_lock<std::mutex> lock ( m_localMapMutex );
        map_pose = getMapPose();
    }
    const size_t num_scenes = one_scene_map_centers.size();
    for ( size_t scene_idx = 0; scene_idx < num_scenes; ++scene_idx )
    {
        const Sophus::SE3d pose = map_pose * scene_poses[scene_idx];
        one_scene_map_centers[scene_idx]->transform_points(pose);
        one_scene_map_centers[scene_idx]->setPose(Sophus::SE3d());
        scene_map_centers->add(one_scene_map_centers[scene_idx]);
    }
    return scene_map_centers;
}

Eigen::Vector6d MarsSplineRegistrator::getOdomVelocity() const
{
    if ( m_scene_velocity.empty() )
        return Eigen::Vector6d::Zero();
    else
        return m_scene_velocity.back();
}

std::vector<Sophus::SE3d> MarsSplineRegistrator::getSplinePoses() const
{
    return m_scene_poses;
}

Sophus::SE3d MarsSplineRegistrator::getInterpolatedSplinePose() const
{
    return m_interpolated_scene_pose;
}

Sophus::SE3d MarsSplineRegistrator::getMapPose() const
{
    return m_first_pose_world_sensor * m_localMarsMap->local_map_pose;
}

void MarsSplineRegistrator::setFirstPose ( const Sophus::SE3d & first_pose_world_sensor )
{
    m_first_pose_world_sensor = first_pose_world_sensor;
}
void MarsSplineRegistrator::setFirstOrientation ( const Eigen::Quaterniond & first_ori_world_sensor )
{
    m_first_ori_world_sensor = first_ori_world_sensor;
}
