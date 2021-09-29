/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration
Copyright (c) 2016-2021, Computer Science Institute VI, University of Bonn.

This file contains Code adapted from the Basalt project.
https://gitlab.com/VladyslavUsenko/basalt-headers.git
Copyright (c) 2019, Vladyslav Usenko and Nikolaus Demmel.

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

#include "MarsSplineRegistration.h"

#include "MarsAssociations.h"
#include "MarsAssociator.h"
#include "basalt/spline/se3_spline.h"
#include "MarsMap.h"

#define LOGURU_REPLACE_GLOG 1
#define LOGURU_WITH_STREAMS 1
#include <loguru.hpp>

#include <Tracy.hpp>
#include "StopWatch.h"

#include <deque>
#include <fstream>
#include <numeric>

typedef float NormalEquationScalarType;
constexpr bool print_info = false;
//constexpr int N = MarsSplineType::N;
constexpr int D = Sophus::SE3d::DoF;
//constexpr int ND = N*D;
typedef Eigen::Matrix<NormalEquationScalarType,3,3> Mat3;
typedef Eigen::Matrix<NormalEquationScalarType,6,6> Mat6;
typedef Eigen::Matrix<NormalEquationScalarType,3,6> Mat36;
//typedef Eigen::Matrix<NormalEquationScalarType,3,ND> Mat3K;
//typedef Eigen::Matrix<NormalEquationScalarType,ND,3> MatK3;
//typedef Eigen::Matrix<NormalEquationScalarType,ND,Eigen::Dynamic> MatKX;
//typedef Eigen::Matrix<NormalEquationScalarType,ND,ND> MatK;
//typedef Eigen::Matrix<NormalEquationScalarType,ND,1> VecK;

using MatD = Eigen::Matrix<NormalEquationScalarType,D,D>;
using VecD = Eigen::Matrix<NormalEquationScalarType,D,1>;
//typedef Eigen::Matrix<NormalEquationScalarType,D,D> MatD;
//typedef Eigen::Matrix<NormalEquationScalarType,D,1> VecD;
//typedef Eigen::Matrix<NormalEquationScalarType,D,ND> MatDK;

struct MarsRegistrationDataPair
{
    MarsMap * m_model = nullptr;
    MarsMap * m_scene = nullptr;
    Sophus::SE3d* m_transform = nullptr;
    size_t m_model_num_points = 0, m_scene_num_points = 0;
    SurfelInfoConstPtrVector m_scene_cells;
};

typedef ModelAssociationData ModelAssociationT;

template <typename Scalar>
inline void updateSingleNormalEquations3( const Eigen::Matrix<Scalar,3,3> & wW,
                                         const Mat36 & Dx,
                                         MatD & accumC )
{
    accumC.noalias() += Dx.transpose() * wW.template cast<NormalEquationScalarType>() * Dx;
}

template <typename Scalar, int N>
inline void updateSingleNormalEquations2( const Eigen::Matrix<Scalar,3,3> & wW,
                                         const Eigen::Matrix<Scalar,3,1> & wdiff,
                                         const Eigen::Matrix<NormalEquationScalarType,3,N*D> & Dx,
                                         Eigen::Matrix<NormalEquationScalarType,N*D,N*D> & accumHL,
                                         Eigen::Matrix<NormalEquationScalarType,N*D,1> & accumBL )
{
    accumBL.noalias() += Dx.transpose() * wdiff.template cast<NormalEquationScalarType>();
    // calculate the right upper blocks
    for ( int i = 0; i < N; ++i )
    {
        const Eigen::Matrix<NormalEquationScalarType,D,3> wDtW = Dx.template block<3,D>(0,D*i).transpose() * wW.template cast<NormalEquationScalarType>();
        for ( int j = i; j < N; ++j )
        {
            // upper:
            accumHL.template block<D,D>(D*i, D*j).noalias() += wDtW * Dx.template block<3,D>(0,D*j);
        }
    }
}

template <typename ModelAssociationT>
inline void associate ( MarsMap * model, const SurfelInfoConstPtrVector& sceneCells,
                        std::vector<SceneAssociationData> & scene_associations,
                        std::vector<ModelAssociationT> & model_associations,
                        const Sophus::SE3d & transform, const int & neighbors )
{
    StopWatch watch;
    scene_associations.reserve( sceneCells.size());
    model_associations.reserve(27 * sceneCells.size());
    if ( print_info ) LOG(INFO) << "surfel_registration_debug: associateMaps size() " << sceneCells.size() << " tf:\n"<<transform.params().transpose();
    ZoneScopedN("MarsSplineRegistration::Associate");
    int num_valid_scene = 0;
    int num_valid_model = 0;
    MarsAssociator<ModelAssociationT> af ( model, &sceneCells, &scene_associations, &model_associations, &num_valid_scene, &num_valid_model, transform, neighbors );
    af.for_each();
    if ( print_info ) LOG(INFO) << "surfel_registration_timing: AssociateFunctor took: " << watch.getTime() << " for " << model_associations.size() << " associations";
}

#include "MarsSplineRegistrationOMP.hpp"

struct SceneAssociationsVec
{
    explicit SceneAssociationsVec(const size_t & s)
    {
        m_vec.resize(s);
    }
    std::vector<std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>>> m_vec;
};

MarsSplineRegistration::MarsSplineRegistration()
{
}

bool MarsSplineRegistration::getCovariance( const int & scene_num, const Sophus::SE3d & pose_ws, Eigen::Matrix6d & covar ) const
{
    if ( !m_prev_scene_assocs || scene_num < 0 || scene_num >= m_prev_scene_assocs->m_vec.size() )
        return false;
    covar.setZero();
    const std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>> & associations = m_prev_scene_assocs->m_vec[scene_num];
    Mat6 info = Mat6::Zero();

    const std::vector<SceneAssociationData> & sceneAssocs = associations.first;
    const std::vector<ModelAssociationT> & modelAssocs = associations.second;
    MarsInfoWeightor<ModelAssociationT,1> weightor ( &sceneAssocs, &modelAssocs, &info, pose_ws, m_model_num_points, m_scene_num_points[scene_num]);
    weightor.info();

    const Eigen::Matrix6d sym_info = info.cast<double>().selfadjointView<Eigen::Upper>();
    covar = sym_info.inverse();

    //LOG(INFO) << "info:" << info.determinant() <<"\n" << info;
    //LOG(INFO) << "inf2:" << sym_info.determinant() <<"\n" << sym_info;
    //LOG(INFO) << "cov:" << covar.determinant() <<"\n" << covar;
    //if (! covar.allFinite() ) LOG(FATAL) << "cov not finite !";
    return covar.allFinite();
}

bool MarsSplineRegistration::estimate( MarsMap* model, const std::vector<MarsMap*> & sceneVec, const Eigen::VectorXt & times, MarsSplineTypePtr & spline, const int & maxIterations )
{
    //constexpr bool print_info = true;
    if ( print_info ) LOG(INFO) << "estimateWindowedSplineLM: sv: " << sceneVec.size();
    StopWatch watch;

    constexpr int N = MarsSplineType::N;
    using MatK = Eigen::Matrix<NormalEquationScalarType,N*D,N*D>;
    //using Mat3K = Eigen::Matrix<NormalEquationScalarType,3,N*D>;
    using VecK = Eigen::Matrix<NormalEquationScalarType,N*D,1>;

    std::vector<Sophus::SE3d> currentTransforms ( sceneVec.size(), Sophus::SE3d() );
    std::vector<MarsSplineType::PosePosSO3JacobianStruct> currentJac ( sceneVec.size() );

    // set up the minimization algorithm
    std::vector<MarsRegistrationDataPair> paramVec; paramVec.reserve(sceneVec.size() * sceneVec.size());
    std::unordered_map<size_t,size_t> modelToSceneIds;
    std::unordered_map<size_t,int> paramToModelVec;

    m_model_num_points = 0;

    std::vector<MarsMap*> sceneModellVec({model});
    sceneModellVec.insert(sceneModellVec.end(),sceneVec.begin(),sceneVec.end());
    std::vector<Sophus::SE3d> currentModellTransforms ( sceneModellVec.size(), Sophus::SE3d() );

    std::vector<MarsMap::Ptr> transformedSceneMaps(sceneModellVec.size());

    constexpr bool dontFullyConnect = true;
    constexpr bool onlyConnectSeq = false;
    constexpr bool onlyConnectToModel = true;

    size_t paramIdx = 0;
    {
    ZoneScopedN("MarsSplineRegistration::get_prep");
    for ( size_t otherSceneIdx = 0; otherSceneIdx < sceneModellVec.size(); ++otherSceneIdx )
    {
        if ( !sceneModellVec[otherSceneIdx] ) continue;
        MarsMap * curModel = sceneModellVec[otherSceneIdx];
        MarsMap * transformed_scene = curModel;
        if ( otherSceneIdx > 0 && ! onlyConnectToModel ) // not the original model!
        {
            currentModellTransforms[otherSceneIdx] = spline->pose ( times(otherSceneIdx-1) );
            transformedSceneMaps[otherSceneIdx] = MarsMap::create();
            transformed_scene = transformedSceneMaps[otherSceneIdx].get();
            transformed_scene->initParams ( *curModel );

            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            curModel->getCells( scene_cells, pose );

            Sophus::SE3d pose_w2 = transformed_scene->getMapPose() * currentModellTransforms[otherSceneIdx] * (model->getMapPose().inverse() * curModel->getMapPose());
            transformed_scene->addCells ( scene_cells, pose_w2 );
            transformed_scene->setPose ( model->getMapPose() );

            if (print_info) LOG(INFO) << "Oidx: " << otherSceneIdx << " modelId: " << model->m_id << " curModel: " << curModel->m_id << " TransformedScene: " << transformed_scene->m_id;
            if (print_info) LOG(INFO) << "OSI: " << otherSceneIdx << " sc: " << scene_cells.size()<<"\nT:\n"<<currentModellTransforms[otherSceneIdx].params().transpose() << "\nP:\n"<<curModel->getMapPose().params().transpose()
                      << "\nP2:\n"<<pose_w2.params().transpose() << "\nPP2:\n"<<(transformed_scene->getMapPose().inverse()*pose_w2).params().transpose();
            if (print_info) LOG(INFO) << "tp:\n" << transformed_scene->getMapPose().params().transpose() << "\npw:\n" << pose_w2.params().transpose() << "\nps1s2:\n"<< (transformed_scene->getMapPose().inverse() * pose_w2).params().transpose();
        }
        const size_t num_scenes = sceneVec.size();
        for ( size_t sceneIdx = 0; sceneIdx < num_scenes; ++sceneIdx )
        {
            if ( sceneIdx+1 == otherSceneIdx ) continue;
            if ( dontFullyConnect && sceneIdx+1 < otherSceneIdx ) continue;
            if ( onlyConnectSeq && sceneIdx != otherSceneIdx ) continue;
            if ( !sceneVec[sceneIdx] ) continue;
            MarsMap * scene = sceneVec[sceneIdx];

            paramVec.emplace_back();
            modelToSceneIds[paramIdx] = sceneIdx;
            paramToModelVec[paramIdx] = otherSceneIdx-1;
            MarsRegistrationDataPair & params = paramVec[paramIdx];
            Sophus::SE3d & currentTransform = currentTransforms[sceneIdx];
            //LOG(INFO) << "times: " << times.transpose() << " dt: " << spline->getDtNs() << " minT: " << spline->minTimeNs() << " maxT: " << spline->maxTimeNs() ;
            //spline->print_knots();
            currentTransform = spline->pose ( times(sceneIdx) );
            params.m_scene = scene;
            params.m_model = transformed_scene;
            if (print_info) LOG(INFO) << "surfel_registration_debug: current_transform["<< sceneIdx<<"]=["<<currentTransform.params().transpose()<<"]";

            Sophus::SE3d pose;
            SurfelInfoConstPtrVector scene_cells;
            StopWatch newWatch;
            if ( m_use_adaptive )
                scene->getCellsAdaptive( scene_cells, pose );
            else
                scene->getCells( scene_cells, pose );
            if (print_info) LOG(INFO) << "surfel_registration_timing: getting cells took : " << newWatch.getTime() << " got cells: " << scene_cells.size();
            params.m_scene_cells.reserve(scene_cells.size());
            int num_invalid = 0, num_less_pts = 0;
            for ( const SurfelInfoConstPtr & cell : scene_cells )
            {
                if ( ! cell->m_surfel ) continue;
                if ( ! cell->m_surfel->evaluated_ ) LOG(FATAL) << "this should not have happend!";
                num_invalid += !cell->m_surfel->valid_;
                num_less_pts += (cell->m_surfel->getNumPoints() < 10.0);
                if ( !cell->m_surfel->valid_ || cell->m_surfel->getNumPoints() < 10.0 ) continue;
                params.m_scene_cells.emplace_back(cell);
            }
            if ( print_info ) LOG(INFO) << "SceneCells: " << scene_cells.size() << " invalid: " << num_invalid << " less: " << num_less_pts << " valid: " << params.m_scene_cells.size();

            if ( params.m_scene_cells.empty() )
            {
                LOG(INFO)  << "MarsSplineRegistration::estimateTransformationLM(): no valid surfels in scene! Forgot to evaluate()?";
            }
            params.m_scene_num_points = scene->getNumPoints();

            if (print_info) LOG(INFO) << "surfel_registration_timing: estimateTransformationLM took : " << watch.getTime() << " ( " << newWatch.getTime()<< " ) "
                      << " after scene cells with " << scene_cells.size() << " cells and " << params.m_scene_cells.size() << " valid surfels.";
            params.m_model_num_points = transformed_scene->getNumPoints();
            params.m_transform = &currentTransform;
            ++paramIdx;
            m_model_num_points = params.m_model_num_points;
            m_scene_num_points.emplace_back(params.m_scene_num_points);
        }
        if ( onlyConnectToModel ) break;
    }
    }
    const size_t num_params = paramVec.size();
    m_prev_scene_assocs = std::make_shared<SceneAssociationsVec>(num_params);
    std::vector<std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>>> & associations = m_prev_scene_assocs->m_vec;
    Eigen::VectorXi num_cells = Eigen::VectorXi::Zero(num_params,1);
    for ( size_t paramIdx = 0; paramIdx < num_params; ++paramIdx )
    {
        const MarsRegistrationDataPair & params = paramVec[paramIdx];
        std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>> & assocs = associations[paramIdx];
        assocs.first.resize(params.m_scene_cells.size());
        assocs.second.resize(params.m_scene_cells.size()*27);
        num_cells(paramIdx) = params.m_scene_cells.size();
    }
    LOG(INFO) << "registering num_cells: " << num_cells.transpose() << " map: " << model->m_surfels.size();

    bool retVal = true;

    constexpr double init_lambda = 1e-6;
    constexpr double min_lambda = 1e-18;
    constexpr double max_lambda= 100;
    constexpr double stop_thresh = 1e-5;
    constexpr double err_thresh = 1e-4;
    double lambda_vee = 2.;
    double lambda = init_lambda;
    double old_error = std::numeric_limits<double>::max();
    MatK accumH = MatK::Zero();
    VecK accumB = VecK::Zero();
    double opt_error = 0;

    constexpr int num_neighbors = 1;
    int iter = 0;

    while (iter < maxIterations)
    {
        if (print_info) LOG(INFO) << "iteration " << iter << " last error: " << old_error;
        const size_t num_params = paramVec.size();
        for ( size_t paramIdx = 0; paramIdx < num_params; ++paramIdx )
        {
            // relative pose = model pose inverse * scene pose
            // err = model pt - T * relative pose * scene pt
            // new: model pt - relative pose * T * scene pt
            const size_t & sceneIdx = modelToSceneIds[paramIdx];
            MarsSplineType::PosePosSO3JacobianStruct & jac = currentJac[sceneIdx];
            currentTransforms[sceneIdx] = spline->pose ( times(sceneIdx), &jac );
            if (print_info) LOG(INFO) << "paramIdx: " << paramIdx << " sceneIdx: " << sceneIdx << " time: "<< times(sceneIdx) << " pose: " << currentTransforms[sceneIdx].params().transpose();
        }

        {
            ZoneScopedN("MarsSplineRegistration::estimate::associateScenesToModel");
            associateScenesToModel(paramVec, associations, modelToSceneIds, currentTransforms, num_neighbors );
        }

        accumH.setZero();
        accumB.setZero();
        NormalEquationScalarType accumE = 0;
        NormalEquationScalarType accumW = 0;
        //retVal = true
        int num_valid_assocs = 0;
        {
            ZoneScopedN("MarsSplineRegistration::estimate::setupNormalEquations");
            num_valid_assocs = setUpNormalEquations<MarsSplineType>(paramVec,associations,modelToSceneIds,currentTransforms,currentJac,accumH,accumB,accumE,accumW);
        }
        if ( num_valid_assocs == 0 ) LOG(FATAL) << "Could not establish valid associations?";

        opt_error = accumE;
        //if (print_info) LOG(INFO) << "surfel_registration_debug: ret: " << retVal << " err: " << accumE << " took: " << watch.getTime();


        // following Code is modified from Basalt Project
        bool converged = false;
        bool step = false;
        int remaining_step_iter = 10;
        VecK Hdiag = accumH.diagonal();
        double max_inc = std::numeric_limits<double>::max();
        {
        ZoneScopedN("MarsSplineRegistration::solve_update");
        while (!step && remaining_step_iter > 0 && !converged) {
          ZoneScopedN("MarsSplineRegistration::solve_one_iter_step");

          VecK Hdiag_lambda = Hdiag * lambda;
          for (int i = 0; i < Hdiag_lambda.size(); ++i)
            Hdiag_lambda[i] = std::max<NormalEquationScalarType>(Hdiag_lambda[i], min_lambda);

          const MatK Hd = Hdiag_lambda.asDiagonal();
          const VecK inc_full = -(accumH+Hd).selfadjointView<Eigen::Upper>().ldlt().solve(accumB);
          max_inc = inc_full.array().abs().maxCoeff();
          if (max_inc < stop_thresh) converged = true;

          MarsSplineTypePtr new_spline = std::make_shared<MarsSplineType>( spline->getDtNs(), spline->minTimeNs() );
          *new_spline = *spline;
          const size_t num_knots = new_spline->numKnots();
          for (size_t knot_idx = 0; knot_idx < num_knots; ++knot_idx) {
            const Eigen::Vector6d inc = inc_full.template segment<D>(D * knot_idx).template cast<double>();
            if (print_info) LOG(INFO) << " remaining_step_iter : " << remaining_step_iter << " knotInc["<<knot_idx<<"]=" << inc.transpose();
            if ( !inc.allFinite() ) LOG(FATAL) << "Increment was not valid: " << inc.transpose();
            new_spline->applyInc(knot_idx, inc);
          }

          NormalEquationScalarType newAccumE = 0;
          NormalEquationScalarType newAccumW = 0;
          {
              evalOnce<MarsSplineType> (paramVec, associations, modelToSceneIds, new_spline, times, newAccumE, newAccumW );
          }
          newAccumE /= newAccumW;

          const double new_spline_error = newAccumE;
          const double f_diff = (opt_error - new_spline_error);
          const double l_diff = 0.5 * inc_full.dot(inc_full * lambda - accumB);
          if ( std::abs(f_diff) < err_thresh ) converged = true;

          const double step_quality = f_diff / l_diff;
          if (step_quality < 0 ) {
              //if (print_info)
              //              LOG(INFO) << "it: " << iter << " [REJECTED] lambda:" << lambda
              //                        << " step_quality: " << step_quality
              //                        << " max_inc: " << max_inc << " Error: " << new_spline_error << " oldError: " << opt_error;
              lambda = std::min(max_lambda, lambda_vee * lambda);
              lambda_vee *= 2;
          } else {
              if (print_info)
                LOG(INFO) << "it: " << iter << " [ACCEPTED] lambda:" << lambda
                          << " step_quality: " << step_quality
                          << " max_inc: " << max_inc << " Error: " << new_spline_error << " oldError: " << opt_error;
              lambda = std::max( min_lambda, lambda * std::max(1.0 / 3, 1 - std::pow(2 * step_quality - 1, 3.0)));
              lambda_vee = 2;
              spline = new_spline;
              opt_error = new_spline_error;
              step = true;
          }
          --remaining_step_iter;
        }
        }
        // end

        if ( ! std::isfinite(opt_error) )
        {
            LOG(INFO) << "surfel_registration_debug: registration failed " << old_error << " new_error "<< opt_error;
            return false;
        }
        if (print_info)
            LOG(INFO) << "surfel_registration_debug: after update eval took : " << watch.getTime() << " last_error: " << old_error << " new_error: " << opt_error << " conv: " << converged << " errorDiff: " << std::abs(old_error - opt_error) << " maxInc: " << max_inc;
        if ( converged || std::abs(old_error - opt_error) < err_thresh || max_inc < stop_thresh )
        {
            old_error = opt_error;
            break;
        }
        old_error = opt_error;

        const size_t num_scenes = sceneVec.size();
        for ( size_t sceneIdx = 0; sceneIdx < num_scenes; ++sceneIdx )
        {
            MarsSplineType::PosePosSO3JacobianStruct & curJac = currentJac[sceneIdx];
            currentTransforms[sceneIdx] = spline->pose( times(sceneIdx), &curJac );
        }

        if ( ! onlyConnectToModel )
        {
            for ( size_t otherSceneIdx = 0; otherSceneIdx < sceneModellVec.size(); ++otherSceneIdx )
            {
                if ( !sceneModellVec[otherSceneIdx] ) continue;
                MarsMap * curModel = sceneModellVec[otherSceneIdx];
                if ( otherSceneIdx > 0 && ! onlyConnectToModel) // not the original model!
                {
                    currentModellTransforms[otherSceneIdx] = spline->pose( times(otherSceneIdx-1) );
                    MarsMap * transformed_scene = transformedSceneMaps[otherSceneIdx].get();
                    transformed_scene->clear();
                    transformed_scene->setPose ( model->getMapPose() );

                    Sophus::SE3d pose;
                    SurfelInfoConstPtrVector scene_cells;
                    if ( m_use_adaptive )
                        curModel->getCellsAdaptive(scene_cells, pose);
                    else
                        curModel->getCells( scene_cells, pose );
                    Sophus::SE3d pose_w2 = transformed_scene->getMapPose() * currentModellTransforms[otherSceneIdx] * (model->getMapPose().inverse() * curModel->getMapPose());
                    transformed_scene->addCells ( scene_cells, pose_w2 );
                }
            }
        }
        ++iter;
    }
    static float watch_max = watch.getTime();
    if ( watch.getTime() > watch_max ) watch_max = watch.getTime();
    if ( print_info )
    LOG(INFO) << "Timing: overall: "<< watch.getTime()<< " (max: " << watch_max << " ) ";
//              << " paramConnectTime: " << paramConnectTime << " iterPrepTime: "<< iterPrepTime
//              << " iterAssocTime: " << iterAssocTime << " iterGradTime: " << iterGradTime << " ( dk: " << iterDkTime << " cov: " << iterCovTime << ", gradMul: " << iterGradMulTime << " Dx: " << iterGradMulDxTime  << " Add: " << iterGradMulAddTime << " )"
//              << " iterLMTime: " << iterLMTime << " ( Cov: " << iterLMCovTime << " ) iterPostTime: " << iterPostTime << " ( spline: " << iterPostSplineTime << ", tf: " << iterPostTransformTime << ")";

    static int64_t sum_iters = 0;
    static int64_t num_estims = 0;
    ++num_estims;
    sum_iters += iter;

    //if ( print_info )
    LOG(INFO) << "surfel_registration_timing: estimateTransformationLM took: " << watch.getTime()
              << " after " << iter << " iterations. ( avg= " << (sum_iters/float(num_estims)) << " )";
    if (iter == maxIterations)
        if ( print_info ) LOG(INFO) << "estimateTransformationLM: maximum number of iterations (" << maxIterations << ") reached ";
    return std::isfinite(opt_error);
}

// following Code is modified from Basalt Project
void MarsSplineRegistration::optimizeSplineKnotsFromPoses ( MarsSplineTypePtr & spline, const Eigen::VectorXt & times, const std::vector<Sophus::SE3d> & poses ) const
{
    if ( print_info ) LOG(INFO) << "opt knots from poses. t: " << times.transpose();
    constexpr double init_lambda = 1e-6;
    constexpr double min_lambda = 1e-18;
    constexpr double max_lambda= 100;
    constexpr double stop_thresh = 1e-8;
    constexpr double err_thresh = 1e-6;
    constexpr int maxIter = 100;


    typedef double NormalEquationScalarType;
    constexpr int N = MarsSplineType::N;
    constexpr int D = Sophus::SE3d::DoF;
    constexpr int ND = N*D;
    using MatK = Eigen::Matrix<NormalEquationScalarType,ND,ND>;
    using VecK = Eigen::Matrix<NormalEquationScalarType,ND,1>;
    StopWatch watch;

    MatK accumH = MatK::Zero();
    VecK accumB = VecK::Zero();
    int iter = 0;
    double lambda = init_lambda;
    double lambda_vee = 2.;
    double accumE = 0;
    double old_error = std::numeric_limits<double>::max();
    std::vector<Sophus::SE3d> meas_pose_inv ( poses.size() );
    for ( size_t posesIdx = 0; posesIdx < poses.size(); ++posesIdx )
        meas_pose_inv[posesIdx] = poses[posesIdx].inverse();
    MarsSplineTypePtr new_spline = std::make_shared<MarsSplineType> ( spline->getDtNs(), spline->minTimeNs() );
    new_spline->setKnots( spline->getKnot(0), N );
    for (size_t knot_idx = 0; knot_idx < new_spline->numKnots(); ++knot_idx)
        new_spline->setKnot(spline->getKnot(knot_idx),knot_idx);
    MarsSplineType::PosJacobianStruct Jp;
    MarsSplineType::SO3JacobianStruct Jr;
    Eigen::Matrix<double,D,ND> Dx;

    for ( iter = 0; iter < maxIter; ++iter )
    {
        accumH.setZero();
        accumB.setZero();
        accumE = 0;
        Dx.setZero();
        Eigen::Matrix<double,D,1> diff;
        for ( size_t posesIdx = 0; posesIdx < poses.size(); ++posesIdx )
        {
            const int64_t & t = times[posesIdx];
            const Sophus::Vector3d diff_pos = spline->positionResidual( t, poses[posesIdx].translation(), &Jp );
            const Sophus::Vector3d diff_so3 = spline->orientationResidual( t, poses[posesIdx].so3(), &Jr );
            diff.head<3>() = diff_pos;
            diff.tail<3>() = diff_so3;
            for (int i = 0; i < N; ++i) {
                Dx.template block<3,3>(0,D*i) = Jp.d_val_d_knot[i] * Eigen::Matrix3d::Identity(); // pos
                Dx.template block<3,3>(3,D*i+3) = Jr.d_val_d_knot[i]; // rot
            }
            const double error = .5 * diff.dot(diff);
            accumH += (Dx.transpose()* Dx);
            accumB += (Dx.transpose()* diff);
            accumE += error;
        }

        bool converged = accumE < err_thresh;
        bool step = false;
        int remaining_step_iter = 10;
        Eigen::VectorXd Hdiag = accumH.diagonal();
        double max_inc = std::numeric_limits<double>::max();
        while (!step && remaining_step_iter > 0 && !converged) {
            Eigen::VectorXd Hdiag_lambda = Hdiag * lambda;
            for (int i = 0; i < Hdiag_lambda.size(); ++i)
                Hdiag_lambda[i] = std::max(Hdiag_lambda[i], min_lambda);

            MatK Hd = Hdiag_lambda.asDiagonal();
            Eigen::VectorXd inc_full = -(accumH+Hd).ldlt().solve(accumB);
            max_inc = inc_full.array().abs().maxCoeff();
            if (max_inc < stop_thresh) converged = true;

            for (size_t knot_idx = 0; knot_idx < new_spline->numKnots(); ++knot_idx) {
                new_spline->setKnot(spline->getKnot(knot_idx),knot_idx);
                Eigen::Vector6d inc = inc_full.template segment<D>(D * knot_idx);
                if ( print_info ) LOG(INFO) << " remaining_step_iter : " << remaining_step_iter << " knotInc["<<knot_idx<<"]=" << inc.transpose();
                new_spline->applyInc(knot_idx, inc);
            }

            double newAccumE = 0;
            Eigen::Matrix<double,D,1> diff;
            for ( size_t posesIdx = 0; posesIdx < poses.size(); ++posesIdx )
            {
                const int64_t & t = times[posesIdx];
                const Sophus::Vector3d diff_pos = new_spline->positionResidual( t, poses[posesIdx].translation() );
                const Sophus::Vector3d diff_so3 = new_spline->orientationResidual( t, poses[posesIdx].so3() );
                diff.head<3>() = diff_pos;
                diff.tail<3>() = diff_so3;
                if ( print_info ) LOG(INFO) << "newErr["<<posesIdx<<"]: " << diff.transpose();
                const double error = .5 * diff.dot(diff);
                newAccumE += error;
            }

            const double f_diff = (accumE - newAccumE);
            const double l_diff = 0.5 * inc_full.dot(inc_full * lambda - accumB);
            if ( std::abs(f_diff) < err_thresh ) converged = true;

            const double step_quality = f_diff / l_diff;
            if (step_quality < 0 ) {
                //if (print_info)
                //              LOG(INFO) << "it: " << iter << " [REJECTED] lambda:" << lambda << " step_quality: " << step_quality
                //                        << " max_inc: " << max_inc << " Error: " << newAccumE << " oldError: " << accumE;
                lambda = std::min(max_lambda, lambda_vee * lambda);
                lambda_vee *= 2;
            } else {
                if (print_info)
                LOG(INFO) << "it: " << iter << " [ACCEPTED] lambda:" << lambda << " step_quality: " << step_quality
                          << " max_inc: " << max_inc << " Error: " << newAccumE << " oldError: " << accumE;
                lambda = std::max( min_lambda, lambda * std::max(1.0 / 3, 1 - std::pow(2 * step_quality - 1, 3.0)));
                lambda_vee = 2;
                // set spline
                for (size_t knot_idx = 0; knot_idx < new_spline->numKnots(); ++knot_idx)
                    spline->setKnot(new_spline->getKnot(knot_idx),knot_idx);
                accumE = newAccumE;
                step = true;
            }
            --remaining_step_iter;
        }
        if ( converged || std::abs(old_error - accumE) < err_thresh || max_inc < stop_thresh )
        {
            if ( print_info ) LOG(INFO) << "converged: " << converged << " abs: " <<  std::abs(old_error - accumE) << " maxInc: "<< max_inc; old_error = accumE;
            break;
        }
        old_error = accumE;
    }
    LOG(INFO) << "With PoseSpline needed " << iter << " iterations with err: " << old_error << " and took: " << watch.getTime() << " seconds";
}
