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

#include "MarsInfoWeightor.hpp"
#include "MarsEvaluator.hpp"

inline void associateScenesToModel ( std::vector<MarsRegistrationDataPair> & paramVec,
                                     std::vector<std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>>> & associations,
                                     std::unordered_map<size_t,size_t> & modelToSceneIds,
                                     const std::vector<Sophus::SE3d> & currentTransforms, const int & num_neighbors )
{
    ZoneScopedN("MarsSplineRegistration::reg_associate_all");
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( size_t paramIdx = 0; paramIdx < paramVec.size(); ++paramIdx )
    {
        std::vector<SceneAssociationData> & sceneAssocs = associations[paramIdx].first;
        std::vector<ModelAssociationT> & modelAssocs = associations[paramIdx].second;
        const MarsRegistrationDataPair & params = paramVec[paramIdx];
        const size_t & sceneIdx = modelToSceneIds[paramIdx];
        const Sophus::SE3d & currentTransform = currentTransforms[sceneIdx];
        if ( print_info ) LOG(INFO) <<"model has: "<< params.m_model->m_surfels.size() << " surfels";
        associate ( params.m_model, params.m_scene_cells, sceneAssocs, modelAssocs, currentTransform, num_neighbors );
    }
}

template <typename SplineType>
inline int  setUpNormalEquations ( std::vector<MarsRegistrationDataPair> & paramVec,
                                   std::vector<std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>>> & associations,
                                   std::unordered_map<size_t,size_t> & modelToSceneIds,
                                   const std::vector<Sophus::SE3d> & currentTransforms,
                                   const std::vector<typename SplineType::PosePosSO3JacobianStruct> & currentJac,
                                   Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,SplineType::N*D> & accumH,
                                   Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,1> & accumB,
                                   NormalEquationScalarType & accumE,
                                   NormalEquationScalarType & accumW )
{
    using MatK = Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,SplineType::N*D>;
    using Mat3K = Eigen::Matrix<NormalEquationScalarType,3,SplineType::N*D>;
    using VecK = Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,1>;
    ZoneScopedN("MarsSplineRegistration::get_system_one_iter");
    int accumN = 0;
#pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( size_t paramIdx = 0; paramIdx < paramVec.size(); ++paramIdx )
    {
        ZoneScopedN("MarsSplineRegistration::get_system_one_iter::once");

        NormalEquationScalarType accumEL = 0;
        NormalEquationScalarType accumWL = 0;
        MatK accumHL = MatK::Zero();
        VecK accumBL = VecK::Zero();
        int accumNL = 0;

        const std::vector<SceneAssociationData> & sceneAssocs = associations[paramIdx].first;
        const std::vector<ModelAssociationT> & modelAssocs = associations[paramIdx].second;
        const MarsRegistrationDataPair & params = paramVec[paramIdx];
        const size_t & sceneIdx = modelToSceneIds[paramIdx];
        const Sophus::SE3d & currentTransform = currentTransforms[sceneIdx];
        const typename SplineType::PosePosSO3JacobianStruct & curJac = currentJac[sceneIdx];

            Mat3K pDx = Mat3K::Zero();
            for (int i = 0; i < SplineType::N; ++i)
                pDx.template block<3,3>(0,D*i).noalias() = curJac.d_val_d_knot[i].template topLeftCorner<3,3>().template cast<NormalEquationScalarType>(); // rot added in weightor->

            //const Eigen::Matrix3d curRot_d_point_d_so3 = curRot * (-Sophus::SO3d::hat(sceneAssoc.m_scene_mean.cast<double>())); --> done in weightor
            for ( int i = 0; i < SplineType::N; ++i ) {
                pDx.template block<3,3>(0,D*i+3).noalias() = curJac.d_val_d_knot[i].template bottomRightCorner<3,3>().template cast<NormalEquationScalarType>();
            }

            MarsInfoWeightor<ModelAssociationT,SplineType::N> weightor ( &sceneAssocs, &modelAssocs, &accumHL, &accumBL, &accumWL, &accumEL, &accumNL, pDx, currentTransform, params.m_model_num_points, params.m_scene_num_points );
            weightor.for_each();

        #pragma omp critical
        {
            accumH += accumHL.template selfadjointView<Eigen::Upper>();
            accumB += accumBL;
            accumE += accumEL;
            accumW += accumWL;
            accumN += accumNL;
        }
    }
#ifdef DO_NORMALIZE
    accumH /= accumW;
    accumB /= accumW;
    accumE /= accumW;
#endif
    return accumN;
}

template <typename SplineType>
inline void evalOnce ( std::vector<MarsRegistrationDataPair> & paramVec,
                       std::vector<std::pair<std::vector<SceneAssociationData>,std::vector<ModelAssociationT>>> & associations,
                       std::unordered_map<size_t,size_t> & modelToSceneIds,
                       const std::shared_ptr<SplineType> & new_spline,
                       const Eigen::VectorXt & times,
                       NormalEquationScalarType & newAccumE,
                       NormalEquationScalarType & newAccumW )
{
    ZoneScopedN("MarsSplineRegistration::estimate::evalOnce");
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:newAccumE) reduction(+:newAccumW)
    for ( size_t paramIdx = 0; paramIdx < paramVec.size(); ++paramIdx )
    {
        NormalEquationScalarType newerAccumE = 0;
        NormalEquationScalarType newerAccumW = 0;

        const std::vector<SceneAssociationData> & sceneAssocs = associations[paramIdx].first;
        const std::vector<ModelAssociationT> & modelAssocs = associations[paramIdx].second;
        const MarsRegistrationDataPair & params = paramVec[paramIdx];
        const size_t & sceneIdx = modelToSceneIds[paramIdx];
        const Sophus::SE3d currentTransform = new_spline->pose( times(sceneIdx) );

        // accumulates the error too.
        MarsEvaluator<ModelAssociationT,NormalEquationScalarType> evaluator ( &sceneAssocs, &modelAssocs, &newerAccumW, &newerAccumE, currentTransform, params.m_model_num_points, params.m_scene_num_points );
        evaluator.for_each();

        // omp update via reduction
        newAccumE += newerAccumE;
        newAccumW += newerAccumW;
    }
#ifdef DO_NORMALIZE
    newAccumE /= newAccumW;
#endif
}

// TODO: weight according to data constraints direction?
template <typename SplineType, int M>
inline
bool setUpNormalEquationsSoftConstraints( const std::shared_ptr<SplineType> & spline,
                                          const Eigen::Matrix<int64_t,M,1> & times,
                                          const Eigen::Matrix<int,M,1> & used,
                                          const double & sca_weight,
                                          Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,SplineType::N*D> & accumH,
                                          Eigen::Matrix<NormalEquationScalarType,SplineType::N*D,1> & accumB,
                                          NormalEquationScalarType & accumE,
                                          NormalEquationScalarType & accumW )
{
    //LOG(FATAL) << "currently disabled ;-)";
    //constexpr double const_vel_reg_weight = 0.001;
    constexpr int min_num_constrained = 1;
    constexpr int ND = SplineType::N * D;
    using MatDK = Eigen::Matrix<NormalEquationScalarType,D,ND>;
    using VecD = Eigen::Matrix<NormalEquationScalarType,D,1>; // pos, ori
    using MatK = Eigen::Matrix<NormalEquationScalarType,ND,ND>;
    using VecK = Eigen::Matrix<NormalEquationScalarType,ND,1>;
    using AccelPosJac = typename SplineType::PosJacobianStruct;
    using AccelRotJac = typename SplineType::SO3JacobianStruct;

    MatK scaAccumH = MatK::Zero();
    VecK scaAccumB = VecK::Zero();

    bool any_valid = false;
    NormalEquationScalarType scaAccumE = 0, scaAccumW = 0;
    MatDK Dx = MatDK::Zero();
    VecD diff = VecD::Zero();
    AccelPosJac Ja;
    AccelRotJac Jr;

    for ( int idx = 0; idx < M; ++idx )
    {
        const int64_t t = times[idx];
        if ( used[idx] >= min_num_constrained ) continue;

        const Sophus::Vector3d diff_acc = spline->transAccelWorldResidual( t, Eigen::Vector3d::Zero(), &Ja );
        const Sophus::Vector3d diff_rot = spline->rotAccelBodyResidual( t, Eigen::Vector3d::Zero(), &Jr );
        diff.head<3>() = diff_acc.cast<NormalEquationScalarType>();
        diff.tail<3>() = diff_rot.cast<NormalEquationScalarType>();

        for (int i = 0; i < SplineType::N; ++i) {
            Dx.template block<3,3>(0,D*i) = (Ja.d_val_d_knot[i] * Eigen::Matrix3d::Identity()).template cast<NormalEquationScalarType>(); // pos
            Dx.template block<3,3>(3,D*i+3) = Jr.d_val_d_knot[i].template cast<NormalEquationScalarType>(); // rot
        }

        const double error = .5 * diff.dot(diff);
        scaAccumH += (Dx.transpose() * sca_weight * Dx);
        scaAccumB += (Dx.transpose() * sca_weight * diff);
        scaAccumE += sca_weight * error;
        scaAccumW += sca_weight;
        any_valid = true;
    }
    if constexpr ( print_info )
    LOG(INFO) << "SCA: ew: " << (scaAccumE/scaAccumW) << " w: " << scaAccumW << " e: " << scaAccumE << " av: " << any_valid;
    if ( any_valid )
    {
#ifdef DO_NORMALIZE
        accumH += (scaAccumH / scaAccumW).template triangularView<Eigen::Upper>();
        accumB += scaAccumB / scaAccumW;
        accumE += scaAccumE / scaAccumW;
#else
        accumH += (scaAccumH).template triangularView<Eigen::Upper>();
        accumB += scaAccumB;
        accumE += scaAccumE;
#endif
        accumW += scaAccumW;
    }
    return any_valid;
}

template <typename SplineType, int M>
inline
bool evalSoftConstraintsOnce ( const std::shared_ptr<SplineType> & new_spline,
                               const Eigen::Matrix<int64_t,M,1> & times,
                               const Eigen::Matrix<int,M,1> & used,
                               const double & sca_weight,
                               NormalEquationScalarType & newAccumE,
                               NormalEquationScalarType & newAccumW )
{
    constexpr int min_num_constrained = 1;

    bool any_valid = false;
    NormalEquationScalarType scaAccumE = 0, scaAccumW = 0;
    for ( int idx = 0; idx < M; ++idx )
    {
        const int64_t t = times[idx];
        if ( used[idx] >= min_num_constrained ) continue;
        const Sophus::Vector3d diff_acc = new_spline->transAccelWorldResidual( t, Eigen::Vector3d::Zero() );
        const Sophus::Vector3d diff_rot = new_spline->rotAccelBodyResidual( t, Eigen::Vector3d::Zero() );
        const double error_acc = .5 * diff_acc.dot(diff_acc);
        const double error_rot = .5 * diff_rot.dot(diff_rot);
        scaAccumE += sca_weight * (error_acc + error_rot);
        scaAccumW += sca_weight;
        any_valid = true;
    }
    if constexpr ( print_info )
    LOG(INFO) << "SCAA: ew: " << (scaAccumE/scaAccumW) << " w: "<< scaAccumW << " e: " << scaAccumE;
    if ( any_valid )
    {
#ifdef DO_NORMALIZE
        newAccumE += scaAccumE / scaAccumW;
#else
        newAccumE += scaAccumE;
#endif
        newAccumW += scaAccumW;
    }
    return any_valid;
}
