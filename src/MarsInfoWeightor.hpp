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

typedef ModelAssociationInfoOnly InfoType;

template<typename ModelAssociationT, int N>
class MarsInfoWeightor
{
public:
    using MatK = Eigen::Matrix<NormalEquationScalarType,N*D,N*D>;
    using Mat3K = Eigen::Matrix<NormalEquationScalarType,3,N*D>;
    using VecK = Eigen::Matrix<NormalEquationScalarType,N*D,1>;

    MarsInfoWeightor( const std::vector<SceneAssociationData> * scene_associations,
                      const std::vector<ModelAssociationT> * model_associations,
                      Mat6 * accumC,
                      const Sophus::SE3d & currTransform,
                      const size_t & num_points_model, const size_t & num_points_scene );

    MarsInfoWeightor( const std::vector<SceneAssociationData> * scene_associations,
                      const std::vector<ModelAssociationT> * model_associations,
                      MatK * accumH,
                      VecK * accumB,
                      NormalEquationScalarType * accum_W, NormalEquationScalarType * accum_E,
                      int * accum_N,
                      const Mat3K & pDx,
                      const Sophus::SE3d & currTransform,
                      const size_t & num_points_model, const size_t & num_points_scene );
    ~MarsInfoWeightor(){}

    void compute_info ( const SceneAssociationData& assoc ) const;
    void info();
    void for_each( );

    Sophus::SE3d m_currentTransform;
    Mat3 m_currentRelativeRotation;
    Mat3 m_currentRelativeRotationT;
    Mat3K m_pDx;
    Mat3 m_curT_d_point_d_t;

    NormalEquationScalarType m_soft_assoc2_PI_3 = std::numeric_limits<NormalEquationScalarType>::signaling_NaN();
    NormalEquationScalarType m_prior_prob_norm_factor = std::numeric_limits<NormalEquationScalarType>::signaling_NaN();
    NormalEquationScalarType m_soft_assoc_c1 = std::numeric_limits<NormalEquationScalarType>::signaling_NaN();
    NormalEquationScalarType m_soft_assoc_c2 = std::numeric_limits<NormalEquationScalarType>::signaling_NaN();
    NormalEquationScalarType m_soft_assoc_c3 = std::numeric_limits<NormalEquationScalarType>::signaling_NaN();

    Mat6 * m_accumC = nullptr;
    MatK * m_accumH = nullptr;
    VecK * m_accumB = nullptr;
    NormalEquationScalarType * m_accumW = nullptr;
    NormalEquationScalarType * m_accumE = nullptr;
    int * m_accumN = nullptr;
    std::vector<InfoType> * m_errors = nullptr;
    const std::vector<SceneAssociationData> * m_scene_associations = nullptr;
    const std::vector<ModelAssociationT> * m_model_associations = nullptr;
    inline void operator()( const SceneAssociationData& assoc ) const;
};



template<typename ModelAssociationT, int N>
MarsInfoWeightor<ModelAssociationT,N>::MarsInfoWeightor( const std::vector<SceneAssociationData> * scene_associations,
                                                         const std::vector<ModelAssociationT> * model_associations,
                                                         Mat6 * accumC,
                                                         const Sophus::SE3d & currTransform,
                                                         const size_t & num_points_model, const size_t & num_points_scene )
    : m_scene_associations ( scene_associations ), m_model_associations ( model_associations ), m_accumC( accumC ), m_currentTransform ( currTransform )
{
    m_currentRelativeRotation = m_currentTransform.so3().template cast<NormalEquationScalarType>().matrix();
    m_currentRelativeRotationT = m_currentRelativeRotation.transpose();

    m_curT_d_point_d_t = m_currentRelativeRotation.cast<NormalEquationScalarType>();

    m_soft_assoc2_PI_3 = MarsGMMParameters::soft_assoc_c2 * std::pow(M_PI,3);
    m_prior_prob_norm_factor = MarsGMMParameters::prior_prob / (1.0 - MarsGMMParameters::prior_prob) * (num_points_model / (NormalEquationScalarType) num_points_scene) * MarsGMMParameters::soft_assoc_c1;
}

template<typename ModelAssociationT, int N>
MarsInfoWeightor<ModelAssociationT,N>::MarsInfoWeightor( const std::vector<SceneAssociationData> * scene_associations,
                                                       const std::vector<ModelAssociationT> * model_associations,
                                                       MatK * accumH, VecK * accumB,
                                                       NormalEquationScalarType * accumW, NormalEquationScalarType * accumE, int * accumN,
                                                       const Mat3K & pDx,
                                                       const Sophus::SE3d & currTransform,
                                                       const size_t & num_points_model, const size_t & num_points_scene )
    : m_scene_associations ( scene_associations ), m_model_associations ( model_associations ), m_accumH( accumH ), m_accumB ( accumB ), m_accumW ( accumW ), m_accumE ( accumE ), m_accumN ( accumN ), m_currentTransform ( currTransform ), m_pDx ( pDx )
{
    m_currentRelativeRotation = m_currentTransform.so3().template cast<NormalEquationScalarType>().matrix();
    m_currentRelativeRotationT = m_currentRelativeRotation.transpose();

    for (int i = 0; i < N; ++i)
        m_pDx.template block<3,3>(0,D*i).noalias() = m_currentRelativeRotation.cast<NormalEquationScalarType>() * pDx.template block<3,3>(0,D*i);

    m_soft_assoc2_PI_3 = MarsGMMParameters::soft_assoc_c2 * std::pow(M_PI,3);
    m_prior_prob_norm_factor = MarsGMMParameters::prior_prob / (1.0 - MarsGMMParameters::prior_prob) * (num_points_model / (NormalEquationScalarType) num_points_scene) * MarsGMMParameters::soft_assoc_c1;
}

template<typename ModelAssociationT, int N>
inline void MarsInfoWeightor<ModelAssociationT,N>::operator()( const SceneAssociationData& assoc ) const
{
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar,3,3> Mat3;
    typedef Eigen::Matrix<Scalar,3,1> Vec3;

    const auto & mean_scene = assoc.m_scene_mean; // already in relative pose
    const Mat3 cov_scene = assoc.m_scene_cov.template cast<Scalar>() + Mat3(Vec3::Constant(0.0001 * mean_scene.norm()).asDiagonal()); // cov_scene.diagonal() += Eigen::Vector3d::Constant(0.0001 * mean_scene.norm());
    const Mat3 cov_scene_transformed = m_currentRelativeRotation.cast<Scalar>() * cov_scene * m_currentRelativeRotationT.cast<Scalar>();
    const Vec3 scene_mean_transformed = m_currentTransform.template cast<Scalar>() * mean_scene;

    // normalize weights
    Scalar sum_weight = m_prior_prob_norm_factor / sqrt( m_soft_assoc2_PI_3 * symDetAdd(cov_scene_transformed,Scalar(assoc.m_sigma2))); // determinant is equivalent to (cov_scene + diag(sigma2)) since rotation has det 1.

    const Scalar normal_sig = assoc.m_sigma2;
    const Scalar inv_normal_sig = 1./normal_sig;
    const Scalar nf = 1. / sqrt(2*M_PI*normal_sig);

    //ZoneScopedN("cov_weight_cells");

    const int start_idx = assoc.m_indices(0);
    const int end_idx = assoc.m_indices(1);
    const int errorStartIdx = 0;
    const int errorEndIdx = end_idx - start_idx;
    Eigen::Array<Scalar,Eigen::Dynamic,1> err (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> wes (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> dde (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> nde (errorEndIdx,1);
    for ( int idx = start_idx, errorIdx = errorStartIdx; idx < end_idx; ++idx, ++errorIdx  )
    {
        const ModelAssociationT & singleAssoc = (*m_model_associations)[idx];
        InfoType & error = (*m_errors)[errorIdx];
        const Mat3 cov_model = singleAssoc.m_model_cov.template cast<Scalar>();
        const Vec3 diff = scene_mean_transformed - singleAssoc.m_model_mean.template cast<Scalar>();
        const Mat3 cov = cov_model + cov_scene_transformed;
        error.m_W.noalias() = symInverse<Scalar>(cov).template cast<typename InfoType::Scalar>();
        error.m_diff.noalias() = error.m_W*diff.template cast<typename InfoType::Scalar>();
        err(errorIdx) = error.m_diff.dot(diff.template cast<typename InfoType::Scalar>());

        const Vec3 invcov_ps_diff = symInverseMultiplyDiagAdd(cov, Scalar(assoc.m_sigma2), diff);
        wes(errorIdx) = //singleAssoc.m_model_num_points *
                MarsGMMParameters::soft_assoc_c1 / sqrt(m_soft_assoc2_PI_3 * symDetAdd(cov,Scalar(assoc.m_sigma2))) * MarsGMMParameters::soft_assoc_c3 * singleAssoc.m_normal_weight * singleAssoc.m_view_dir_weight;
        dde(errorIdx) = -0.5 * diff.dot(invcov_ps_diff);
        nde(errorIdx) =  -0.5 * std::pow(singleAssoc.m_model_normal.template cast<Scalar>().dot(invcov_ps_diff),2) * inv_normal_sig;
    }
    wes *= dde.exp() * nf * nde.exp();
    sum_weight += wes.sum();
    if (sum_weight > 0.0)
    {
        const Mat3 curRot_d_point_d_so3 = (m_currentRelativeRotation * (-Sophus::SO3<NormalEquationScalarType>::hat(mean_scene.template cast<NormalEquationScalarType>()))).template cast<Scalar>();
        Mat3K Dx = m_pDx;
        for ( int i = 0; i < N; ++i ) {
            Dx.template block<3,3>(0,D*i+3).noalias() = curRot_d_point_d_so3.template cast<NormalEquationScalarType>() * m_pDx.template block<3,3>(0,D*i+3);
        }
        const Scalar invSumWeight = 1.0 / sum_weight;
        const Scalar invWeight = //assoc.m_scene_num_points *
                invSumWeight;
        wes *= invWeight;
        (*m_accumW) += wes.sum();
        (*m_accumE) += .5*(wes*err).sum();
        for ( int idx = errorStartIdx; idx < errorEndIdx; ++idx  )
        {
            const InfoType & error = (*m_errors)[idx];
            const Scalar weight = wes(idx);
            updateSingleNormalEquations2<typename InfoType::Scalar,N>( error.m_W*weight, error.m_diff*weight, Dx, *m_accumH, *m_accumB );
            ++(*m_accumN);
        }
    }
}


template<typename ModelAssociationT, int N>
inline void MarsInfoWeightor<ModelAssociationT,N>::compute_info( const SceneAssociationData& assoc ) const
{
    typedef float Scalar;
    typedef Eigen::Matrix<Scalar,3,3> Mat3;
    typedef Eigen::Matrix<Scalar,3,1> Vec3;

    const auto & mean_scene = assoc.m_scene_mean; // already in relative pose
    const Mat3 cov_scene = assoc.m_scene_cov.template cast<Scalar>() + Mat3(Vec3::Constant(0.0001 * mean_scene.norm()).asDiagonal()); // cov_scene.diagonal() += Eigen::Vector3d::Constant(0.0001 * mean_scene.norm());
    const Mat3 cov_scene_transformed = m_currentRelativeRotation.cast<Scalar>() * cov_scene * m_currentRelativeRotationT.cast<Scalar>();
    const Vec3 scene_mean_transformed = m_currentTransform.template cast<Scalar>() * mean_scene;

    // normalize weights
    Scalar sum_weight = m_prior_prob_norm_factor / sqrt( m_soft_assoc2_PI_3 * symDetAdd(cov_scene_transformed,Scalar(assoc.m_sigma2))); // determinant is equivalent to (cov_scene + diag(sigma2)) since rotation has det 1.

    const Scalar normal_sig = assoc.m_sigma2;
    const Scalar inv_normal_sig = 1./normal_sig;
    const Scalar nf = 1. / sqrt(2*M_PI*normal_sig);

    const int start_idx = assoc.m_indices(0);
    const int end_idx = assoc.m_indices(1);
    const int errorStartIdx = 0;
    const int errorEndIdx = end_idx - start_idx;
    Eigen::Array<Scalar,Eigen::Dynamic,1> err (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> wes (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> dde (errorEndIdx,1);
    Eigen::Array<Scalar,Eigen::Dynamic,1> nde (errorEndIdx,1);
    for ( int idx = start_idx, errorIdx = errorStartIdx; idx < end_idx; ++idx, ++errorIdx  )
    {
        const ModelAssociationT & singleAssoc = (*m_model_associations)[idx];
        InfoType & error = (*m_errors)[errorIdx];
        const Mat3 cov_model = singleAssoc.m_model_cov.template cast<Scalar>();
        const Vec3 diff = scene_mean_transformed - singleAssoc.m_model_mean.template cast<Scalar>();
        const Mat3 cov = cov_model + cov_scene_transformed;
        error.m_W.noalias() = symInverse<Scalar>(cov).template cast<typename InfoType::Scalar>();

        error.m_diff.noalias() = error.m_W*diff.template cast<typename InfoType::Scalar>();
        err(errorIdx) = error.m_diff.dot(diff.template cast<typename InfoType::Scalar>());

        const Vec3 invcov_ps_diff = symInverseMultiplyDiagAdd(cov, Scalar(assoc.m_sigma2), diff);
        wes(errorIdx) = //singleAssoc.m_model_num_points *
                MarsGMMParameters::soft_assoc_c1 / sqrt(m_soft_assoc2_PI_3 * symDetAdd(cov,Scalar(assoc.m_sigma2))) * MarsGMMParameters::soft_assoc_c3 * singleAssoc.m_normal_weight * singleAssoc.m_view_dir_weight;
        dde(errorIdx) = -0.5 * diff.dot(invcov_ps_diff);
        nde(errorIdx) =  -0.5 * std::pow(singleAssoc.m_model_normal.template cast<Scalar>().dot(invcov_ps_diff),2) * inv_normal_sig;
    }
    wes *= dde.exp() * nf * nde.exp();
    sum_weight += wes.sum();
    if (sum_weight > 0.0)
    {
        const Scalar invSumWeight = 1.0 / sum_weight;
        const Scalar invWeight = //assoc.m_scene_num_points *
                invSumWeight;
        wes *= invWeight;

        Mat36 curT_d_point_d_se3;
        curT_d_point_d_se3.leftCols<3>() = m_curT_d_point_d_t;
        curT_d_point_d_se3.rightCols<3>() = m_currentRelativeRotation * (-Sophus::SO3<NormalEquationScalarType>::hat(mean_scene.template cast<NormalEquationScalarType>()));

        for ( int idx = errorStartIdx; idx < errorEndIdx; ++idx  )
        {
            const InfoType & error = (*m_errors)[idx];
            const Scalar weight = wes(idx);
            updateSingleNormalEquations3<typename InfoType::Scalar>( error.m_W*weight, curT_d_point_d_se3, *m_accumC );
        }
    }
}

// further possibility for exp optimization:
// https://stackoverflow.com/questions/412019/math-optimization-in-c-sharp/412102#412102
// https://stackoverflow.com/questions/10552280/fast-exp-calculation-possible-to-improve-accuracy-without-losing-too-much-perfo

template<typename ModelAssociationT, int N>
void MarsInfoWeightor<ModelAssociationT,N>::for_each( )
{
    if ( m_scene_associations->empty() ) return;
    *m_accumW = 0;
    *m_accumE = 0;
    std::vector<InfoType> tmp_errors ( 27 );
    m_errors = &tmp_errors;
    ZoneScopedN("MarsCovWeight::weight_all");
    const size_t num_assocs = m_scene_associations->size();
    for ( size_t idx = 0; idx < num_assocs; ++idx )
    {
        const SceneAssociationData & assoc = (*m_scene_associations)[idx];
        (*this)(assoc);
    }
}

template<typename ModelAssociationT, int N>
void MarsInfoWeightor<ModelAssociationT,N>::info( )
{
    if ( m_scene_associations->empty() ) return;
    m_accumC->setZero();
    std::vector<InfoType> tmp_errors ( 27 );
    m_errors = &tmp_errors;
    ZoneScopedN("MarsCovWeight::cov_all");
    const size_t num_assocs = m_scene_associations->size();
    for ( size_t idx = 0; idx < num_assocs; ++idx )
    {
        const SceneAssociationData & assoc = (*m_scene_associations)[idx];
        compute_info(assoc);
    }
}

// non vectorized
template class MarsInfoWeightor<ModelAssociationData, MarsSplineType::N>;
template class MarsInfoWeightor<ModelAssociationData, 1>;


