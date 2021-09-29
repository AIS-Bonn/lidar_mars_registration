/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration

Copyright (c) 2016-2021, Computer Science Institute VI, University of Bonn.
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

template<typename ModelAssociationT, typename Scalar>
class MarsEvaluator
{
public:

    MarsEvaluator( const std::vector<SceneAssociationData> * scene_associations,
                   const std::vector<ModelAssociationT> * model_associations,
                   NormalEquationScalarType * accum_W, NormalEquationScalarType * accum_E,
                   const Sophus::SE3d & currTransform,
                   const size_t & num_points_model, const size_t & num_points_scene );
    ~MarsEvaluator(){}

    void for_each( );

    Sophus::SE3d m_currentTransform;
    Eigen::Matrix<Scalar,3,3> m_currentRelativeRotation;
    Eigen::Matrix<Scalar,3,3> m_currentRelativeRotationT;

    Scalar m_soft_assoc2_PI_3 = std::numeric_limits<Scalar>::signaling_NaN();
    Scalar m_prior_prob_norm_factor = std::numeric_limits<Scalar>::signaling_NaN();

    NormalEquationScalarType * m_accumW = nullptr;
    NormalEquationScalarType * m_accumE = nullptr;
    const std::vector<SceneAssociationData>* m_scene_associations = nullptr;
    const std::vector<ModelAssociationT>* m_model_associations = nullptr;
    inline void operator()( const SceneAssociationData& assoc ) const;
};


template<typename ModelAssociationT, typename Scalar>
MarsEvaluator<ModelAssociationT,Scalar>::MarsEvaluator( const std::vector<SceneAssociationData> * scene_associations,
                                                 const std::vector<ModelAssociationT> * model_associations,
                                                 NormalEquationScalarType * accumW, NormalEquationScalarType * accumE,
                                                 const Sophus::SE3d & currTransform,
                                                 const size_t & num_points_model, const size_t & num_points_scene )
    : m_scene_associations ( scene_associations ), m_model_associations ( model_associations ), m_accumW ( accumW ), m_accumE ( accumE ), m_currentTransform ( currTransform )
{
    m_currentRelativeRotation = m_currentTransform.so3().template cast<NormalEquationScalarType>().matrix();
    m_currentRelativeRotationT = m_currentRelativeRotation.transpose();

    m_soft_assoc2_PI_3 = MarsGMMParameters::soft_assoc_c2 * std::pow(M_PI,3);
    m_prior_prob_norm_factor = MarsGMMParameters::prior_prob / (1.0 - MarsGMMParameters::prior_prob) * (num_points_model / (NormalEquationScalarType) num_points_scene) * MarsGMMParameters::soft_assoc_c1;
}

template<typename ModelAssociationT, typename Scalar>
inline void MarsEvaluator<ModelAssociationT,Scalar>::operator()( const SceneAssociationData& assoc ) const
{
    typedef Eigen::Matrix<Scalar,3,3> Mat3;
    typedef Eigen::Matrix<Scalar,3,1> Vec3;

    const auto & mean_scene = assoc.m_scene_mean; // already in relative pose
    const Mat3 cov_scene = assoc.m_scene_cov.template cast<Scalar>() + Mat3(Vec3::Constant(0.0001 * mean_scene.norm()).asDiagonal()); // cov_scene.diagonal() += Eigen::Vector3d::Constant(0.0001 * mean_scene.norm());
    const Mat3 cov_scene_transformed = m_currentRelativeRotation * cov_scene * m_currentRelativeRotationT;
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
        const Mat3 cov_model = singleAssoc.m_model_cov.template cast<Scalar>();
        const Vec3 diff = scene_mean_transformed - singleAssoc.m_model_mean.template cast<Scalar>();
        const Mat3 cov = cov_model + cov_scene_transformed;
        err(errorIdx) = diff.dot( symInverseMultiply(cov,diff) );

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
        (*m_accumW) += wes.sum();
        (*m_accumE) += .5*(wes*err).sum();
    }
}
// further possibility for exp optimization:
// https://stackoverflow.com/questions/412019/math-optimization-in-c-sharp/412102#412102
// https://stackoverflow.com/questions/10552280/fast-exp-calculation-possible-to-improve-accuracy-without-losing-too-much-perfo

template<typename ModelAssociationT, typename Scalar>
void MarsEvaluator<ModelAssociationT,Scalar>::for_each( )
{
    if ( m_scene_associations->empty() ) return;
    *m_accumW = 0;
    *m_accumE = 0;
    ZoneScopedN("MarsCovWeight::weight_all");
    const size_t num_assocs = m_scene_associations->size();
    for ( size_t idx = 0; idx < num_assocs; ++idx )
    {
        const SceneAssociationData & assoc = (*m_scene_associations)[idx];
        (*this)(assoc);
    }
}

// non vectorized
template class MarsEvaluator<ModelAssociationData,NormalEquationScalarType>;
