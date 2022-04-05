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
#include "MarsAssociator.h"
#include "MarsAssociations.h"
#include "MarsMap.h"
#include <loguru.hpp>
#include <Tracy.hpp>

template<typename ModelAssociationT>
MarsAssociator<ModelAssociationT>::MarsAssociator( MarsMap * model, const SurfelInfoConstPtrVector* scene_cells,
                                std::vector<SceneAssociationData> * scene_associations,
                                std::vector<ModelAssociationT> * model_associations,
                                int * num_valid_scene,
                                int * num_valid_model,
                                const Sophus::SE3d & transform, const int & neighbors )
    : m_model ( model ), m_scene_cells ( scene_cells ), m_scene_associations ( scene_associations ), m_model_associations( model_associations ), m_neighbors ( neighbors ), m_transform ( transform ),
      m_num_valid_scene ( num_valid_scene ), m_num_valid_model ( num_valid_model )
{ }
template<typename ModelAssociationT>
void MarsAssociator<ModelAssociationT>::for_each( )
{
    ZoneScopedN("MarsAssociator::associate_all");
    *m_num_valid_scene = 0, *m_num_valid_model = 0;
    m_scene_associations->resize(m_scene_cells->size());
    m_model_associations->resize(m_scene_cells->size() * 27);
    for ( const SurfelInfoConstPtr & cell : *m_scene_cells )
        (*this)(cell);
    //LOG(1) << "associated: scene: " << m_scene_associations->size() << " valid: "  << (*m_num_valid_scene) << " model: " << m_model_associations->size() << " valid: " << (*m_num_valid_model);
    m_scene_associations->erase(m_scene_associations->begin()+*m_num_valid_scene, m_scene_associations->end());
    m_model_associations->erase(m_model_associations->begin()+*m_num_valid_model, m_model_associations->end());
}

// https://developer.download.nvidia.com/cg/acos.html
// faster approx acos:
// Absolute error <= 6.7e-5
template <typename ModelAssociationT> template<typename Scalar>
inline Scalar MarsAssociator<ModelAssociationT>::faster_acos( const Scalar & x_ ) const {
  const Scalar negate = Scalar(x_ < 0);
  const Scalar x = abs(x_);
  const Scalar ret = (((Scalar(-0.0187293) * x + Scalar(0.0742610)) * x - Scalar(0.2121144)) * x + Scalar(1.5707288)) * (-x + Scalar(1)).sqrt();
  return negate * M_PI + (ret - 2 * negate * ret);
}
template <typename ModelAssociationT> template<typename Scalar>
inline Eigen::Array<Scalar,2,1> MarsAssociator<ModelAssociationT>::faster_vec_acos( const Eigen::Array<Scalar,2,1> & x_ ) const {
  const Eigen::Array<Scalar,2,1> negate = (x_ < 0).template cast<Scalar>();
  const Eigen::Array<Scalar,2,1> x = x_.abs();
  const Eigen::Array<Scalar,2,1> ret = (((x * Scalar(-0.0187293)  + Scalar(0.0742610)) * x - Scalar(0.2121144)) * x + Scalar(1.5707288)) * (-x + Scalar(1)).sqrt();
  return negate * M_PI + (ret - 2 * negate * ret);
}

template<typename ModelAssociationT>
inline void MarsAssociator<ModelAssociationT>::operator()( const SurfelInfoConstPtr & sceneCell ) const
{
    //ZoneScopedN("association_one");
    const Surfel & sceneSurfel = *sceneCell->m_surfel;
    const Eigen::Vector3f pos = sceneCell->m_center_s + sceneSurfel.mean_;

    std::vector<SurfelInfoConstPtr> cellPtrs;
    m_model->getSensorCell( m_transform.template cast<float>() * pos, sceneCell->m_level, cellPtrs, m_neighbors );

    if ( cellPtrs.empty() ) { return; }
    const size_t num_cell_ptrs = cellPtrs.size();

    SceneAssociationData & assoc = (*m_scene_associations)[(*m_num_valid_scene)];
    assoc.m_cell_scene = sceneCell;
    assoc.m_scene_mean = pos.cast<typename ModelAssociationT::Scalar>();
    assoc.m_indices(0) = *m_num_valid_model;
    assoc.m_sigma2 = std::pow(MarsGMMParameters::sigma_size_factor * m_model->getCellSize(sceneCell->m_level),2);
    assoc.m_scene_num_points = sceneSurfel.getNumPoints();
    assoc.m_scene_cov = sceneSurfel.cov_.cast<typename ModelAssociationT::Scalar>();

    typedef float Scalar;
    const Eigen::Vector3f normal_scene_transformed = m_transform.so3().cast<float>() * sceneSurfel.normal_; // rotate in correct direction;
    const Eigen::Vector3f view_dir_scene_transformed = m_transform.so3().cast<float>() * sceneSurfel.first_view_dir_;

    typedef Eigen::Array<Scalar,2,1> ArrayType;
    Eigen::VectorXi model_points = Eigen::VectorXi::Zero(num_cell_ptrs,1);

    constexpr Scalar min_exponent = -2.0;
    constexpr Scalar min_weight = -4.5;
    constexpr Scalar min_weight_add = (min_weight+5.0);
    constexpr Scalar additional_weight = 5.0;
    constexpr Scalar normal_std_mod =  0.125*M_PI;
    constexpr Scalar minus_inv_half_normal_std_mod_2 = -.5 / ( normal_std_mod*normal_std_mod );
    // iterate through neighbors of the directly associated node to eventually find a better match
    //ZoneScopedN("associate_cells");
    for ( size_t idx = 0; idx < num_cell_ptrs; ++idx)
    {
        const SurfelInfoConstPtr & cell = cellPtrs[idx];
        const Surfel & modelSurfel = *(cell->m_surfel);

        ModelAssociationT & singleAssoc = (*m_model_associations)[*m_num_valid_model];
        singleAssoc.m_model_mean.noalias() = cell->m_center_s.template cast<typename ModelAssociationT::Scalar>() + modelSurfel.mean_.template cast<typename ModelAssociationT::Scalar>();
        singleAssoc.m_model_normal.noalias() = modelSurfel.normal_.cast<typename ModelAssociationT::Scalar>();
        singleAssoc.m_model_cov.noalias() = modelSurfel.cov_.cast<typename ModelAssociationT::Scalar>();
        singleAssoc.m_model_num_points = modelSurfel.getNumPoints();

        const ArrayType dp ( normal_scene_transformed.dot(singleAssoc.m_model_normal.template cast<Scalar>()), view_dir_scene_transformed.dot(modelSurfel.first_view_dir_.template cast<Scalar>()));
        const ArrayType dd = (dp.min(Scalar(1))).max(Scalar(-1));
        const ArrayType expo = faster_vec_acos(dd).square() * minus_inv_half_normal_std_mod_2;
        singleAssoc.m_normal_weight = expo(0) < min_exponent ? min_weight_add : expo(0) + additional_weight;
        singleAssoc.m_view_dir_weight = expo(1) < min_exponent ? min_weight_add : expo(1) + additional_weight;
        model_points(idx) = singleAssoc.m_model_num_points;
        //assoc.m_model_num_points += singleAssoc.m_model_num_points;
        ++(*m_num_valid_model);
    }
    assoc.m_model_num_points += model_points.sum();
    assoc.m_indices(1) = *m_num_valid_model;
    if ( (*m_num_valid_model) > assoc.m_indices(0) )
        ++(*m_num_valid_scene);
}

template class MarsAssociator<ModelAssociationData>;
