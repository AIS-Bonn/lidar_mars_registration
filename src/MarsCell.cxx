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
#include "MarsCell.h"
#include "MarsSurfelInfo.h"
#include "MarsFwdDec.h"
#include "loguru.hpp"

void Cell::create( const Eigen::Vector3f & center_s, const Eigen::Vector3f & first_view_dir, const bool & in_center, const MapParameters::IndexType & lvl, const MapParameters::CellIndexType & idx )
{
    m_surfel = Surfel::create(); m_surfel->first_view_dir_ = first_view_dir.cast<Surfel::Scalar>();
    m_class = SemanticClassType::create();
    m_updated = true;

    m_surfel_info.m_surfel = m_surfel;
    m_surfel_info.m_in_center = in_center;
    m_surfel_info.m_center_s = center_s;
    m_surfel_info.m_level = lvl;
    m_surfel_info.m_scan_id = ( m_cov_per_scan.empty() ? -1 : m_cov_per_scan.front().first );
    for ( size_t i = 0; i < idx.size(); ++i )
        m_surfel_info.m_index[i] = idx[i];
    m_surfel_info.m_class = m_class;
}

void Cell::addSurfel( const SurfelPtr & surfel, const Sophus::SE3d & pose_s1s2, const IndexType & scan_id, const SemanticClassPtr & class_ )
{
    if ( ! surfel ) return;
    m_updated = false;
    SurfelPtr o_surfel = Surfel::create(surfel);
    if ( ! pose_s1s2.params().isApprox(Sophus::SE3d().params()))
        o_surfel->transform(pose_s1s2);
    m_cov_per_scan.emplace_back(scan_id, o_surfel);
    if ( m_surfel )
        *m_surfel += *o_surfel;
    else
    {
        m_surfel = Surfel::create(o_surfel);
        m_surfel_info.m_surfel = m_surfel;
    }
    if ( class_ )
    {
        m_class->add(*class_);
    }
}

bool Cell::updateSurfel( const bool & full )
{
    if ( m_surfel == nullptr ) LOG(FATAL) << "ooh boy!";
    if ( m_surfel == cur_scan_surfel ) LOG(FATAL) << "this should not be happing...";
    if ( full )
    {
        m_surfel->clear();
        if ( m_cov_per_scan.empty() ) { return false; }
        for ( const auto & id_cov : m_cov_per_scan )
        {
            *m_surfel += *id_cov.second;
        }
        m_surfel->evaluate();
    }
    else
    {
        if ( m_surfel )
        {
            if ( cur_scan_surfel )
                *m_surfel += *cur_scan_surfel;
            m_surfel->evaluate();
        }
        else
            return false;
    }
    m_updated = true;
    return true;
}

bool Cell::removeByScanId ( const IndexType & scan_id )
{
    if ( m_cov_per_scan.empty() || m_cov_per_scan.front().first != scan_id ) return false;
    m_cov_per_scan.pop_front();
    updateSurfel(true);
    return true;
}

void Cell::add ( const Cell & other )
{
    m_updated = false;
    m_class->add(*other.m_class);

    for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : other.m_cov_per_scan )
    {
        const IndexType & scan_id = surfelInfo.first;
        const SurfelPtr & surfel = surfelInfo.second;


        bool found = false;
        for ( const std::pair<IndexType,SurfelPtr> & thisSurfelInfo: m_cov_per_scan )
        {
            const IndexType & thisScanId = thisSurfelInfo.first;
            if ( thisScanId == scan_id ) { found = true; *(thisSurfelInfo.second) += *surfel; break; }
        }
        if ( ! found )
        {
            m_cov_per_scan.emplace_back(scan_id,surfel);
        }
    }
}

void Cell::addLowerCell ( const Cell & other )
{
    m_class->add(*other.m_class);

    const Surfel::Vec3 t = m_surfel_info.m_center_s - other.m_surfel_info.m_center_s;
    for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : other.m_cov_per_scan )
    {
        const IndexType & scan_id = surfelInfo.first;
        const SurfelPtr & surfel = surfelInfo.second;
        bool found = false;
        for ( const std::pair<IndexType,SurfelPtr> & thisSurfelInfo: m_cov_per_scan )
        {
            const IndexType & thisScanId = thisSurfelInfo.first;
            if ( thisScanId == scan_id ) { found = true; *(thisSurfelInfo.second) += *surfel; break; }
        }
        if ( ! found )
        {
            m_cov_per_scan.emplace_back(scan_id,surfel);
            m_cov_per_scan.back().second->translate(t);
        }
    }
}
