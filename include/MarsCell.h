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
#include "MarsSurfelInfo.h"
#include "MarsFwdDec.h"

#include <loguru.hpp>
#include <deque>

struct Cell
{
    typedef typename MapParameters::IndexType IndexType;
    typedef typename MapParameters::CellIndexType CellIndexType;
    typedef SurfelInfo::SurfelPtr SurfelPtr;
    typedef SurfelInfo::SemanticClassType SemanticClassType;
    typedef SurfelInfo::SemanticClassType::Ptr SemanticClassPtr;
    typedef std::deque<std::pair<IndexType,SurfelPtr>> CovContainerType; // per scan one Vector

    CovContainerType m_cov_per_scan;

    SurfelPtr m_surfel = Surfel::create();
    SurfelPtr cur_scan_surfel = nullptr;
    SemanticClassPtr m_class = SemanticClassType::create();
    bool m_updated = true;

    SurfelInfo m_surfel_info;
    SurfelInfoConstPtr getSurfelInfo ( ) const
    {
        return &m_surfel_info;
    }

    void create( const Eigen::Vector3f & center_s, const Eigen::Vector3f & first_view_dir, const bool & in_center, const MapParameters::IndexType & lvl, const MapParameters::CellIndexType & idx );

    inline void addNewScan( const int & scan_id, const Eigen::Vector3f & view_dir )
    {
        m_updated = false;
        cur_scan_surfel = Surfel::create();
        cur_scan_surfel->first_view_dir_ = view_dir.cast<Surfel::Scalar>();
        m_cov_per_scan.emplace_back(scan_id,cur_scan_surfel);
    }

    inline void addPoint( const Eigen::Vector3f & pt )
    {
        if ( ! cur_scan_surfel ) LOG(FATAL) << "why is this surfel empty?";
        m_updated = false;
        cur_scan_surfel->add(pt-m_surfel_info.m_center_s);
    }

    inline void addSemantic( const Eigen::VectorXf & semantic )
    {
        m_updated = false;
        m_class->addProb(semantic);
    }

    bool updateSurfel( const bool & full = false);

    inline bool isUpdated() const { return m_updated; }
    inline bool isEmpty() const
    {
        return m_surfel == nullptr || m_cov_per_scan.empty() || ( isUpdated() && m_surfel->getNumPoints() == 0);
    }

    inline int getNumPoints() const
    {
        return (m_surfel == nullptr || ! isUpdated()) ? 0 : m_surfel->getNumPoints();
    }


    Eigen::Vector3f getFirstViewDir ( ) const
    {
        if ( cur_scan_surfel )
            return cur_scan_surfel->first_view_dir_;
        else
            return Eigen::Vector3f::Zero();
    }

    void add ( const Cell & other );
    void addLowerCell ( const Cell & other );

    bool removeByScanId ( const IndexType & scan_id );

    void addSurfel( const SurfelPtr & surfel, const Sophus::SE3d & pose_s1s2, const IndexType & scan_id, const SemanticClassPtr & class_ = nullptr);
};
