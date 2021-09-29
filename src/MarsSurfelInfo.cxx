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
#include "MarsSurfelInfo.h"
#include <loguru.hpp>

SurfelInfo::SurfelInfo ( SurfelPtr surfel, const Eigen::Vector3f & center_s, const LevelIndexType & level, const CellIndexType & idx, const IndexType & scan_id, const bool & in_center, const SemanticClassPtr & class_  )
    : m_level ( level ), m_index ( idx ), m_center_s(center_s), m_surfel(surfel), m_scan_id ( scan_id ), m_in_center ( in_center ), m_class ( class_ )
{}

void SurfelInfo::addClass ( const SemanticClassPtr & class_ )
{
    if ( ! class_ ) return;
    m_class = SemanticClassType::create( class_ );
}

void SurfelInfo::addSurfel( const SurfelPtr & surfel, const Sophus::SE3d & pose_s1s2 ) //const
{
    if ( ! surfel ) return;
    if ( ! m_surfel ) LOG(FATAL) << "hui...";
    Surfel::Ptr o_surfel = Surfel::create(surfel);
    o_surfel->transform(pose_s1s2);
    if ( ! o_surfel ) LOG(FATAL) << "why is o-surfel not there?";

    if ( m_surfel )
        (*m_surfel) += (*o_surfel);
    else
        m_surfel = o_surfel;
}

SurfelInfo SurfelInfo::clone() const
{
    SurfelInfo si;
    si.m_level = m_level;
    si.m_scan_id = m_scan_id;
    si.m_index = m_index;
    si.m_in_center = m_in_center;
    si.m_center_s = m_center_s;
    si.m_surfel = Surfel::create(m_surfel);
    si.m_class = SemanticClassType::create(m_class);
    return si;
}
