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
#include "MarsMapParameters.h"
#include "MarsSurfel.h"
#include "MarsSemanticClass.h"

struct SurfelInfo
{
    //static constexpr int NumSemanticClasses = 15;
    static constexpr int NumSemanticClasses = 20;

    typedef std::shared_ptr<SurfelInfo> Ptr;
    typedef Surfel::Ptr SurfelPtr;
    typedef typename MapParameters::IndexType IndexType;
    typedef typename MapParameters::LevelIndexType LevelIndexType;
    typedef typename MapParameters::CellIndexType CellIndexType;
    typedef SemanticClass<NumSemanticClasses> SemanticClassType;
    typedef SemanticClassType::Ptr SemanticClassPtr;

    LevelIndexType m_level = -1;
    IndexType m_scan_id = -1;
    CellIndexType m_index;
    bool m_in_center = false;
    Eigen::Vector3f m_center_s = Eigen::Vector3f::Constant(std::numeric_limits<float>::signaling_NaN());// center relative to sensor
    SurfelPtr m_surfel = nullptr;

    SurfelInfo clone () const;

    SemanticClassPtr m_class = nullptr;

    SurfelInfo ( SurfelPtr surfel = Surfel::create(), const Eigen::Vector3f & center_s = Eigen::Vector3f::Constant(std::numeric_limits<float>::signaling_NaN()), const LevelIndexType & level = -1, const CellIndexType & idx = CellIndexType(3,0), const IndexType & scan_id = 0, const bool & in_center = false, const SemanticClassPtr & class_ = nullptr );

    void addClass ( const SemanticClassPtr & class_ );
    void addSurfel( const SurfelPtr & surfel, const Sophus::SE3d & pose_s1s2 );
};
typedef const SurfelInfo * SurfelInfoConstPtr;
typedef std::vector<SurfelInfo> SurfelInfoVector;
typedef std::vector<SurfelInfoConstPtr> SurfelInfoConstPtrVector;
