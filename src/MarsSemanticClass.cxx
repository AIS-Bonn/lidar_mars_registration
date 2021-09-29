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
#include "MarsSemanticClass.h"

template <int _N>
typename SemanticClass<_N>::Ptr SemanticClass<_N>::create ( const SemanticClass<_N>::Ptr & other )
{
    if ( ! other )
        return std::make_shared<SemanticClass<_N>>( );

    SemanticClass<_N>::Ptr n = std::make_shared<SemanticClass<_N>>();
    n->m_empty = other->m_empty;
    n->m_class = other->m_class;
    n->m_use_log_prob = other->m_use_log_prob;
    return n;
}

template <int _N>
void SemanticClass<_N>::addLogProb ( const SemanticClass<_N>::VecN & new_log_prob )
{
    m_empty = false;
    m_use_log_prob = true;
    m_class += new_log_prob;
}

template <int _N>
void SemanticClass<_N>::addProb ( const SemanticClass<_N>::VecN & new_prob )
{
    m_empty = false;
    m_use_log_prob = false;
    m_class.array() += new_prob.array();
}

template <int _N>
void SemanticClass<_N>::add ( const SemanticClass<_N> & other )
{
    m_empty = false;
    m_class += other.m_class;
}

template <int _N>
int SemanticClass<_N>::getArgMaxClass ( ) const
{
    typename SemanticClass<_N>::VecN::Index idx = 0; // assuming background is first class.
    if ( m_empty ) return idx;
    if ( m_use_log_prob )
        m_class.array().exp().maxCoeff(&idx);
    else
        m_class.array().maxCoeff(&idx);
    return idx;
}

template class SemanticClass<15>;
template class SemanticClass<20>;
