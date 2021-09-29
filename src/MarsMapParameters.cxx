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
#include "MarsMapParameters.h"
#include "loguru.hpp"

double MapParameters::plane_scale_factor = 0.1;
double MapParameters::first_ev_planarity_threshold = 0.1;
double MapParameters::second_ev_degenerate_threshold = 0.1;

Eigen::Vector3i MapParameters::getMiddleIndexVector( const int & num_cells )
{
    return Eigen::Vector3i( num_cells>>1, num_cells>>1, num_cells>>1 );
}

int MapParameters::getMiddleIndex( const int & num_cells )
{
    return getMiddleIndexVector(num_cells).dot(Eigen::Vector3i(1, num_cells, num_cells * num_cells ));
}

const MapParameters MapParameters::getMapParameters(const MapParameters & params )
{
    return MapParameters ( params.m_size, params.m_num_cells, params.m_num_levels, params.m_min_origin_bounds, params.m_max_origin_bounds, params.m_use_adaptive, params.m_omit_center );
}

const MapParameters MapParameters::getMapParameters(const float & size, const IndexType & num_cells, const LevelIndexType & levels, const bool & use_adaptive, const bool & omit_center )
{
    const IndexType used_num_cells = ((num_cells>>1)<<1);
    //const IndexType used_num_cells = (num_cells|0x1);// TODO: disallow this, but fix the other stuff first!
    Eigen::VectorXf min_origin_bounds = Eigen::VectorXf::Zero(levels,1);
    Eigen::VectorXf max_origin_bounds = Eigen::VectorXf::Zero(levels,1);
    for ( LevelIndexType lvl = 0; lvl < levels; ++lvl )
    {
        const float halfSideLength = getMapSideLength(size, lvl) / 2;
        max_origin_bounds(lvl) = halfSideLength;
        min_origin_bounds(lvl) = -halfSideLength;
    }
    return MapParameters ( size, used_num_cells, levels, min_origin_bounds, max_origin_bounds, use_adaptive, omit_center );
}

// Block Sparse
// CHECK
template<typename Scalar>
Eigen::Vector3i MapParameters::toBlockIndex ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    return (point * getBlockScale(level)).array().floor().template cast<IndexType>(); // ((numcells>>1)<<1) // even
}
template Eigen::Vector3i MapParameters::toBlockIndex<double>(const Eigen::Matrix<double,3,1> &, const LevelIndexType & ) const;
template Eigen::Vector3i MapParameters::toBlockIndex<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;

template<typename Scalar>
void MapParameters::toBlockIndexAndIndexVector ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level, Eigen::Vector3i & block_idx_vec, Eigen::Vector3i & block_idx_inside_vec ) const
{
    const Eigen::Vector3i cell_idx_vec = toCellIndexVector ( point, level );
    block_idx_vec = cell_idx_vec / block_size;
    block_idx_inside_vec = cell_idx_vec - block_idx_vec * block_size;
}
template void MapParameters::toBlockIndexAndIndexVector<double>(const Eigen::Matrix<double,3,1> &, const LevelIndexType &, Eigen::Vector3i &, Eigen::Vector3i & ) const;
template void MapParameters::toBlockIndexAndIndexVector<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType &, Eigen::Vector3i &, Eigen::Vector3i & ) const;


void MapParameters::fromCellIndexVectorToBlockIndices  ( const Eigen::Vector3i & cell_idx_vec, IndexType & blockIdx, IndexType & blockInsideIdx ) const
{
    const Eigen::Vector3i block_idx_vec = fromCellIndexVectorToBlockIdx( cell_idx_vec );
    const Eigen::Vector3i block_idx_inside_vec = fromCellIndexVectorToBlockInsideIdx( cell_idx_vec, block_idx_vec );

    blockIdx = fromBlockIndexVectorToBlockIdx ( block_idx_vec );
    blockInsideIdx = fromBlockIndexInsideVectorToBlockIdx ( block_idx_inside_vec );
}

MapParameters::IndexType MapParameters::fromBlockIndexInsideVectorToBlockIdx ( const Eigen::Vector3i & block_idx_inside_vec ) const
{
    const Eigen::Vector3i num_blocks_per_side ( 1, block_size, block_size * block_size );
    return block_idx_inside_vec.dot(num_blocks_per_side);
}

MapParameters::IndexType MapParameters::fromBlockIndexVectorToBlockIdx ( const Eigen::Vector3i & block_idx_vec ) const
{
    const Eigen::Vector3i p = block_idx_vec + blockMiddleIndexVector;
    const Eigen::Vector3i num_cells_per_side ( 1, m_num_cells/block_size, m_num_cells/block_size * m_num_cells/block_size );
    return p.dot(num_cells_per_side);
}

template<typename Scalar>
Eigen::Matrix<Scalar,3,1> MapParameters::fromCellIdxToCenter ( const Eigen::Vector3i & index, const LevelIndexType & level ) const
{
    return (index.cast<Scalar>() * getCellSizeOnLevel(level)).array() + getCellSizeOnLevel(level) / 2; // even
    //return index.cast<double>() * getCellSizeOnLevel(level); // uneven
}
template Eigen::Matrix<double,3,1> MapParameters::fromCellIdxToCenter<double>(const Eigen::Vector3i &, const LevelIndexType & ) const;
template Eigen::Matrix<float,3,1> MapParameters::fromCellIdxToCenter<float>(const Eigen::Vector3i&, const LevelIndexType & ) const;


template<typename Scalar>
Eigen::Matrix<Scalar,3,1> MapParameters::fromCellIndexToCenter ( const IndexType & idx, const LevelIndexType & level ) const
{
    return (fromCellIndexToCellIdx ( idx ).template cast<Scalar>() * getCellSizeOnLevel(level)).array() + getCellSizeOnLevel(level) / 2; // even
    //return fromCellIndexToCellIdx ( idx ).cast<double>() * getCellSizeOnLevel(level); // uneven
}

template Eigen::Matrix<double,3,1> MapParameters::fromCellIndexToCenter<double>(const IndexType &, const LevelIndexType & ) const;
template Eigen::Matrix<float,3,1> MapParameters::fromCellIndexToCenter<float>(const IndexType &, const LevelIndexType & ) const;

Eigen::Vector3i MapParameters::fromCellIndexToCellIdx( const IndexType & idx) const
{
    IndexType pm = idx;
    Eigen::Vector3i center;
    center(0) =  pm % m_num_cells;
    pm /= m_num_cells;
    center(1) = pm % m_num_cells;
    pm /= m_num_cells;
    center(2) = pm;
    return center-middleIndexVector;
}

template<typename Scalar>
Eigen::Vector3i MapParameters::toCellIndexVector ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    //return (point * getCellScale(level)).array().round().cast<IndexType>(); // (numcells|0x1) // uneven
    return (point * getCellScale(level)).array().floor().template cast<IndexType>(); // ((numcells>>1)<<1) // even
}
template Eigen::Vector3i MapParameters::toCellIndexVector<double>(const Eigen::Matrix<double,3,1> &, const LevelIndexType & ) const;
template Eigen::Vector3i MapParameters::toCellIndexVector<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;

MapParameters::IndexType MapParameters::toCellIndexFromIndexVector ( const Eigen::Vector3i & idx_vec ) const
{
    const Eigen::Vector3i p = idx_vec + middleIndexVector;
    const Eigen::Vector3i num_cells_per_side ( 1, m_num_cells, m_num_cells * m_num_cells );
    return p.dot(num_cells_per_side);
}

template<typename Scalar>
MapParameters::IndexType MapParameters::toCellIndexOffset ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    const Eigen::Vector3i num_cells_per_side ( 1, m_num_cells, m_num_cells * m_num_cells );
    return toCellIndexVector(point, level).dot(num_cells_per_side);
}

template MapParameters::IndexType MapParameters::toCellIndexOffset<double>(const Eigen::Matrix<double,3,1>  &, const LevelIndexType & ) const;
template MapParameters::IndexType MapParameters::toCellIndexOffset<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;


template<typename Scalar>
MapParameters::IndexType MapParameters::toCellIndex ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    return toCellIndexFromIndexVector ( toCellIndexVector( point, level ) );
}
template MapParameters::IndexType MapParameters::toCellIndex<double>(const Eigen::Matrix<double,3,1>  &, const LevelIndexType & ) const;
template MapParameters::IndexType MapParameters::toCellIndex<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;


template<typename Scalar>
Eigen::VectorXi MapParameters::toCellIndices ( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pts, const LevelIndexType & level ) const
{
    const Eigen::Matrix3Xi m = (pts * getCellScale(level)).array().round().template cast<IndexType>();
    const Eigen::Matrix3Xi m2 = m.colwise() + middleIndexVector;
    const Eigen::Vector3i num_cells_per_side ( 1, m_num_cells, m_num_cells * m_num_cells );
    const Eigen::VectorXi cell_indices = (num_cells_per_side.transpose() * m2).transpose();
    return cell_indices;
}
template Eigen::VectorXi MapParameters::toCellIndices<double>(const Eigen::Matrix<double,3,Eigen::Dynamic>  &, const LevelIndexType & ) const;
template Eigen::VectorXi MapParameters::toCellIndices<float>(const Eigen::Matrix<float,3,Eigen::Dynamic> &, const LevelIndexType & ) const;

template<typename Scalar>
bool MapParameters::isInBounds ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    return level >= 0 && level < m_min_origin_bounds.rows()
            && (point.array() > m_min_origin_bounds(level)).all() && (point.array() < m_max_origin_bounds(level)).all();
}
template bool MapParameters::isInBounds<double>(const Eigen::Matrix<double,3,1> &, const LevelIndexType & ) const;
template bool MapParameters::isInBounds<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;

template<typename Scalar>
MapParameters::ArrX MapParameters::areInBounds( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pt_s ) const
{
    ArrX inBounds = ArrX::Constant(m_min_origin_bounds.rows(),pt_s.cols(),true);
    for ( int level = 0; level < m_min_origin_bounds.rows(); ++level)
        inBounds.row(level) = (pt_s.array() > m_min_origin_bounds(level)).colwise().all() && (pt_s.array() < m_max_origin_bounds(level)).colwise().all();
    return inBounds;
}
template MapParameters::ArrX MapParameters::areInBounds<double>(const Eigen::Matrix<double,3,Eigen::Dynamic> & ) const;
template MapParameters::ArrX MapParameters::areInBounds<float>(const Eigen::Matrix<float,3,Eigen::Dynamic> & ) const;

template<typename Scalar>
MapParameters::ArrX MapParameters::areInBoundsT( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pt_s ) const
{
    ArrX inBounds = ArrX::Constant(pt_s.cols(),m_min_origin_bounds.rows(),true);
    for ( int level = 0; level < m_min_origin_bounds.rows(); ++level)
        inBounds.col(level) = ((pt_s.array() > m_min_origin_bounds(level)).colwise().all() && (pt_s.array() < m_max_origin_bounds(level)).colwise().all());
    return inBounds;
}
template MapParameters::ArrX MapParameters::areInBoundsT<double>(const Eigen::Matrix<double,3,Eigen::Dynamic> & ) const;
template MapParameters::ArrX MapParameters::areInBoundsT<float>(const Eigen::Matrix<float,3,Eigen::Dynamic> &) const;

template<typename Scalar>
bool MapParameters::isInBoundsExceptCenter ( const Eigen::Matrix<Scalar,3,1> & point, const LevelIndexType & level ) const
{
    return ( level >= 0 && level+1 < m_min_origin_bounds.rows() && isInBounds(point, level) && !isInBounds(point,level+1) );
}
template bool MapParameters::isInBoundsExceptCenter<double>(const Eigen::Matrix<double,3,1> &, const LevelIndexType & ) const;
template bool MapParameters::isInBoundsExceptCenter<float>(const Eigen::Matrix<float,3,1> &, const LevelIndexType & ) const;

int MapParameters::getNumCells ( const LevelIndexType & level ) const
{
    return m_num_cells;
}

MapParameters::MapParameters( const float & size, const int & num_cells, const LevelIndexType & levels, const Eigen::VectorXf & min_bounds, const Eigen::VectorXf & max_bounds, const bool & use_adaptive, const bool & omit_center )
    : m_size(size)
    , m_num_cells(num_cells)
    , m_num_levels( std::max<LevelIndexType>(0,std::min<LevelIndexType>(levels,maxAllowedMapLevels)) )
    , middleIndexVector ( MapParameters::getMiddleIndexVector(num_cells) )
    , middleIndex ( MapParameters::getMiddleIndex(num_cells) )
    , m_min_origin_bounds ( min_bounds )
    , m_max_origin_bounds ( max_bounds )
    , m_use_adaptive ( use_adaptive )
    , m_omit_center ( omit_center )
{
    if ( num_cells > 1024) LOG(FATAL) << "This number of cells (" << num_cells << ") is not supported as it would overflow the int index arithmetik.";
}
