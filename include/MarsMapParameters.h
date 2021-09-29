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
#include <memory>
#include <Eigen/Dense>
#include <absl/container/inlined_vector.h>

struct MapParameters
{
    typedef std::shared_ptr<const MapParameters> Ptr;
    typedef int LevelIndexType;
    typedef int IndexType;
    typedef absl::InlinedVector<IndexType,3> CellIndexType;

    static constexpr size_t maxAllowedMapLevels = 4;
    static constexpr int NUM_CELLS = ((64>>1)<<1);
    static constexpr float MAP_SIDE = 256;
    static constexpr int NUM_LEVELS = 3;

    static constexpr int block_size = 4;
    static constexpr int block_volume = block_size * block_size * block_size;


    static double plane_scale_factor;
    static double first_ev_planarity_threshold;
    static double second_ev_degenerate_threshold;

    bool m_use_adaptive = true;
    float m_size = MAP_SIDE; // map size in m
    int m_num_cells = NUM_CELLS;
    bool m_omit_center = true;

    LevelIndexType m_num_levels = 1;
    Eigen::Vector3i middleIndexVector = getMiddleIndexVector( NUM_CELLS );
    Eigen::Vector3i blockMiddleIndexVector = getMiddleIndexVector( NUM_CELLS/block_size );
    Eigen::Vector3i blockMiddleInsideIndexVector = getMiddleIndexVector( block_size );
    int middleIndex = getMiddleIndex( NUM_CELLS );
    int blockMiddleIndex = getMiddleIndex( NUM_CELLS/block_size );
    int blockInsideMiddleIndex = getMiddleIndex( block_size );

    Eigen::VectorXf m_min_origin_bounds;
    Eigen::VectorXf m_max_origin_bounds;

    typedef Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> ArrX;
    template <typename Scalar>
    ArrX areInBounds( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pt_s ) const;
    template <typename Scalar>
    ArrX areInBoundsT( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pt_s ) const;

    int m_scan_window_size_{25};

    constexpr float getCellSizeOnLevel ( const int & level ) const
    {
        return m_size / ((0x1<<level) * m_num_cells);
    }
    constexpr float getCellScale ( const int & level ) const
    {
        return ((0x1<<level) * m_num_cells) / m_size;
    }

    constexpr float getBlockScale ( const int & level ) const
    {
        return ((0x1<<level) * m_num_cells) / m_size / block_size;
    }

    constexpr float getMapSideLength( const int & level ) const
    {
        return getMapSideLength( m_size, level );
    }
    constexpr static float getMapSideLength( const float & size, const int & level )
    {
        return size / (0x1 << level);
    }

    int getNumCells ( const int & level ) const;

    static Eigen::Vector3i getMiddleIndexVector( const int & num_cells );

    static int getMiddleIndex( const int & num_cells );

    template<typename Scalar>
    Eigen::Matrix<Scalar,3,1> fromCellIndexToCenter ( const IndexType & idx, const int & level ) const;
    template<typename Scalar>
    Eigen::Matrix<Scalar,3,1> fromCellIdxToCenter ( const Eigen::Vector3i & index, const int & level ) const;

    Eigen::Vector3i fromCellIndexToCellIdx( const IndexType & idx) const;

    template<typename Scalar>
    Eigen::Vector3i toCellIndexVector ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;
    IndexType toCellIndexFromIndexVector ( const Eigen::Vector3i & idx_vec ) const;
    template<typename Scalar>
    IndexType toCellIndexOffset ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;

    template<typename Scalar>
    IndexType toCellIndex ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;
    template<typename Scalar>
    Eigen::VectorXi toCellIndices( const Eigen::Matrix<Scalar,3,Eigen::Dynamic> & pts, const int & level ) const;

    template<typename Scalar>
    bool isInBounds ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;
    template<typename Scalar>
    bool isInBoundsExceptCenter ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;


    // for Block Sparse

    template<typename Scalar>
    Eigen::Vector3i toBlockIndex ( const Eigen::Matrix<Scalar,3,1> & point, const int & level ) const;
    template<typename Scalar>
    void toBlockIndexAndIndexVector ( const Eigen::Matrix<Scalar,3,1> & point, const int & level, Eigen::Vector3i & block_idx_vec, Eigen::Vector3i & block_idx_inside_vec ) const;
    void fromCellIndexVectorToBlockIndices  ( const Eigen::Vector3i & cell_idx_vec, IndexType & block_idx, IndexType & block_idx_inside ) const;

    IndexType fromBlockIndexInsideVectorToBlockIdx ( const Eigen::Vector3i & block_idx_inside_vec ) const;
    IndexType fromBlockIndexVectorToBlockIdx ( const Eigen::Vector3i & block_idx_vec ) const;


    inline Eigen::Vector3i fromCellIndexVectorToBlockIdx( const Eigen::Vector3i & cell_idx_vec ) const
    {
        return (cell_idx_vec.cast<float>() / block_size).array().floor().cast<IndexType>();
    }

    inline Eigen::Vector3i fromCellIndexVectorToBlockInsideIdx( const Eigen::Vector3i & cell_idx_vec, const Eigen::Vector3i & block_idx_vec ) const
    {
        return cell_idx_vec - block_idx_vec * block_size;
    }

    static const MapParameters getMapParameters(const MapParameters & params );

    static const MapParameters getMapParameters(const float & size = MAP_SIDE, const IndexType & num_cells = NUM_CELLS, const LevelIndexType & levels = NUM_LEVELS, const bool & use_adaptive = true, const bool & omit_center = true );

private:
    MapParameters( const float & size, const int & num_cells, const LevelIndexType & levels, const Eigen::VectorXf & min_bounds, const Eigen::VectorXf & max_bounds, const bool & use_adaptive, const bool & omit_center );
};
