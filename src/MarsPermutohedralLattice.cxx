/**
BSD 3-Clause License

This file is part of the LiDAR MARS registration project.
https://github.com/AIS-Bonn/lidar_mars_registration
Copyright (c) 2021, Computer Science Institute VI, University of Bonn.

This file contains adapted code from:
https://graphics.stanford.edu/papers/permutohedral/
Copyright (c) 2010, Andrew Adams, Jongmin Baek and Abe Davis

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
#include "MarsPermutohedralLattice.h"
#include <absl/container/flat_hash_set.h>
#include <Tracy.hpp>

//#define OPT_THREAD_NUM 3

template <int P>
typename PermutohedralLattice<P>::MatHP PermutohedralLattice<P>::m_E = PermutohedralLattice<P>::getE();
template <int P>
typename PermutohedralLattice<P>::MatPH PermutohedralLattice<P>::m_Einv = PermutohedralLattice<P>::getE().transpose();

template <int P>
PermutohedralLattice<P>::PermutohedralLattice ( const MapParameters & params, const bool & scale_positions, const bool & use_barycentric  )
    : m_use_barycentric ( use_barycentric ), m_scale_positions ( scale_positions ), m_map_params ( params ), m_sigmas(params.m_num_levels, 0.), m_maps ( params.m_num_levels ), m_maps_surfel_info ( params.m_num_levels )
{
    //lattice properties
    m_scale_factor = compute_scale_factor();

    const Scalar spatial_dist_between_vertices = m_map_params.getCellSizeOnLevel(0); // X in meters, with 'nearest' -> each has influence/size of X meters
    const Scalar spatial_sigma = (spatial_dist_between_vertices/2.) / std::sqrt(3.);
    for ( int lvl = 0; lvl < m_map_params.m_num_levels; ++lvl )
    {
        const Scalar lvl_scale = std::pow<Scalar>(2.,Scalar(-lvl));
        set_sigmas(lvl_scale*spatial_sigma, lvl);
    }
}

template <int P>
typename PermutohedralLattice<P>::VecP PermutohedralLattice<P>::compute_scale_factor() const {
    VecP scaleFactor;
    constexpr float invStdDev = 1;
    // Compute parts of the rotation matrix E. (See pg.4-5 of paper.)
    for (int i = 0; i < pos_dim; ++i) {
        // the diagonal entries for normalization
        scaleFactor[i] = 1.0 / (sqrt((i + 1) * (i + 2))) * invStdDev;
    }
    return scaleFactor;
}

template <int P>
typename PermutohedralLattice<P>::MatHP PermutohedralLattice<P>::getE()
{
    MatHP E_left = MatHP::Zero(); //(m_pos_dim+1, m_pos_dim );
    MatP E_right = MatP::Zero(); //(m_pos_dim, m_pos_dim );
    //E left is has at the bottom a square matrix which has an upper triangular part of ones. Afterwards the whole E_left gets appended another row on top of all ones
    E_left.template bottomRows<pos_dim>().template triangularView<Eigen::Upper>().setOnes();
    //the diagonal of the bottom square is linearly incresing from [-1, -m_pos_dim]
    E_left.template bottomRows<pos_dim>().template diagonal().setLinSpaced(pos_dim,1,pos_dim);
    E_left.template bottomRows<pos_dim>().template diagonal()= -E_left.template bottomRows<pos_dim>().diagonal();
    //E_left has the first row all set to ones
    E_left.row(0).setOnes();
    //E right is just a diagonal matrix with entried in the diag set to 1/sqrt((d+1)(d+2)). Take into account that the d in the paper starts at 1 and we start at 0 so we add a +1 to diag_idx
    for(int diag_idx=0; diag_idx<pos_dim; ++diag_idx){
        E_right(diag_idx, diag_idx) =  1.0 / (sqrt((diag_idx + 1) * (diag_idx + 2))) ;
    }
    //rotate into H_d
    MatHP E = E_left*E_right;
    return E;
}

template <int P>
typename PermutohedralLattice<P>::MatPH PermutohedralLattice<P>::getEinv(const MatHP & E)
{
    return (E.transpose()*E).inverse() * E.transpose();
}

template <int P>
void PermutohedralLattice<P>::find_enclosing_simplex( const VecH & elevated, VecH & rem0, VecHi & rank ) const {
    // Find the closest 0-colored simplex through rounding
    // greedily search for the closest zero-colored lattice point
    int sum = 0;
    for (int i = 0; i <= pos_dim; ++i) {
        const float v = elevated[i] * (1.0 / hom_dim);
        const float up = std::ceil(v) * hom_dim;
        const float down = std::floor(v) * hom_dim;
        if (up - elevated[i] < elevated[i] - down) {
            rem0[i] = (short) up;
        } else {
            rem0[i] = (short) down;
        }
        sum += rem0[i];
    }
    sum /= hom_dim;

    // Find the simplex we are in and store it in rank (where rank describes what position coordinate i has in the sorted order of the features values)
    rank.setZero();
    for (int i = 0; i < pos_dim; ++i) {
        float di = elevated[i] - rem0[i];
        for (int j = i + 1; j <= pos_dim; ++j)
            if (di < elevated[j] - rem0[j])
                ++rank[i];
            else
                ++rank[j];
    }

    // If the point doesn't lie on the plane (sum != 0) bring it back
    for (int i = 0; i <= pos_dim; ++i) {
        rank[i] += sum;
        if (rank[i] < 0) {
            rank[i] += hom_dim;
            rem0[i] += hom_dim;
        } else if (rank[i] > pos_dim) {
            rank[i] -= hom_dim;
            rem0[i] -= hom_dim;
        }
    }
}

template <int P>
void PermutohedralLattice<P>::compute_barycentric_coordinates( const VecH & elevated, const VecH & rem0, const VecHi & rank, VecP2 & barycentric ) const {
    barycentric.setZero();
    // Compute the barycentric coordinates (p.10 in [Adams etal 2010])
    for (int i = 0; i <= pos_dim; ++i) {
        const float delta = (elevated[i] - rem0[i]) *  (1.0 / (hom_dim));
        barycentric[pos_dim - rank[i]] += delta;
        barycentric[hom_dim - rank[i]] -= delta;
    }
    // Wrap around
    barycentric[0] += 1.0 + barycentric[hom_dim];
}

template <int P>
int PermutohedralLattice<P>::update( const SurfelInfoVector & surfels )
{
    int num_points = 0;
    for ( auto & map : m_maps_surfel_info )
    {
        map.clear();
        map.reserve(surfels.size());
    }
    for ( const SurfelInfo & surfel : surfels )
    {
        if ( ! surfel.m_surfel ) continue;
        m_maps_surfel_info[surfel.m_level][surfel.m_index] = &surfel;
        num_points += surfel.m_surfel->getNumPoints();
    }
    return num_points;
}

template <int P>
void PermutohedralLattice<P>::removeCloudByScanId( const IndexType & scan_id, int & cellsWithRemovedPoints, int & pointsRemoved, int & cellsRemoved)
{
    cellsRemoved = 0;
    pointsRemoved = 0;
    cellsWithRemovedPoints = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)  reduction(+:cellsWithRemovedPoints) reduction(+:pointsRemoved) reduction(+:cellsRemoved)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        int cellsWithRemoved = 0;
        MapLevelType & map_level = m_maps[lvlIdx];
        std::vector<CellIndexType> cellIndicesToBeRemoved;
        for ( auto & cell_pair : map_level )
        {
            CellType & cell = cell_pair.second;
            const int prevPoints = cell.getNumPoints();
            const bool removedSome = cell.removeByScanId ( scan_id );
            cellsWithRemoved+=removedSome;
            if ( cell.isEmpty() )
                cellIndicesToBeRemoved.emplace_back(cell_pair.first);
            else
                cell.updateSurfel( true );
            const int afterPoints = cell.getNumPoints();
            pointsRemoved += prevPoints - afterPoints;
        }
        int num_erased = 0;
        for ( const CellIndexType & cellRemoveIdx : cellIndicesToBeRemoved )
            num_erased += map_level.erase( cellRemoveIdx );

        cellsWithRemovedPoints+=cellsWithRemoved;
        cellsRemoved += num_erased;
    }
}

template <int P>
void PermutohedralLattice<P>::getCellsAdaptive ( SurfelInfoVector & surfels, int & normal_level_cnter, int & coarser_level_cnter ) const
{
    normal_level_cnter = 0;
    coarser_level_cnter = 0;

    static Eigen::Matrix3Xi neighborAdd = Eigen::Matrix3Xi::Zero(3,8);
    static bool neighborAdded = false;
    if ( ! neighborAdded )
    {
        constexpr int neighbors = 2; // due to subdivision by 2, there are only 2 entries per axis
        int vi = 0;
        for ( IndexType iz = 0; iz < neighbors; ++iz )
        {
            for ( IndexType iy = 0; iy < neighbors; ++iy )
            {
                for ( IndexType ix = 0; ix < neighbors; ++ix )
                {
                    neighborAdd.col(vi) << ix,iy,iz;
                    ++vi;
                }
            }
        }
        neighborAdded = true;
    }

    const auto checkPlanarity = [](const CellType & cell ){
        return ( ( cell.m_surfel->eigen_values_(0) < MapParameters::first_ev_planarity_threshold )
              || ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold )
              || ( cell.m_surfel->eigen_values_(0) < MapParameters::plane_scale_factor * cell.m_surfel->eigen_values_(1)));};
    const auto checkIsValid = [] ( const CellType & cell ) { return cell.m_surfel->valid_ && cell.m_surfel->getNumPoints() >= 10.; };
    const auto checkNotInCenter = [] ( const CellType & cell ) { return ! cell.m_surfel_info.m_in_center; };

    absl::flat_hash_set<CellIndexType> prevCheckingCoarserFailed, prevUseCoarserLevel;
    for ( LevelIndexType lvlIdx = m_map_params.m_num_levels-1; lvlIdx >= 1; --lvlIdx ) // bottom up // top down
    {
        const auto & mapLevel = m_maps[lvlIdx];
        absl::flat_hash_set<CellIndexType> useCoarserLevel, checkingCoarserFailed;
        useCoarserLevel.reserve(mapLevel.size());
        checkingCoarserFailed.reserve(mapLevel.size());
        for ( const auto & pa : mapLevel )
        {
            const CellType & cell = pa.second;
            if ( ! cell.m_surfel ) LOG(FATAL) << "no surfel set."; // Should never happen
            // get coarser idx:
            const CellIndexType coarse_index ({{pa.first[0]>>1,pa.first[1]>>1,pa.first[2]>>1}});
            //     0     |     1     |     2 =>     a     |     a     |     x
            //  0  |  1  |  2  |  3  |  4  | =>  m  |  m  |  m  |  c  |  x  |  a
            // 0 1 | 2 3 | 4 5 | 6 7 | 8 9 | => c c | c c | c - | - - | c x | c c

            if ( prevCheckingCoarserFailed.count(pa.first) > 0 ){ // left case (c x -> x) on the right
                // this one was rejected on previous level
                checkingCoarserFailed.insert(coarse_index); // thus the one above also fails.
                continue;
            }

            //top-down: If this one is invalid, we do not need to check smaller levels since they will also be invalid!
            const bool isValid = checkIsValid ( cell );
            bool coarser_is_valid = true;
            bool checkCoarser = true;
            if ( isValid )
                checkCoarser = checkPlanarity( cell );

            if ( useCoarserLevel.count(coarse_index) > 0 )  // corresponding coarser was already added and should be used, skip this one
                continue;

            // finishes last case ( c c -> a ) on the right
            if ( checkingCoarserFailed.count(coarse_index) > 0 ) // corresponding coarser was already added and this one was checked!
            {
                checkCoarser = false;
            }


            if ( checkCoarser ) // should check upper one
            {
                // check neighbors:
                const VecI key3 (coarse_index[0]<<1,coarse_index[1]<<1,coarse_index[2]<<1);
                const VecHi key = (VecHi() << key3[0], key3[1], key3[2], -key3.sum()).finished();

                constexpr int max_num_neighbors = int(1<<pos_dim)+1;
                absl::InlinedVector<CellIndexType,max_num_neighbors> neighborKeys ( max_num_neighbors, CellIndexType({{0,0,0}}));
                neighborKeys[0] = {key[0],key[1],key[2]};

                VecHi n1_key, n2_key;
                // along each axis
                for (int remainder = 0, nIdx = 1; remainder <= pos_dim; ++remainder, nIdx+=2){
                    for (int k = 0; k < pos_dim; ++k){
                        n1_key[k] = key[k] + 1;
                        n2_key[k] = key[k] - 1;
                    }
                    // normally we can ignore n1_key[pos_dim], just in the last case not.
                    n1_key[remainder] = key[remainder] - pos_dim;
                    n2_key[remainder] = key[remainder] + pos_dim;

                    //const CellIndexType v1 ();
                    neighborKeys[nIdx] = {n1_key[0],n1_key[1],n1_key[2]};
                    neighborKeys[nIdx+1] = {n2_key[0],n2_key[1],n2_key[2]};
                }

                Eigen::Vector3f fine_normal = cell.m_surfel->normal_;
                int numChecked = 0;
                for ( int nIdx = 0; nIdx < max_num_neighbors; ++nIdx)
                {
                    const CellIndexType & nI = neighborKeys[nIdx];

                    if ( prevCheckingCoarserFailed.count(nI) > 0 ){ // inverse of left case ( a x -> x ) on the right
                        // this one previously failed!
                        checkCoarser = false;
                        break;
                    }

                    if ( nI != pa.first ) // only neighbors, not this one.
                    {
                        const auto it = mapLevel.find(nI);
                        if ( it == mapLevel.end() ) continue;
                        const auto & other_cell = it->second;
                        if ( ! other_cell.m_surfel ) continue;
                        if ( ! checkIsValid ( other_cell ) ) continue;
                        checkCoarser &= checkPlanarity ( other_cell );
                        fine_normal+= other_cell.m_surfel->normal_;
                        ++numChecked;
                    }
                    if ( ! checkCoarser ) break; // need to add all.
                }


                if ( checkCoarser )
                {
                    // check the actual coarser one
                    const auto it = m_maps[lvlIdx-1].find(coarse_index);
                    if ( it != m_maps[lvlIdx-1].end() )
                    {
                        const auto & coarser_cell = it->second;
                        coarser_is_valid = checkIsValid ( coarser_cell );
                        if ( coarser_cell.m_surfel && coarser_is_valid && checkPlanarity ( coarser_cell )
                             && ( numChecked == 0 || (numChecked > 0 && std::abs(fine_normal.normalized().dot(coarser_cell.m_surfel->normal_)) > 0.8) )
                             )
                        {
                            // use coarser one.
                            useCoarserLevel.insert(coarse_index);
                        }
                        else
                            checkCoarser = false; // should add the neighbors too
                    }
                }

                if ( ! checkCoarser ) // upper failed, add all neighbors
                {
                    // adding all available neighbors:
                    for ( int nIdx = 0; nIdx < max_num_neighbors; ++nIdx)
                    {
                        const CellIndexType & nI = neighborKeys[nIdx];
                        if ( prevCheckingCoarserFailed.count(nI) > 0 ) continue; // this one failed before, thus should not be added.
                        if ( nI != pa.first ) // only neighbors, not this one.
                        {
                            const auto it = mapLevel.find(nI);
                            if ( it == mapLevel.end() ) continue;
                            const auto & other_cell = it->second;
                            if ( ! other_cell.m_surfel ) continue;
                            if ( ! checkIsValid ( other_cell ) ) continue;
                            if ( other_cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                            surfels.emplace_back(other_cell.m_surfel, other_cell.m_surfel_info.m_center_s, lvlIdx, pa.first, other_cell.m_cov_per_scan.empty() ? -1 : other_cell.m_cov_per_scan.front().first, other_cell.m_surfel_info.m_in_center, other_cell.m_class);
                            if ( ! other_cell.m_surfel_info.m_in_center )
                                ++coarser_level_cnter;
                            else
                                ++normal_level_cnter;
                        }
                    }
                }
            }
            if ( ! checkCoarser && coarser_is_valid )
                checkingCoarserFailed.insert(coarse_index);

            if ( !checkCoarser && isValid ) // add current one
            {
                if ( cell.m_surfel->eigen_values_(1) < MapParameters::second_ev_degenerate_threshold ) continue;
                surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
                if ( ! cell.m_surfel_info.m_in_center )
                    ++coarser_level_cnter;
                else
                    ++normal_level_cnter;
            }
        }
        prevCheckingCoarserFailed = checkingCoarserFailed;
        prevUseCoarserLevel = useCoarserLevel;
    }

    int notInCenterL1 = 0;
    for( const auto & s : surfels )
        if ( s.m_level == 1 && !s.m_in_center )
            ++notInCenterL1;

    // add outer rim if available
    const auto & mapLevel = m_maps[0];
    for ( const auto & pa : mapLevel )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) LOG(FATAL) << "nullptr on surfel?"; // Should never happen

        //top-down: If this one is invalid, we do not need to check smaller levels since they will also be invalid!
        const bool isInvalid = ! checkIsValid ( cell );
        if ( isInvalid ) continue;

        if ( checkNotInCenter(cell) || prevUseCoarserLevel.count(pa.first) > 0 )
        {
            if ( ! cell.m_surfel_info.m_in_center )
                ++coarser_level_cnter;
            else
                ++normal_level_cnter;

            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, 0, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center,  cell.m_class);
            //continue; // top down, outside can be added directly
        }
    }

    std::vector<int> surfelsPerLevel(m_maps.size(),0);
    for ( const auto & s : surfels )
        ++surfelsPerLevel[s.m_level];
    LOG(INFO) << "Adaptive: per lvl: " << Eigen::Map<Eigen::VectorXi>(surfelsPerLevel.data(),surfelsPerLevel.size(), 1).transpose() << " prCoarser: " << prevUseCoarserLevel.size() << " notCL1: " << notInCenterL1;
}

template <int P>
void PermutohedralLattice<P>::getCellsOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

template <int P>
void PermutohedralLattice<P>::getCellsFusedOnLevel ( const IndexType & lvlIdx, SurfelInfoConstPtrVector & surfels, const bool & omit_center ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        surfels.emplace_back(cell.getSurfelInfo());
    }
}

template <int P>
void PermutohedralLattice<P>::getCellsScanWiseOnLevel ( const IndexType & lvlIdx, SurfelInfoVector & surfels, const bool & omit_center, const LevelIndexType & level_modifier  ) const
{
    for ( const auto & pa : m_maps[lvlIdx] )
    {
        const CellType & cell = pa.second;
        if ( ! cell.m_surfel ) continue;
        if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
        for ( const std::pair<IndexType,SurfelPtr> & surfelInfo : cell.m_cov_per_scan)
            surfels.emplace_back(surfelInfo.second, cell.m_surfel_info.m_center_s, lvlIdx+level_modifier, pa.first, surfelInfo.first, cell.m_surfel_info.m_in_center, cell.m_class);
    }
}

template <int P>
void PermutohedralLattice<P>::updateCells()
{
    ZoneScopedN("PermutohedralLattice::updateChangedCells");
    #pragma omp parallel for num_threads(OPT_THREAD_NUM)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        for ( auto & cell : map_level )
        {
            cell.second.updateSurfel(true);
        }
    }
}

template <int P>
void PermutohedralLattice<P>::compute_key( const VecH & rem0, const VecHi & rank, const int & closest_vertex_idx, VecI & key ) const
{
    for (int i = 0; i < pos_dim; ++i) {
        key[i] = static_cast<short>(rem0[i] + closest_vertex_idx);
        if (rank[i] > pos_dim - closest_vertex_idx)
            key[i] -= hom_dim;
    }
}

template <int P>
typename PermutohedralLattice<P>::VecI PermutohedralLattice<P>::compute_splat_key ( const Eigen::Vector3f & pt, const LevelIndexType & lvlIdx ) const
{
    VecI key;
    VecH rem0;
    VecHi rank;
    VecP2 barycentric;
    compute_splat_key(pt, lvlIdx, rem0, rank, barycentric);
    int closest_vertex_idx=-1;
    barycentric.maxCoeff(&closest_vertex_idx);
    compute_key ( rem0, rank, closest_vertex_idx, key );
    return key;
}

template <int P>
void PermutohedralLattice<P>::compute_splat_key ( const Eigen::Vector3f & pt, const LevelIndexType & lvlIdx, VecH & rem0, VecHi & rank, VecP2 & barycentric ) const
{
    Eigen::Vector3f pt_l = pt;
    if ( m_scale_positions )
        pt_l.array() /= m_sigmas[lvlIdx];
    const VecH elevated = m_E*pt_l.cast<Scalar>();
    find_enclosing_simplex( elevated, rem0, rank );
    compute_barycentric_coordinates( elevated, rem0, rank, barycentric );
}

template <int P>
typename PermutohedralLattice<P>::VecP PermutohedralLattice<P>::compute_vertex_position ( const VecI & key, const LevelIndexType & lvlIdx ) const
{
    VecHi vKey = VecHi::Zero();
    vKey.template head<3>() = key;
    vKey[pos_dim] = -key.template sum();

    VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
    if ( m_scale_positions )
        new_cell_center *= m_sigmas[lvlIdx];
    return new_cell_center;
}

template <int P>
void PermutohedralLattice<P>::getOriginCellShift( const Eigen::Vector3f & rel_translation, const LevelIndexType & lvl, Eigen::Vector3f & map_translation, Eigen::Vector3i & num_cells ) const
{
    const Eigen::Vector3f origin = Eigen::Vector3f::Zero();
    const VecI origin_key = compute_splat_key ( origin, lvl );
    const VecP origin_vertex_pos = compute_vertex_position ( origin_key, lvl );

    const VecI rel_key = compute_splat_key ( rel_translation, lvl );
    const VecP rel_vertex_pos = compute_vertex_position ( rel_key, lvl );

    num_cells = (rel_key-origin_key).template head<3>();
    map_translation = (rel_vertex_pos-origin_vertex_pos).template head<3>().template cast<float>();
}


template <int P> template<typename SurfelInfoT>
int PermutohedralLattice<P>::addCells ( const std::vector<SurfelInfoT> & surfels, const Sophus::SE3d & orig_pose_s1s2, const LevelIndexType & level_modifier )
{
    ZoneScopedN("PermutohedralLattice::addCells");
    Sophus::SE3d pose_s1s2 = orig_pose_s1s2;
    if ( pose_s1s2.unit_quaternion().isApprox(Eigen::Quaterniond::Identity()) )
        pose_s1s2.so3().setQuaternion(Eigen::Quaterniond::Identity());
    if ( pose_s1s2.translation().isApprox(Eigen::Vector3d::Zero()))
        pose_s1s2.translation().setZero();

    int num_points = 0;
    Sophus::SE3d n_pose_s1s2 = pose_s1s2;
    const size_t num_surfels = surfels.size();
    for ( size_t idx = 0; idx < num_surfels; ++idx)
    {
        if constexpr ( std::is_pointer<SurfelInfoT>::value )
            num_points += addCellSurfel ( *(surfels[idx]), pose_s1s2, level_modifier, n_pose_s1s2);
        else
            num_points += addCellSurfel ( surfels[idx], pose_s1s2, level_modifier, n_pose_s1s2);
    }
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(INFO) << "num_points: " << num_points << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose(); // << " cp: " << cp;
    return num_points;
}

template <int P>
inline int PermutohedralLattice<P>::addCellSurfel ( const SurfelInfo & surfel_info, const Sophus::SE3d & pose_s1s2, const LevelIndexType & level_modifier, Sophus::SE3d & n_pose_s1s2 )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;

    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & center_s = surfel_info.m_center_s;
    const Eigen::Vector3f meanPt = (pose_s1s2 * (surfel->mean_ + center_s).template cast<double>()).template cast<float>();

    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    Eigen::Vector3f pt_l = meanPt;
    if ( m_scale_positions )
        pt_l.array() /= m_sigmas[lvlIdx];

    VecI key;
    VecHi vKey = VecHi::Zero();
    VecH rem0;
    VecHi rank;
    VecP2 barycentric;

    const VecH elevated = m_E*pt_l.cast<Scalar>();
    find_enclosing_simplex( elevated, rem0, rank );
    compute_barycentric_coordinates( elevated, rem0, rank, barycentric );

    int num_points = 0;
    if ( m_use_barycentric )
    {
        for (int remainder = 0; remainder <= pos_dim; ++remainder) {
            // Compute the location of the lattice point explicitly (all but
            // the last coordinate - it's redundant because they sum to zero)
            vKey[pos_dim] = 0;
            for (int i = 0; i < pos_dim; ++i) {
                key[i] = static_cast<short>(rem0[i] + remainder);
                if (rank[i] > pos_dim - remainder)
                    key[i] -= hom_dim;
                vKey[i] = key[i];
                vKey[pos_dim] -= key[i];
            }

            // Retrieve pointer to the value at this vertex.
            const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
            CellType & cell = map_level[vertex_idx];
            if ( cell.isEmpty() )
            {
                VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
                if ( m_scale_positions )
                    new_cell_center *= m_sigmas[lvlIdx];
                cell.create( new_cell_center, (pose_s1s2.translation().template cast<Scalar>()-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1 ), lvlIdx, vertex_idx );
            }
            n_pose_s1s2.translation() =  (meanPt - surfel->mean_ - cell.m_surfel_info.m_center_s).template cast<double>();
            // TODO: weight with barycentric coord.
            cell.addSurfel(surfel, n_pose_s1s2, surfel_info.m_scan_id, surfel_info.m_class);
            num_points += surfel->getNumPoints();
        }
    }
    else
    {
        int closest_vertex_idx=-1;
        barycentric.maxCoeff(&closest_vertex_idx);
        //splat
        vKey[pos_dim] = 0;
        for (int i = 0; i < pos_dim; ++i) {
            key[i] = static_cast<short>(rem0[i] + closest_vertex_idx);
            if (rank[i] > pos_dim - closest_vertex_idx)
                key[i] -= hom_dim;
            vKey[i] = key[i];
            vKey[pos_dim] -= key[i];
        }
        const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
        CellType & cell = map_level[vertex_idx];
        if ( cell.isEmpty() )
        {
            VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
            if ( m_scale_positions )
                new_cell_center.array() *= m_sigmas[lvlIdx];
            // add this one.
            cell.create( new_cell_center, (pose_s1s2.translation().template cast<Scalar>()-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1 ), lvlIdx, vertex_idx );
        }

        // meanPt = pose_s1s2 * (mean_old + center_old) = mean_new + center_new
        // mean_new = meanPt - center_new - surfel mean
        n_pose_s1s2.translation() =  (meanPt - surfel->mean_ - cell.m_surfel_info.m_center_s).template cast<double>();
        cell.addSurfel(surfel, n_pose_s1s2, surfel_info.m_scan_id, surfel_info.m_class);
        num_points = surfel->getNumPoints();
    }
    return num_points;
}

template<int P>
inline int PermutohedralLattice<P>::addCellSurfelOnGrid( const SurfelInfo & surfel_info, const Eigen::Vector3f & t_s1s2, const LevelIndexType & level_modifier, const Sophus::SE3d & id )
{
    const SurfelPtr & surfel = surfel_info.m_surfel;
    if ( ! surfel ) return 0;
    const LevelIndexType lvlIdx = surfel_info.m_level + level_modifier;
    MapLevelType & map_level = m_maps[lvlIdx];
    const Eigen::Vector3f & old_center = surfel_info.m_center_s;
    const Eigen::Vector3f surfel_mean = surfel->mean_ + old_center;
    const Eigen::Vector3f meanPt = surfel_mean + t_s1s2;
    if ( ! m_map_params.isInBounds( meanPt, lvlIdx ) ) return 0;

    Eigen::Vector3f pt_l = meanPt;
    if ( m_scale_positions )
        pt_l.array() /= m_sigmas[lvlIdx];

    VecI key;
    VecHi vKey = VecHi::Zero();
    VecH rem0;
    VecHi rank;
    VecP2 barycentric;

    const VecH elevated = m_E*pt_l.cast<Scalar>();
    find_enclosing_simplex( elevated, rem0, rank );
    compute_barycentric_coordinates( elevated, rem0, rank, barycentric );
    int num_points = 0;
    if ( m_use_barycentric )
    {
        for (int remainder = 0; remainder <= pos_dim; ++remainder) {
            // Compute the location of the lattice point explicitly (all but
            // the last coordinate - it's redundant because they sum to zero)

            vKey[pos_dim] = 0;
            for (int i = 0; i < pos_dim; ++i) {
                key[i] = static_cast<short>(rem0[i] + remainder);
                if (rank[i] > pos_dim - remainder)
                    key[i] -= hom_dim;
                vKey[i] = key[i];
                vKey[pos_dim] -= key[i];
            }

            // Retrieve pointer to the value at this vertex.
            const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
            CellType & cell = map_level[vertex_idx];
            if ( cell.isEmpty() )
            {
                VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
                if ( m_scale_positions )
                    new_cell_center *= m_sigmas[lvlIdx];
                cell.create( new_cell_center, (t_s1s2-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, vertex_idx );
            }
            // TODO: weight with barycentric coord.
            cell.addSurfel(surfel, id, surfel_info.m_scan_id, surfel_info.m_class);
            num_points += surfel->getNumPoints();
        }
    }
    else
    {
        int closest_vertex_idx=-1;
        barycentric.maxCoeff(&closest_vertex_idx);
        //splat
        vKey[pos_dim] = 0;
        for (int i = 0; i < pos_dim; ++i) {
            key[i] = static_cast<short>(rem0[i] + closest_vertex_idx);
            if (rank[i] > pos_dim - closest_vertex_idx)
                key[i] -= hom_dim;
            vKey[i] = key[i];
            vKey[pos_dim] -= key[i];
        }
        const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
        CellType & cell = map_level[vertex_idx];
        if ( cell.isEmpty() )
        {
            VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
            if ( m_scale_positions )
                new_cell_center.array() *= m_sigmas[lvlIdx];
            // add this one.
            cell.create( new_cell_center, (t_s1s2-meanPt).normalized(), m_map_params.isInBounds( new_cell_center, lvlIdx+1), lvlIdx, vertex_idx );
        }
        cell.addSurfel(surfel, id, surfel_info.m_scan_id, surfel_info.m_class);
        num_points = surfel->getNumPoints();
    }
    return num_points;
}

template <int P> template<typename SurfelInfoT>
int PermutohedralLattice<P>::addCellsOnGrid ( const std::vector<SurfelInfoT> & surfels, const Eigen::Vector3d & orig_t_s1s2, const LevelIndexType & level_modifier )
{
    // Assumes surfels are on the same grid with a shift of n_i per dimension. -> no rotation or fractions within translation
    ZoneScopedN("PermutohedralLattice::addCells");

    Eigen::Vector3f t_s1s2 = orig_t_s1s2.cast<float>();//Eigen::round(orig_t_s1s2.template cast<float>().array());
    const Sophus::SE3d id;
    int num_points = 0;
    const size_t num_surfels = surfels.size();
    for ( size_t idx = 0; idx < num_surfels; ++idx)
    {
        if constexpr ( std::is_pointer<SurfelInfoT>::value )
            num_points += addCellSurfelOnGrid ( *(surfels[idx]), t_s1s2, level_modifier, id);
        else
            num_points += addCellSurfelOnGrid ( surfels[idx], t_s1s2, level_modifier, id);
    }
    Eigen::VectorXi cellsFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    Eigen::VectorXi cellsMaxFilled = Eigen::VectorXi::Zero(m_maps.size(),1);
    for ( size_t lvlIdx = 0; lvlIdx < m_maps.size(); ++lvlIdx )
    {
        cellsFilled[lvlIdx] = m_maps[lvlIdx].size();
        cellsMaxFilled[lvlIdx] = m_maps[lvlIdx].capacity();
    }
    LOG(INFO) << "num_points: " << num_points << " cellsPerLevel: " << cellsFilled.transpose() << " max: " << cellsMaxFilled.transpose();
    return num_points;
}

template <int P>
int PermutohedralLattice<P>::addCloud(MarsMapPointCloud::Ptr cloud, const Sophus::SE3d & pose_s1s2, const Eigen::Vector3d & origin )
{
    ZoneScopedN("PermutohedralLattice::addCloud");
    absl::InlinedVector<absl::flat_hash_set<CellIndexType>,MapParameters::maxAllowedMapLevels> changedCells ( m_map_params.m_num_levels );
    const Eigen::Vector3f originf = origin.cast<float>();
    const Eigen::Matrix3Xf pts_s = (pose_s1s2.so3().cast<float>().matrix() * cloud->m_points).colwise() + pose_s1s2.translation().cast<float>();
    const int scan_id = cloud->single_scan() ? cloud->m_scan_id(0) : -1;
    const int num_pts = cloud->size();
    int num_points = 0;
    #pragma omp parallel for num_threads(OPT_THREAD_NUM) reduction(+:num_points)
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        absl::flat_hash_set<CellIndexType> & changed_level_cells = changedCells[lvlIdx];
        changed_level_cells.reserve(0.1*num_pts);

        VecI key;
        VecHi vKey = VecHi::Zero();
        VecH rem0;
        VecHi rank;
        VecP2 barycentric;
        Eigen::Vector3f pt_l;
        for ( int pointIdx = 0; pointIdx < num_pts; ++pointIdx )
        {
            const Eigen::Vector3f pt_s = pts_s.col(pointIdx);
            const bool isInside = m_map_params.isInBounds( pt_s, lvlIdx );
            if ( ! isInside ) continue;

            pt_l = pt_s;
            if ( m_scale_positions )
                pt_l.array() /= m_sigmas[lvlIdx];
            vKey.setZero();

            const VecH elevated = m_E*pt_l.cast<Scalar>();
            find_enclosing_simplex( elevated, rem0, rank );
            compute_barycentric_coordinates( elevated, rem0, rank, barycentric );
            if ( m_use_barycentric)
            {
                // splat to all the vertices in the simplex

                for (int remainder = 0; remainder <= pos_dim; ++remainder) {
                    // Compute the location of the lattice point explicitly (all but
                    // the last coordinate - it's redundant because they sum to zero)

                    vKey[pos_dim] = 0;
                    for (int i = 0; i < pos_dim; ++i) {
                        key[i] = static_cast<short>(rem0[i] + remainder);
                        if (rank[i] > pos_dim - remainder)
                            key[i] -= hom_dim;
                        vKey[i] = key[i];
                        vKey[pos_dim] -= key[i];
                    }

                    // Retrieve pointer to the value at this vertex.
                    const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
                    CellType & cell = map_level[vertex_idx];
                    if ( cell.isEmpty() )
                    {
                        VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
                        if ( m_scale_positions )
                            new_cell_center.array() *= m_sigmas[lvlIdx];
                        cell.create( new_cell_center, (pt_s-originf).normalized(), m_map_params.isInBounds<float>( new_cell_center, lvlIdx+1 ), lvlIdx, vertex_idx );
                    }
                    if ( !changed_level_cells.contains(vertex_idx) )
                    {
                        changed_level_cells.insert(vertex_idx);
                        cell.addNewScan( scan_id, (pt_s-originf).normalized() );
                    }
                    cell.addPoint( pt_s );
                    ++num_points;
                }
            } else {
                //splat only on the closest vertex
                int closest_vertex_idx=-1;
                barycentric.maxCoeff(&closest_vertex_idx);

                //splat
                vKey[pos_dim] = 0;
                for (int i = 0; i < pos_dim; ++i) {
                    key[i] = static_cast<short>(rem0[i] + closest_vertex_idx);
                    if (rank[i] > pos_dim - closest_vertex_idx)
                        key[i] -= hom_dim;
                    vKey[i] = key[i];
                    vKey[pos_dim] -= key[i];
                }

                // Retrieve pointer to the value at this vertex.
                const CellIndexType vertex_idx ({{key[0],key[1],key[2]}});
                CellType & cell = map_level[vertex_idx];
                if ( cell.isEmpty() )
                {
                    VecP new_cell_center = (m_Einv * vKey.template cast<Scalar>());
                    if ( m_scale_positions )
                        new_cell_center.array() *= m_sigmas[lvlIdx];
                    cell.create( new_cell_center, (pt_s-originf).normalized(), m_map_params.isInBounds<float>( new_cell_center, lvlIdx+1 ), lvlIdx, vertex_idx );
                }

                if ( !changed_level_cells.contains(vertex_idx) )
                {
                    changed_level_cells.insert(vertex_idx);
                    cell.addNewScan( scan_id, (pt_s-originf).normalized() );
                }
                cell.addPoint( pt_s );
                ++num_points;
            }
        }
    }
    int updated = 0;
    {
        ZoneScopedN("PermutohedralLattice::addCloud::UpdateSurfel");
        for ( size_t lvlIdx = 0; lvlIdx < changedCells.size(); ++lvlIdx )
        {
            const absl::flat_hash_set<CellIndexType> & changed_level_cells = changedCells[lvlIdx];
            auto & map_level = m_maps[lvlIdx];
            for ( const CellIndexType & cellIdx : changed_level_cells )
            {
                const auto & cellIt = map_level.find(cellIdx);
                if ( cellIt != map_level.end() )
                {
                    CellType & cell = cellIt->second;
                    cell.updateSurfel( );
                    ++updated;
                }
            }
        }
    }
    LOG(INFO) << "num_points: " << num_points << " ( " << num_pts << " ) updatedCells: " << updated;
    return num_points;
}

template <int P>
bool PermutohedralLattice<P>::getSensorCell ( const Eigen::Vector3f & pt_s, const LevelIndexType & search_lvl, SurfelInfoConstPtrVector & cellPtrs, const IndexType & neighbors ) const
{
    ZoneScopedN("PermutohedralLattice::getSensorCell");
    cellPtrs.clear();
    if ( neighbors > 1 ) LOG(FATAL) << "Currently not supported for higher hop neigbors.";

    if ( search_lvl < 0 || search_lvl >= LevelIndexType(m_maps.size()) ) { LOG(INFO) << "ooL? " << pt_s.transpose() << " lvl: " << search_lvl; return false; } // check lvl bounds
    if ( ! m_map_params.isInBounds ( pt_s , search_lvl ) ) {  LOG(INFO) << "oob: " << pt_s.transpose() << " lvl: " << search_lvl; return false; }
    cellPtrs.reserve(9);

    // find the corresponding simplex:
    Eigen::Vector3f pt_l = pt_s;
    if ( m_scale_positions )
        pt_l.array() /= m_sigmas[search_lvl];

    VecHi key;
    VecH rem0;
    VecHi rank;
    VecP2 barycentric;
    const VecH elevated = m_E*pt_l.cast<Scalar>();
    find_enclosing_simplex( elevated, rem0, rank );
    compute_barycentric_coordinates( elevated, rem0, rank, barycentric );

    //splat only on the closest vertex
    int closest_vertex_idx=-1;
    barycentric.maxCoeff(&closest_vertex_idx);
    //splat
    for (int i = 0; i < pos_dim; ++i) {
        key[i] = static_cast<short>(rem0[i] + closest_vertex_idx);
        if (rank[i] > pos_dim - closest_vertex_idx)
            key[i] -= hom_dim;
    }
    key[pos_dim] = - key.template head<3>().sum();
    const auto & map_level = m_maps_surfel_info[search_lvl];
    VecHi n1_key;
    VecHi n2_key;
    // along each axis
    for (int remainder = 0; remainder <= pos_dim; ++remainder){
        for (int k = 0; k < pos_dim; ++k) {
            n1_key[k] = key[k] + 1;
            n2_key[k] = key[k] - 1;
        }
        n1_key[remainder] = key[remainder] - pos_dim;
        n2_key[remainder] = key[remainder] + pos_dim;

        const CellIndexType v1 ({{n1_key[0],n1_key[1],n1_key[2]}});
        const auto it1 = map_level.find( v1 );
        if ( it1 != map_level.end() )
        {
            cellPtrs.emplace_back(it1->second);
        }
        const CellIndexType v2 ({{n2_key[0],n2_key[1],n2_key[2]}});
        const auto it2 = map_level.find( v2 );
        if ( it2 != map_level.end() )
        {
            cellPtrs.emplace_back(it2->second);
        }
    }
    return !cellPtrs.empty();
}


template <int P>
int PermutohedralLattice<P>::getSensorCells ( const SurfelInfoVector & data, const IndexType & neighbors ) const
{
    ZoneScopedN("PermutohedralLattice::getSensorCells");
    int numNeighbors = 0;
    std::vector<SurfelInfoConstPtrVector> cnPtrs (data.size());
    for ( size_t data_idx = 0; data_idx < data.size(); ++data_idx )
    {
        const LevelIndexType & search_lvl = data[data_idx].m_level;
        const Eigen::Vector3f pt_s = data[data_idx].m_surfel->mean_ + data[data_idx].m_center_s;
        SurfelInfoConstPtrVector & cellPtrs = cnPtrs[data_idx];

        getSensorCell(pt_s, search_lvl, cellPtrs, neighbors );
        numNeighbors+= cellPtrs.size();
    }
    return numNeighbors;
}

template <int P>
void PermutohedralLattice<P>::allocate()
{
    ZoneScopedN("PermutohedralLattice::allocate");
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        MapLevelType & map_level = m_maps[lvlIdx];
        if ( map_level.empty() )
            map_level.reserve(std::min(1e5,std::pow(m_map_params.getNumCells(lvlIdx),3)));
        if ( m_maps_surfel_info[lvlIdx].empty() )
            m_maps_surfel_info[lvlIdx].reserve(std::min(1e5,std::pow(m_map_params.getNumCells(lvlIdx),3)));
    }
}

template <int P>
void PermutohedralLattice<P>::clear()
{
    for ( LevelIndexType lvl = 0; lvl < LevelIndexType(m_maps.size()); ++lvl )
    {
        m_maps[lvl].clear();
        m_maps_surfel_info[lvl].clear();
    }
}

template <int P>
SurfelInfoVector PermutohedralLattice<P>::iterate() const
{
    const bool should_omit_center = true;
    SurfelInfoVector surfels;
    surfels.clear();
    surfels.reserve(1e5);
    for ( LevelIndexType lvlIdx = 0; lvlIdx < m_map_params.m_num_levels; ++lvlIdx )
    {
        const bool omit_center = should_omit_center && lvlIdx < m_map_params.m_num_levels-1;
        for ( const auto & pa : m_maps[lvlIdx] )
        {
            const CellType & cell = pa.second;
            if ( ! cell.m_surfel ) continue;
            if ( omit_center && cell.m_surfel_info.m_in_center ) continue;
            surfels.emplace_back(cell.m_surfel, cell.m_surfel_info.m_center_s, lvlIdx, pa.first, cell.m_cov_per_scan.empty() ? -1 : cell.m_cov_per_scan.front().first, cell.m_surfel_info.m_in_center, cell.m_class);
        }
    }
    return surfels;
}

template <int P>
void PermutohedralLattice<P>::set_sigmas( const double & sigma, const int & level ) //std::initializer_list<std::pair<double,int>> & sigmas_list, const int & level)
{
    if ( level < 0 || level >= int(m_maps.size()) ) return;
    m_sigmas[level] = sigma;
}

template <int P>
typename PermutohedralLattice<P>::MatX PermutohedralLattice<P>::get_key_vertices( const int & level )
{
    if ( level < 0 || level >= int(m_maps.size()) ) return MatX::Zero(pos_dim,0);
    int idx = 0;
    MatX key_vertices = MatX::Zero(pos_dim, m_maps[level].size());
    for ( const auto & p : m_maps[level] )
    {
        key_vertices.col(idx) = p.second.m_surfel_info.m_center_s.cast<float>();
        ++idx;
    }
    return key_vertices;
}

template <int P>
std::vector<Surfel::Ptr> PermutohedralLattice<P>::get_key_surfels( const int & level )
{
    if ( level < 0 || level >= int(m_maps.size()) ) return std::vector<Surfel::Ptr>();
    std::vector<Surfel::Ptr> surfels; surfels.reserve(m_maps[level].size());
    for ( const auto & p : m_maps[level] )
    {
        surfels.emplace_back( p.second.m_surfel );
    }
    return surfels;
}

template class PermutohedralLattice<3>;
template int Lattice::addCellsOnGrid<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Eigen::Vector3d & , const LevelIndexType &);
template int Lattice::addCellsOnGrid<SurfelInfo>( const SurfelInfoVector & , const Eigen::Vector3d & , const LevelIndexType & );
template int Lattice::addCells<SurfelInfoConstPtr>( const SurfelInfoConstPtrVector & , const Sophus::SE3d & , const LevelIndexType & );
template int Lattice::addCells<SurfelInfo>( const SurfelInfoVector & , const Sophus::SE3d & , const LevelIndexType & );
