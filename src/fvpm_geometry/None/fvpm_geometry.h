#ifndef SWIFT_FVPM_GEOMETRY_NONE_H
#define SWIFT_FVPM_GEOMETRY_NONE_H

/**
 * @file None/fvpm_geometry.h
 * @brief Functions related to the Gizmo FVPM geometry struct collection,
 * in particular the collection of the data required for the matrix needed
 * for gradients. Empty definitions for when FVPM are unused.
 */

#include "part.h"
#include "swift_intrinsics.h"

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void fvpm_reset_centroids(
    size_t pind) {}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param p Particle.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void fvpm_normalise_centroid(
    size_t pind, const float wcount) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param p Particle (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_left(
    size_t pind, const float *dx, const float w) {}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param p Particle (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_right(
    size_t pind, const float *dx, const float w) {}

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
fvpm_part_geometry_well_behaved(size_t pind) {
  return 0;
}

/**
 * @brief Collect the data needed for the matrix construction.
 */
__attribute__((always_inline)) INLINE static void
fvpm_accumulate_geometry_and_matrix(size_t pind, const float wi, const float dx[3]) {}

__attribute__((always_inline)) INLINE static void fvpm_geometry_init(size_t pind) {}

/**
 * @brief Finish the computation of the matrix.
 *
 * @param p the particle to work on
 * @param ihdim 1/h^{dim}
 */
__attribute__((always_inline)) INLINE static void
fvpm_compute_volume_and_matrix(size_t pind, const float ihdim) {}

#endif /* SWIFT_FVPM_GEOMETRY_NONE_H */
