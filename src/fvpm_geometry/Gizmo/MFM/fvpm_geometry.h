#ifndef SWIFT_FVPM_GEOMETRY_GIZMO_MFM_H
#define SWIFT_FVPM_GEOMETRY_GIZMO_MFM_H

#include "hydro_part.h"
#include "inline.h"

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void fvpm_reset_centroids(
    size_t pind) {}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param pind Particle index.
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
 * @param pind Particle index (pindi).
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
 * @param pind Particle index (pindj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_right(
    size_t pind, const float *dx, const float w) {}

#endif /* SWIFT_FVPM_GEOMETRY_GIZMO_MFM_H */
