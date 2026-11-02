#ifndef SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H
#define SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H

#include "const.h"
#include "hydro_part.h"
#include "inline.h"
#include "kernel_hydro.h"

/**
 * @brief Reset the variables used to store the centroid; used for the velocity
 * correction.
 */
__attribute__((always_inline)) INLINE static void fvpm_reset_centroids(size_t pind) {

  struct fvpm_geometry_struct *restrict geometry = part_get_fvpm_geometry_p(pind);
  geometry->centroid[0] = 0.0f;
  geometry->centroid[1] = 0.0f;
  geometry->centroid[2] = 0.0f;
}

/**
 * @brief Normalise the centroids after the density loop.
 *
 * @param pind Particle index.
 * @param wcount Wcount for the particle. This is an explicit argument, so that
 * it is clear from the code that wcount needs to be normalised by the time it
 * is used here.
 */
__attribute__((always_inline)) INLINE static void fvpm_normalise_centroid(size_t pind, const float wcount) {

  struct fvpm_geometry_struct *restrict geometry = part_get_fvpm_geometry_p(pind);
  const float norm = kernel_norm / wcount;
  geometry->centroid[0] *= norm;
  geometry->centroid[1] *= norm;
  geometry->centroid[2] *= norm;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the left particle in the neighbour interaction.
 *
 * @param pind Particle index. (pi).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pj->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_left(
    size_t pind, const float *dx, const float w) {

  struct fvpm_geometry_struct *restrict geometry = part_get_fvpm_geometry_p(pind);
  geometry->centroid[0] -= dx[0] * w;
  geometry->centroid[1] -= dx[1] * w;
  geometry->centroid[2] -= dx[2] * w;
}

/**
 * @brief Update the centroid with the given contribution, assuming the particle
 * acts as the right particle in the neighbour interaction.
 *
 * @param pind Particle index. (pj).
 * @param dx Distance vector between the particle and its neighbour (dx = pi->x
 * - pj->x).
 * @param w Kernel value at position pi->x.
 */
__attribute__((always_inline)) INLINE static void fvpm_update_centroid_right(
    size_t pind, const float *dx, const float w) {

  struct fvpm_geometry_struct *restrict geometry = part_get_fvpm_geometry_p(pind);
  geometry->centroid[0] += dx[0] * w;
  geometry->centroid[1] += dx[1] * w;
  geometry->centroid[2] += dx[2] * w;
}

#endif /* SWIFT_FVPM_GEOMETRY_GIZMO_MFV_H */
