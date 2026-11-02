#ifndef HYDRO_DATA_H
#define HYDRO_DATA_H

#include "inline.h"
#include "hydro_part.h"

/**
 * @brief swap the data of two particles.
 *
 * @param pindi index of first particle
 * @param pindj index of second particle
 */
__attribute__((__always_inline__)) INLINE void hydro_data_swap_particles(const size_t pindi, const size_t pindj) {

  struct part tmp = global_part_arrays._part[pindi];
  global_part_arrays._part[pindi] = global_part_arrays._part[pindj];
  global_part_arrays._part[pindj] = tmp;

}

#endif

