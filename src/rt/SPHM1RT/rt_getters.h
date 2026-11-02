/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Tsang Keung Chan (chantsangkeung@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_SPHM1RT_RT_GETTERS_H
#define SWIFT_SPHM1RT_RT_GETTERS_H

/**
 * @file src/rt/SPHM1RT/rt_getters.h
 * @brief Independent getter functions for the SPHM1RT scheme
 */

#include "inline.h"

#include <stddef.h>

/**
 * @brief Returns the comoving radiation speed of a particle
 * @param pind Index of the particle.
 * @param a Current scale factor.
 * @return comoving reduced speed of light
 */
__attribute__((always_inline)) INLINE static float rt_get_comoving_cred(
    const size_t pind, float a) {
  return part_get_rt_data_p(pind)->params.cred_phys / a;
}

/**
 * @brief Returns the physical radiation speed of a particle
 *
 * @param pind Index of the particle.
 * @param a Current scale factor.
 * @return physical reduced speed of light
 *
 */
__attribute__((always_inline)) INLINE static float rt_get_physical_cred(
    const size_t pind, float a) {
  return part_get_rt_data_p(pind)->params.cred_phys;
}

/**
 * @brief Returns the comoving radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param pind Index of the particle.
 * @param urad The comoving radiation energy per mass.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_comoving_urad_multifrequency(const size_t pind,
                                    float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = part_get_const_rt_data_p(pind)->conserved[g].urad;
  }
}

/**
 * @brief Returns the physical radiation energy per mass of a particle
 * (note that the comoving and physical energy per mass are the same in our
 * convention)
 *
 * @param pind Index of the particle.
 * @param cosmo Cosmology data structure.
 * @param urad The physical radiation energy.
 *
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_urad_multifrequency(const size_t pind,
                                    const struct cosmology *cosmo,
                                    float urad[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    urad[g] = part_get_const_rt_data_p(pind)->conserved[g].urad;
  }
}

/**
 * @brief Returns the comoving radiation flux per gas density of a particle
 * (note that the comoving and physical flux per gas density are the same in our
 * convention)
 *
 * @param pind Index of the particle.
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_comoving_frad_multifrequency(const size_t pind,
                                    float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = part_get_const_rt_data_p(pind)->conserved[g].frad[0];
    fradtemp[g][1] = part_get_const_rt_data_p(pind)->conserved[g].frad[1];
    fradtemp[g][2] = part_get_const_rt_data_p(pind)->conserved[g].frad[2];
  }
}

/**
 * @brief Returns the physical radiation flux per gas density of a particle
 * (note that the comoving and physical flux per gas density are the same in our
 * convention)
 *
 * @param pind Index of the particle.
 * @param cosmo Cosmology data structure
 * @param fradtemp The comoving radiation flux per gas density
 */
__attribute__((always_inline)) INLINE static void
rt_get_physical_frad_multifrequency(const size_t pind,
                                    const struct cosmology *cosmo,
                                    float fradtemp[RT_NGROUPS][3]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[g][0] = part_get_const_rt_data_p(pind)->conserved[g].frad[0] * cosmo->a;
    fradtemp[g][1] = part_get_const_rt_data_p(pind)->conserved[g].frad[1] * cosmo->a;
    fradtemp[g][2] = part_get_const_rt_data_p(pind)->conserved[g].frad[2] * cosmo->a;
  }
}

#endif /* SWIFT_SPHM1RT_RT_GETTERS_H */
