/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_CHEMISTRY_GEAR_DIFFUSION_DEBUG_H
#define SWIFT_CHEMISTRY_GEAR_DIFFUSION_DEBUG_H

#include "error.h"
#include "inline.h"

__attribute__((always_inline)) INLINE static void chemistry_debug_particle(const size_t pind) {

  const struct chemistry_part_data *restrict *cpd = part_get_const_chemistry_data_p(pind);

  warning("[PID%lld] chemistry_part_data:", part_get_id(pind));
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] metal_mass[%i]: %.3e", part_get_id(pind), i,
           cpd->metal_mass[i]);
  }
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] smoothed_metal_mass_fraction[%i]: %.3e", part_get_id(pind), i,
           cpd->smoothed_metal_mass_fraction[i]);
  }
  warning("[PID%lld] diff_coef: %.3e",part_get_id(pind), cpd->diff_coef);
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] metal_mass_dt[%i]: %.3e", part_get_id(pind), i,
           cpd->metal_mass_dt[i]);
  }
  warning("[PID%lld S: [[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e]]",
         part_get_id(pind),
         cpd->S[0][0], cpd->S[0][1], cpd->S[0][2],
         cpd->S[1][0], cpd->S[1][1], cpd->S[1][2],
         cpd->S[2][0], cpd->S[2][1], cpd->S[2][2]);
}

#endif /* SWIFT_CHEMISTRY_GEAR_DIFFUSION_DEBUG_H */
