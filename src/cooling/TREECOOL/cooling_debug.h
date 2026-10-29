/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (mschaller@lorentz.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_TREECOOL_DEBUG_H
#define SWIFT_COOLING_TREECOOL_DEBUG_H

/**
 * @file src/cooling/TREECOOL/cooling_debug.h
 * @brief Debugging routines related to the TREECOOL cooling function (Katz,
 * Weinberg & Hernquist 1996).
 */

#include "error.h"
#include "inline.h"

#include <stddef.h>

/**
 * @brief Print the cooling data of a particle.
 *
 * @param pind Index of the hydro particle.
 */
__attribute__((always_inline)) INLINE static void cooling_debug_particle(
    const size_t pind) {

  if (pind != -1ull) {
    warning("[PID%lld] cooling_xpart_data:", part_get_id(pind));
    warning("[PID%lld] radiated_energy = %.3e, electron_fraction = %.3e", part_get_id(pind),
            part_get_const_extra_cooling_data(pind)->radiated_energy,
            part_get_const_extra_cooling_data(pind)->electron_fraction);
  }
}

#endif /* SWIFT_COOLING_TREECOOL_DEBUG_H */
