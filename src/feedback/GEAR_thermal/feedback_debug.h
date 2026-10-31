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
#ifndef SWIFT_FEEDBACK_GEAR_DEBUG_H
#define SWIFT_FEEDBACK_GEAR_DEBUG_H

#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void feedback_debug_particle(
    const size_t pind) {

  const struct feedback_xpart_data *restrict fd = part_get_const_extra_feedback_data_p(pind);

  warning("[PID%lld] feedback_xpart_data:", part_get_id(pind));
  warning(
      "[PID%lld] delta_mass = %.3e, delta_u = %.3e, delta_p = [%.3e, %.3e, %.3e]",
      part_get_id(pind), fd->delta_mass, fd->delta_u, fd->delta_p[0], fd->delta_p[1], fd->delta_p[2]);
}

#endif /* SWIFT_FEEDBACK_GEAR_DEBUG_H */
