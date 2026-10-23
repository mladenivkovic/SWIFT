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
#ifndef SWIFT_BLACK_HOLES_EAGLE_DEBUG_H
#define SWIFT_BLACK_HOLES_EAGLE_DEBUG_H

#include "inline.h"
#include "error.h"
#include "hydro_part.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void black_holes_debug_particle(
    const size_t pind) {

  warning("[PID%lld] black_holes_part_data:", part_get_id(pind));
  warning("[PID%lld] swallow_id = %lld, potential = %.3e", part_get_id(pind),
          part_get_black_holes_data_p(pind)->swallow_id, part_get_black_holes_data_p(pind)->potential);
}

#endif /* SWIFT_BLACK_HOLES_EAGLE_DEBUG_H */
