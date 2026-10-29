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
#ifndef SWIFT_COOLING_GRACKLE_DEBUG_H
#define SWIFT_COOLING_GRACKLE_DEBUG_H

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void cooling_debug_particle(
    const size_t pind) {

  const struct xpart_cooling_data *restrict cd = part_get_const_extra_cooling_data_p(pind);

  if (pind != -1ull) {
    warning("[PID%lld] cooling_xpart_data:",part_get_id(pind));
    warning("[PID%lld] radiated_energy = %.3e, time_last_event = %.3e",part_get_id(pind),
            cd->radiated_energy, cd->time_last_event);

#if COOLING_GRACKLE_MODE >= 1
    warning(
        "[PID%lld] HI_frac = %.3e, HII_frac = %.3e, HeI_frac = %.3e, HeII_frac "
        "= "
        "%.3e, "
        "HeIII_frac = %.3e, e_frac = %.3e",
       part_get_id(pind), cd->HI_frac, cd->HII_frac,
        cd->HeI_frac, cd->HeII_frac,
        cd->HeIII_frac, cd->e_frac);
#if COOLING_GRACKLE_MODE >= 2
    warning("[PID%lld] HM_frac = %.3e, H2I_frac = %.3e, H2II_frac = %.3e",
           part_get_id(pind), cd->HM_frac, cd->H2I_frac,
            cd->H2II_frac);
#if COOLING_GRACKLE_MODE >= 3
    warning("[PID%lld] DI_frac = %.3e, DII_frac = %.3e, HDI_frac = %.3e",part_get_id(pind),
            cd->DI_frac, cd->DII_frac,
            cd->HDI_frac);
#endif
#endif
#endif
    warning("[PID%lld] metal_frac = %.3e",part_get_id(pind), cd->metal_frac);
  }
}

#endif /* SWIFT_COOLING_GRACKLE_DEBUG_H */
