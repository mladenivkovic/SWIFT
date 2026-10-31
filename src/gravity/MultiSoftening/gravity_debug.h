/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MULTI_SOFTENING_GRAVITY_DEBUG_H
#define SWIFT_MULTI_SOFTENING_GRAVITY_DEBUG_H

#include "inline.h"

#include <stddef.h>
#include <stdio.h>

__attribute__((always_inline)) INLINE static void gravity_debug_particle(
    const struct gpart *gp) {
  printf(
      "mass=%.3e time_bin=%d\n"
      "x=[%.5e,%.5e,%.5e], v_full=[%.5e,%.5e,%.5e], type=%d, "
      "a=[%.5e,%.5e,%.5e]\n",
      gp->mass, gp->time_bin, gp->x[0], gp->x[1], gp->x[2], gp->v_full[0],
      gp->v_full[1], gp->v_full[2], (int)gp->type, gp->a_grav[0], gp->a_grav[1],
      gp->a_grav[2]);
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  printf("pot=%e\n", gp->potential);
#endif
#ifdef SWIFT_DEBUG_CHECKS
  printf("num_interacted=%lld ti_drift=%lld ti_kick=%lld\n", gp->num_interacted,
         gp->ti_drift, gp->ti_kick);
#endif
}

#endif /* SWIFT_MULTI_SOFTENING_GRAVITY_DEBUG_H */
