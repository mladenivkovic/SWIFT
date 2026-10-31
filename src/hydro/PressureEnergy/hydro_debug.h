/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl) &
 *                    Josh Borrow (joshua.borrow@durham.ac.uk)
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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_DEBUG_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_DEBUG_H
/**
 * @file PressureEnergy/hydro_debug.h
 * @brief P-U conservative implementation of SPH (Debugging routines)
 */

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {

  warning("[PID%lld] part:",part_get_id(pind));
  warning(
      "[PID%lld] x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e],a=[%.3e,%.3e,%.3e], "
      "u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e, "
      "h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, rho=%.3e, "
      "p_dh=%.3e, p_bar=%.3e, "
      "time_bin=%d wakeup=%d",
     part_get_id(pind),part_get_x_ind(pind, 0),part_get_x_ind(pind, 1),part_get_x_ind(pind, 2),part_get_v_ind(pind, 0),part_get_v_ind(pind, 1),part_get_v_ind(pind, 2),
     part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),part_get_a_hydro_ind(pind, 2), part_get_u(pind), part_get_u_dt(pind),
      part_get_v_sig(pind), hydro_get_comoving_pressure(pind), part_get_h(pind), part_get_h_dt(pind),
      (int)part_get_wcount(pind),part_get_mass(pind),part_get_rho(pind), part_get_pressure_bar_dh(pind),
      part_get_pressure_bar(pind), part_get_time_bin(pind), part_get_const_limiter_data_p(pind)->wakeup);
  warning("[PID%lld] xpart:",part_get_id(pind));
  warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]",part_get_id(pind), part_get_v_full_ind(pind, 0),
          part_get_v_full_ind(pind, 1), part_get_v_full_ind(pind, 2));
}

#endif /* SWIFT_MINIMAL_HYDRO_DEBUG_H */
