/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifndef SWIFT_GASOLINE_HYDRO_DEBUG_H
#define SWIFT_GASOLINE_HYDRO_DEBUG_H

/**
 * @file Gasoline/hydro_debug.h
 * @brief Density-Energy conservative implementation of SPH,
 *        with added Gasoline physics (Wadsley+ 2017) (Debugging routines)
 */

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {
  warning("[PID%lld] part:",part_get_id(pind));
  warning("[PID%lld] x=[%.3e,%.3e,%.3e]",part_get_id(pind),part_get_x_ind(pind, 0),part_get_x_ind(pind, 1),part_get_x_ind(pind, 2));
  warning("[PID%lld] v=[%.3e,%.3e,%.3e]",part_get_id(pind),part_get_v_ind(pind, 0),part_get_v_ind(pind, 1),part_get_v_ind(pind, 2));
  warning("[PID%lld] a=[%.3e,%.3e,%.3e]",part_get_id(pind),part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),
         part_get_a_hydro_ind(pind, 2));
  warning("[PID%lld] u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e",part_get_id(pind), part_get_u(pind),
          part_get_u_dt(pind), part_get_v_sig(pind), hydro_get_comoving_pressure(pind));
  warning("[PID%lld] h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, dh_drho=%.3e",part_get_id(pind),
          part_get_h(pind), part_get_h_dt(pind), (int)part_get_wcount(pind),part_get_mass(pind),
          part_get_rho_dh(pind));
  warning(
      "[PID%lld] alpha=%.3e, time_bin=%d, rho=%.3e, velocity_gradient=["
      "[%.3e,%.3e,%.3e],"
      "[%.3e,%.3e,%.3e],"
      "[%.3e,%.3e,%.3e]]",
     part_get_id(pind), part_get_alpha_av(pind), part_get_time_bin(pind),part_get_rho(pind),
      part_get_velocity_gradient_ind(pind, 0, 0),
      part_get_velocity_gradient_ind(pind, 0, 1),
      part_get_velocity_gradient_ind(pind, 0, 2),
      part_get_velocity_gradient_ind(pind, 1, 0),
      part_get_velocity_gradient_ind(pind, 1, 1),
      part_get_velocity_gradient_ind(pind, 1, 2),
      part_get_velocity_gradient_ind(pind, 2, 0),
      part_get_velocity_gradient_ind(pind, 2, 1),
      part_get_velocity_gradient_ind(pind, 2, 2));
  warning("[PID%lld] smooth_pressure_gradient=[%.3e,%.3e,%.3e]",part_get_id(pind),
          part_get_smooth_pressure_gradient[0], part_get_smooth_pressure_gradient[1],
          part_get_smooth_pressure_gradient[2]);
  warning("[PID%lld] weighted_wcount=%.3e",part_get_id(pind), p->weighted_wcount);

  warning("[PID%lld] xpart:",part_get_id(pind));
  warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]",part_get_id(pind), part_get_v_full_ind(pind, 0),
          part_get_v_full_ind(pind, 1), part_get_v_full_ind(pind, 2));
}

#endif /* SWIFT_GASOLINE_HYDRO_DEBUG_H */
