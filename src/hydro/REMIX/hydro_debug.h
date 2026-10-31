/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_REMIX_HYDRO_DEBUG_H
#define SWIFT_REMIX_HYDRO_DEBUG_H

/**
 * @file REMIX/hydro_debug.h
 * @brief Debugging routines for REMIX SPH.
 */

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {
  warning("[PID%lld] part:",part_get_id(pind));
  warning(
      "[PID%lld] "
      "x=[%.6g, %.6g, %.6g], v=[%.3g, %.3g, %.3g], "
      "a=[%.3g, %.3g, %.3g], "
      "m=%.3g, u=%.3g, du/dt=%.3g, P=%.3g, c_s=%.3g, "
      "v_sig=%.3g, h=%.3g, dh/dt=%.3g, wcount=%.3g, rho=%.3g, "
      "dh_drho=%.3g, rho_evol=%.3g, time_bin=%d, wakeup=%d, mat_id=%d",
     part_get_id(pind),part_get_x_ind(pind, 0),part_get_x_ind(pind, 1),part_get_x_ind(pind, 2),part_get_v_ind(pind, 0),part_get_v_ind(pind, 1),part_get_v_ind(pind, 2),
     part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),part_get_a_hydro_ind(pind, 2),part_get_mass(pind), part_get_u(pind), part_get_u_dt(pind),
      hydro_get_comoving_pressure(pind), part_get_soundspeed(pind), part_get_v_sig(pind), part_get_h(pind),
      part_get_h_dt(pind), part_get_wcount(pind),part_get_rho(pind), part_get_rho_dh(pind), part_get_rho_evol(pind),
      part_get_time_bin(pind), part_get_const_limiter_data_p(pind)->wakeup, part_get_mat_id(pind));
  warning("[PID%lld] xpart:",part_get_id(pind));
  warning(
      "[PID%lld] "
      "v_full=[%.3g, %.3g, %.3g]",
     part_get_id(pind), part_get_v_full_ind(pind, 0), part_get_v_full_ind(pind, 1), part_get_v_full_ind(pind, 2));
}

#endif /* SWIFT_REMIX_HYDRO_DEBUG_H */
