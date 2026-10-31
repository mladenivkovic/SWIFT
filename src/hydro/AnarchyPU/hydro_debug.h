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
#ifndef SWIFT_ANARCHY_PU_HYDRO_DEBUG_H
#define SWIFT_ANARCHY_PU_HYDRO_DEBUG_H
/**
 * @file AnarchyPU/hydro_debug.h
 * @brief P-U conservative implementation of SPH,
 *        with added ANARCHY physics (Cullen & Denhen 2011 AV,
 *        Price 2008 thermal diffusion (Debugging routines).
 */

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {
  warning("[PID%lld] part:",part_get_id(pind));
  warning(
      "[PID%lld] x=[%.3e,%.3e,%.3e], "
      "v=[%.3e,%.3e,%.3e]",
     part_get_id(pind),part_get_x_ind(pind, 0),part_get_x_ind(pind, 1),part_get_x_ind(pind, 2),part_get_v_ind(pind, 0),part_get_v_ind(pind, 1),part_get_v_ind(pind, 2));
  warning(
      "[PID%lld] a=[%.3e,%.3e,%.3e], "
      "u=%.3e, du/dt=%.3e v_sig=%.3e, P=%.3e",
     part_get_id(pind),part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),part_get_a_hydro_ind(pind, 2), part_get_u(pind), part_get_u_dt(pind),
      part_get_v_sig(pind), hydro_get_comoving_pressure(pind));
  warning(
      "[PID%lld] h=%.3e, dh/dt=%.3e wcount=%d, m=%.3e, dh_drho=%.3e, rho=%.3e",
     part_get_id(pind), part_get_h(pind), part_get_h_dt(pind), (int)part_get_wcount(pind),part_get_mass(pind),
      part_get_rho_dh(pind),part_get_rho(pind));
  warning("[PID%lld] p_dh=%.3e, p_bar=%.3e, alpha=%.3e",part_get_id(pind),
          part_get_pressure_bar_dh(pind), part_get_pressure_bar(pind), part_get_alpha_av(pind));
  warning("[PID%lld] time_bin=%d",part_get_id(pind), part_get_time_bin(pind));
  warning("[PID%lld] xpart:",part_get_id(pind));
  warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]",part_get_id(pind), part_get_v_full_ind(pind, 0),
          part_get_v_full_ind(pind, 1), part_get_v_full_ind(pind, 2));
}

#endif /* SWIFT_ANARCHY_PU_HYDRO_DEBUG_H */
