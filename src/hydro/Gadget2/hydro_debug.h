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
#ifndef SWIFT_GADGET2_HYDRO_DEBUG_H
#define SWIFT_GADGET2_HYDRO_DEBUG_H

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
      "[PID%lld] a=[%.3e,%.3e,%.3e] h=%.3e, wcount=%.3f, wcount_dh=%.3e, "
      "m=%.3e, dh_drho=%.3e, rho=%.3e, "
      "P=%.3e, P_over_rho2=%.3e, S=%.3e, dS/dt=%.3e, c=%.3e",
     part_get_id(pind),part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),part_get_a_hydro_ind(pind, 2), part_get_h(pind),
      part_get_wcount(pind), part_get_rho_dh(pind),part_get_mass(pind), part_get_rho_dh(pind),
     part_get_rho(pind), hydro_get_comoving_pressure(pind), part_get_p_over_rho2(pind), part_get_entropy(pind),
      part_get_entropy_dt(pind), part_get_soundspeed(pind));
  warning("[PID%lld] divV=%.3e, rotV=[%.3e,%.3e,%.3e], balsara=%.3e",part_get_id(pind),
          part_get_div_v(pind), part_get_rot_v(pind, 0), part_get_rot_v(pind, 1),
          part_get_rot_v(pind, 2), part_get_balsara(pind));
  warning("[PID%lld] v_sig=%e dh/dt=%.3e time_bin=%d wakeup=%d",part_get_id(pind),
          p->force.v_sig, part_get_h_dt(pind), part_get_time_bin(pind), part_get_const_limiter_data_p(pind)->wakeup);
  warning("[PID%lld] xpart:",part_get_id(pind));
  warning("[PID%lld] v_full=[%.3e,%.3e,%.3e]",part_get_id(pind), part_get_v_full_ind(pind, 0),
          part_get_v_full_ind(pind, 1), part_get_v_full_ind(pind, 2));
}

#endif /* SWIFT_GADGET2_HYDRO_DEBUG_H */
