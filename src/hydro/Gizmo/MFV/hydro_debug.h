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
#ifndef SWIFT_GIZMO_MFV_HYDRO_DEBUG_H
#define SWIFT_GIZMO_MFV_HYDRO_DEBUG_H

#include "error.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {

  const struct fvpm_geometry_struct *restrict geometry =
      part_get_const_fvpm_geometry_p(pind);
  warning(
      "[PID%lld] x=[%.16e,%.16e,%.16e], "
      "v=[%.3e,%.3e,%.3e], "
      "a=[%.3e,%.3e,%.3e], "
      "h=%.3e, "
      "time_bin=%d, "
      "wakeup=%d, "
      "primitives={"
      "v=[%.3e,%.3e,%.3e], "
      "rho=%.3e, "
      "P=%.3e, "
      "gradients={"
      "rho=[%.3e,%.3e,%.3e], "
      "v=[[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e]], "
      "P=[%.3e,%.3e,%.3e]}, "
      "limiter={"
      "rho=[%.3e,%.3e], "
      "v=[[%.3e,%.3e],[%.3e,%.3e],[%.3e,%.3e]], "
      "P=[%.3e,%.3e], "
      "maxr=%.3e}}, "
      "conserved={"
      "momentum=[%.3e,%.3e,%.3e], "
      "mass=%.3e, "
      "energy=%.3e}, "
      "geometry={"
      "volume=%.3e, "
      "matrix_E=[[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e],[%.3e,%.3e,%.3e]]}, "
      "timestepvars={"
      "vmax=%.3e},"
      "density={"
      "wcount_dh=%.3e, "
      "wcount=%.3e}",
     part_get_id(pind),part_get_x_ind(pind, 0),part_get_x_ind(pind, 1),part_get_x_ind(pind, 2),part_get_v_ind(pind, 0),part_get_v_ind(pind, 1),part_get_v_ind(pind, 2),
     part_get_a_hydro_ind(pind, 0),part_get_a_hydro_ind(pind, 1),part_get_a_hydro_ind(pind, 2), part_get_h(pind), part_get_time_bin(pind),
      part_get_const_limiter_data_p(pind)->wakeup, part_get_fluid_v(pind,0), part_get_fluid_v(pind,1), part_get_fluid_v(pind,2),
     part_get_rho(pind), part_get_P(pind), part_get_grad_rho_ind(pind, 0), part_get_grad_rho_ind(pind, 1),
      part_get_grad_rho_ind(pind, 2), part_set_grad_v_ind(pind, 0, 0), part_set_grad_v_ind(pind, 0, 1),
      part_set_grad_v_ind(pind, 0, 2), part_set_grad_v_ind(pind, 1, 0), part_set_grad_v_ind(pind, 1, 1),
      part_set_grad_v_ind(pind, 1, 2), part_set_grad_v_ind(pind, 2, 0), part_set_grad_v_ind(pind, 2, 1),
      part_set_grad_v_ind(pind, 2, 2), part_set_grad_P_ind(pind, 0), part_set_grad_P_ind(pind, 1),
      part_set_grad_P_ind(pind, 2), part_get_limiter_rho(pind, 0), part_get_limiter_rho(pind, 1),
      part_get_limiter_v(pind, 0, 0), part_get_limiter_v(pind, 0, 1), part_get_limiter_v(pind, 1, 0),
      part_get_limiter_v(pind, 1, 1), part_get_limiter_v(pind, 2, 0), part_get_limiter_v(pind, 2, 1),
      part_set_limiter_P_ind(pind, 0), part_set_limiter_P_ind(pind, 1), part_get_limiter_maxr(pind),
      part_set_conserved_momentum(pind, 0), part_set_conserved_momentum(pind, 1),
      part_set_conserved_momentum(pind, 2), part_get_conserved_mass(pind), part_get_conserved_energy(pind),
      part_get_geometry_volume(pind),
      geometry->matrix_E[0][0],
      geometry->matrix_E[0][1], geometry->matrix_E[0][2],
      geometry->matrix_E[1][0], geometry->matrix_E[1][1],
      geometry->matrix_E[1][2], geometry->matrix_E[2][0],
      geometry->matrix_E[2][1], geometry->matrix_E[2][2],
      part_get_timestepvars_vmax(pind), part_get_rho_dh(pind), part_get_wcount(pind));
}

#endif /* SWIFT_GIZMO_MFV_HYDRO_DEBUG_H */
