/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Abouzied M. A. Nasar (abouzied.nasar@manchester.ac.uk)
 *                    Mladen Ivkovic (mladen.ivkovic@durham.ac.uk)
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
#ifndef GPU_PART_PACK_FUNCTIONS_H
#define GPU_PART_PACK_FUNCTIONS_H

/**
 * @file cuda/gpu_part_pack_functions.h
 * @brief Functions related to packing and unpacking particles to/from a cell.
 * This needs to be specific to any hydro scheme/SPH flavour.
 */

#include "active.h"
#include "cell.h"
#include "gpu_offload_data.h"

/**
 * @brief Unpacks the density data from GPU buffers into cell's particle arrays
 *
 * @param c the #cell
 * @param parts_buffer the "receive" buffer to copy from
 * @param unpack_ind the index in parts_buffer to start copying from
 * @param count number of particles to unpack
 * @param e the #engine
 */
__attribute__((always_inline)) INLINE static void gpu_unpack_part_density(
    struct cell *restrict c,
    const struct gpu_part_recv_d *restrict parts_buffer, const int unpack_ind,
    const size_t count, const struct engine *e) {

  const struct gpu_part_recv_d *parts_recv = &parts_buffer[unpack_ind];
  size_t first = c->hydro.first;

  for (size_t i = 0; i < count; i++) {
    const size_t pind = first + i;

    if (!part_is_active(pind, e)) continue;

    struct gpu_part_recv_d pr = parts_recv[i];

    float rho = part_get_rho(pind) + pr.rho_rhodh_wcount_wcount_dh.x;
    part_set_rho(pind, rho);

    float rho_dh = part_get_rho_dh(pind) + pr.rho_rhodh_wcount_wcount_dh.y;
    part_set_rho_dh(pind, rho_dh);

    float wcount = part_get_wcount(pind) + pr.rho_rhodh_wcount_wcount_dh.z;
    part_set_wcount(pind, wcount);

    float wcount_dh = part_get_wcount_dh(pind) + pr.rho_rhodh_wcount_wcount_dh.w;
    part_set_wcount_dh(pind, wcount_dh);

    part_set_rot_v_ind(pind, 0, part_get_rot_v_ind(pind, 0) + pr.rot_vx_div_v.x);
    part_set_rot_v_ind(pind, 1, part_get_rot_v_ind(pind, 1) + pr.rot_vx_div_v.y);
    part_set_rot_v_ind(pind, 2, part_get_rot_v_ind(pind, 2) + pr.rot_vx_div_v.z);

    float div_v = part_get_div_v(pind) + pr.rot_vx_div_v.w;
    part_set_div_v(pind, div_v);
  }
}

/**
 * @brief Unpacks the gradient data from GPU buffers into cell's particle arrays
 *
 * @param c the #cell
 * @param parts_buffer the "receive" buffer to copy from
 * @param unpack_ind the index in parts_buffer to start copying from
 * @param count number of particles to unpack
 * @param e the #engine
 */
__attribute__((always_inline)) INLINE static void gpu_unpack_part_gradient(
    struct cell *restrict c,
    const struct gpu_part_recv_g *restrict parts_buffer, const int unpack_ind,
    const size_t count, const struct engine *e) {

  const struct gpu_part_recv_g *parts_recv = &parts_buffer[unpack_ind];
  const size_t first = c->hydro.first;

  for (size_t i = 0; i < count; i++) {
    const size_t pind = first + i;

    if (!part_is_active(pind, e)) continue;

    struct gpu_part_recv_g pr = parts_recv[i];

    float avisc_old = part_get_alpha_visc_max_ngb(pind);
    float avisc = fmaxf(avisc_old, pr.aviscmax_vsig_lapu.x);
    part_set_alpha_visc_max_ngb(pind, avisc);

    float vsig = fmaxf(pr.aviscmax_vsig_lapu.y, part_get_v_sig(pind));
    part_set_v_sig(pind, vsig);

    float lu = pr.aviscmax_vsig_lapu.z + part_get_laplace_u(pind);
    part_set_laplace_u(pind, lu);
  }
}

/**
 * @brief Unpacks the force data from GPU buffers into cell's particle arrays
 *
 * @param c the #cell
 * @param parts_buffer the "receive" buffer to copy from
 * @param unpack_ind the index in parts_buffer to start copying from
 * @param count number of particles to unpack
 * @param e the #engine
 */
__attribute__((always_inline)) INLINE static void gpu_unpack_part_force(
    struct cell *restrict c,
    const struct gpu_part_recv_f *restrict parts_buffer, const int unpack_ind,
    const int count, const struct engine *e) {

  const struct gpu_part_recv_f *parts_recv = &parts_buffer[unpack_ind];
  const size_t first = c->hydro.first;

  for (int i = 0; i < count; i++) {
    const size_t pind = first + i;

    if (!part_is_active(pind, e)) continue;

    struct gpu_part_recv_f pr = parts_recv[i];

    part_set_a_hydro_ind(pind, 0, part_get_a_hydro_ind(pind, 0) + pr.a_hydro.x);
    part_set_a_hydro_ind(pind, 1, part_get_a_hydro_ind(pind, 1) + pr.a_hydro.y);
    part_set_a_hydro_ind(pind, 2, part_get_a_hydro_ind(pind, 2) + pr.a_hydro.z);

    float u_dt = pr.udt_hdt_minngbtb.x + part_get_u_dt(pind);
    part_set_u_dt(pind, u_dt);

    float h_dt = pr.udt_hdt_minngbtb.y + part_get_h_dt(pind);
    part_set_h_dt(pind, h_dt);

    timebin_t mintbin = (timebin_t)(pr.udt_hdt_minngbtb.z + 0.5f);
    part_set_timestep_limiter_min_ngb_time_bin(pind, mintbin);
  }
}

/**
 * @brief Packs the cell particle data for density interactions into the
 * CPU-side buffers.
 *
 * @param c the #cell
 * @param parts_buffer the buffer to pack into
 * @param pack_ind the first free index in the buffer arrays to copy data into
 * @param shift periodic boundary shift
 * @param cjstart start index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 * @param cjend end index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 */
__attribute__((always_inline)) INLINE static void gpu_pack_part_density(
    const struct cell *restrict c,
    struct gpu_part_send_d *restrict parts_buffer, const int pack_ind,
    const double shift[3], const int cjstart, const int cjend) {

  /* Grab handles */
  const size_t count = c->hydro.count;
  const size_t first = c->hydro.first;
  struct gpu_part_send_d *ps = &parts_buffer[pack_ind];

  for (size_t i = 0; i < count; i++) {
    const size_t pind = first + i;

    const double x[3] = {part_get_x_ind(pind, 0), part_get_x_ind(pind,1), part_get_x_ind(pind,2)};
    ps[i].x_h.x = x[0] - shift[0];
    ps[i].x_h.y = x[1] - shift[1];
    ps[i].x_h.z = x[2] - shift[2];
    ps[i].x_h.w = part_get_h(pind);

    ps[i].vx_m.x = part_get_v_ind(pind, 0);
    ps[i].vx_m.y = part_get_v_ind(pind, 1);
    ps[i].vx_m.z = part_get_v_ind(pind, 2);
    ps[i].vx_m.w = part_get_mass(pind);

    ps[i].pjs_pje.x = cjstart;
    ps[i].pjs_pje.y = cjend;
  }
}

/**
 * @brief Packs the cell particle data for gradient interactions into
 * the CPU-side buffers.
 *
 * @param ci the #cell
 * @param parts_buffer the buffer to pack into
 * @param pack_ind the first free index in the buffer arrays to copy data into
 * @param shift periodic boundary shift
 * @param cjstart start index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 * @param cjend end index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 */
__attribute__((always_inline)) INLINE static void gpu_pack_part_gradient(
    const struct cell *restrict ci,
    struct gpu_part_send_g *restrict parts_buffer, const int pack_ind,
    const double shift[3], const int cjstart, const int cjend) {

  /* Grab handles */
  const size_t count = ci->hydro.count;
  const size_t first = ci->hydro.first;
  struct gpu_part_send_g *ps = &parts_buffer[pack_ind];

  for (size_t i = 0; i < count; i++) {
    const size_t pind = first + i;

    const double x[3] = {part_get_x_ind(pind, 0), part_get_x_ind(pind,1), part_get_x_ind(pind,2)};
    ps[i].x_h.x = x[0] - shift[0];
    ps[i].x_h.y = x[1] - shift[1];
    ps[i].x_h.z = x[2] - shift[2];
    ps[i].x_h.w = part_get_h(pind);

    ps[i].vx_m.x = part_get_v_ind(pind, 0);
    ps[i].vx_m.y = part_get_v_ind(pind, 1);
    ps[i].vx_m.z = part_get_v_ind(pind, 2);
    ps[i].vx_m.w = part_get_mass(pind);

    ps[i].u_rho_c_aviscmax.x = part_get_u(pind);
    ps[i].u_rho_c_aviscmax.y = part_get_rho(pind);
    ps[i].u_rho_c_aviscmax.z = part_get_soundspeed(pind);
    ps[i].u_rho_c_aviscmax.w = part_get_alpha_visc_max_ngb(pind);

    ps[i].avisc_vsig_lapu.x = part_get_alpha_av(pind);
    ps[i].avisc_vsig_lapu.y = part_get_v_sig(pind);
    ps[i].avisc_vsig_lapu.z = part_get_laplace_u(pind);

    ps[i].pjs_pje.x = cjstart;
    ps[i].pjs_pje.y = cjend;
  }
}

/**
 * @brief Packs the cell particle data for force interactions into the
 * CPU-side buffers.
 *
 * @param ci the #cell
 * @param parts_buffer the buffer to pack into
 * @param pack_ind the first free index in the buffer arrays to copy data into
 * @param shift periodic boundary shift
 * @param cjstart start index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 * @param cjend end index of cell cj's particles (which cell ci is to be
 * interacted with) in buffer
 */
__attribute__((always_inline)) INLINE static void gpu_pack_part_force(
    const struct cell *restrict ci,
    struct gpu_part_send_f *restrict parts_buffer, const int pack_ind,
    const double shift[3], const int cjstart, const int cjend) {

  const size_t count = ci->hydro.count;
  const size_t first = ci->hydro.first;
  struct gpu_part_send_f *ps = &parts_buffer[pack_ind];

  for (size_t i = 0; i < count; i++) {

    const size_t pind = first + i;

    const double x[3] = {part_get_x_ind(pind, 0), part_get_x_ind(pind,1), part_get_x_ind(pind,2)};
    ps[i].x_h.x = x[0] - shift[0];
    ps[i].x_h.y = x[1] - shift[1];
    ps[i].x_h.z = x[2] - shift[2];
    ps[i].x_h.w = part_get_h(pind);

    ps[i].vx_m.x = part_get_v_ind(pind, 0);
    ps[i].vx_m.y = part_get_v_ind(pind, 1);
    ps[i].vx_m.z = part_get_v_ind(pind, 2);
    ps[i].vx_m.w = part_get_mass(pind);

    ps[i].u_rho_f_p.x = part_get_u(pind);
    ps[i].u_rho_f_p.y = part_get_rho(pind);
    ps[i].u_rho_f_p.z = part_get_f_gradh(pind);
    ps[i].u_rho_f_p.w = part_get_pressure(pind);

    ps[i].bals_c_avisc_adiff.x = part_get_balsara(pind);
    ps[i].bals_c_avisc_adiff.y = part_get_soundspeed(pind);
    ps[i].bals_c_avisc_adiff.z = part_get_alpha_av(pind);
    ps[i].bals_c_avisc_adiff.w = part_get_alpha_diff(pind);

    ps[i].timebin_minngbtimebin_pjs_pje.x = (int)part_get_time_bin(pind);
    int mintbin = (int)part_get_timestep_limiter_min_ngb_time_bin(pind);
    ps[i].timebin_minngbtimebin_pjs_pje.y = mintbin;
    ps[i].timebin_minngbtimebin_pjs_pje.z = cjstart;
    ps[i].timebin_minngbtimebin_pjs_pje.w = cjend;
  }
}

#endif /* GPU_PART_PACK_FUNCTIONS_H */
