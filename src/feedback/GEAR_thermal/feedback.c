/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "feedback.h"

/* Local includes */
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

#include <strings.h>

/**
 * @brief Update the properties of the particle due to a supernovae.
 *
 * @param pind Index of the hydro particle to consider.
 * @param e The #engine.
 */
void feedback_update_part(size_t pind, const struct engine *e) {

  struct feedback_xpart_data *restrict fd = part_get_extra_feedback_data_p(pind);

  /* Did the particle receive a supernovae */
  if (!fd->hit_by_SN && !fd->hit_by_winds) return;

  const struct cosmology *cosmo = e->cosmology;
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;

  /* Turn off the cooling only when SN are involved */
  if (fd->hit_by_SN) {
    cooling_set_part_time_cooling_off(pind, e->time);
  }

  /* Update mass */
  const float old_mass = hydro_get_mass(pind);
  const float new_mass = old_mass + fd->delta_mass;

  if (fd->delta_mass < 0.) {
    error("Delta mass smaller than 0");
  }

  hydro_set_mass(pind, new_mass);

  fd->delta_mass = 0;

  /* Update the density */
  part_set_rho(pind, part_get_rho(pind) * new_mass / old_mass);

  /* Update internal energy */
  const float u = hydro_get_physical_internal_energy(pind, cosmo) * old_mass / new_mass;
  const float u_new = u + fd->delta_u;

  hydro_set_physical_internal_energy(pind, cosmo, u_new);
  hydro_set_drifted_physical_internal_energy(pind, cosmo, pressure_floor, u_new);

  fd->delta_u = 0.;

  /* Update the velocities */
  for (int i = 0; i < 3; i++) {
    const float dv = fd->delta_p[i] / new_mass;

    part_set_v_full_ind(pind, i, part_get_v_full_ind(pind, i) + dv);
    part_set_v_ind(pind, i, part_get_v_ind(pind, i) + dv);

    fd->delta_p[i] = 0;
  }

  fd->hit_by_SN = 0;
  fd->hit_by_winds = 0;
}

/**
 * @brief Finishes the #part density calculation.
 *
 * Nothing to do here.
 *
 * @param pind Index of the particle to act upon
 */
__attribute__((always_inline)) INLINE void feedback_end_density(size_t pind) {}

/**
 * @brief Reset the gas particle-carried fields related to feedback at the
 * start of a step.
 *
 * Nothing to do here in the GEAR model.
 *
 * @param pind Index of the particle.
 */
void feedback_reset_part(size_t pind) {}

/**
 * @brief Should this particle be doing any feedback-related operation?
 *
 * @param sp The #spart.
 * @param e The #engine.
 */
int feedback_is_active(const struct spart *sp, const struct engine *e) {

  /* the particle is inactive if its birth_scale_factor or birth_time is
   * negative */
  if (sp->birth_scale_factor < 0.0 || sp->birth_time < 0.0) return 0;

  return sp->feedback_data.will_do_feedback;
}

/**
 * @brief Prepares a s-particle for its feedback interactions
 *
 * @param sp The particle to act upon
 */
void feedback_init_spart(struct spart *sp) {

  sp->feedback_data.enrichment_weight = 0.f;
}

/**
 * @brief Prepares a star's feedback field before computing what
 * needs to be distributed.
 *
 * This is called in the stars ghost.
 */
void feedback_reset_feedback(struct spart *sp,
                             const struct feedback_props *feedback_props) {}

/**
 * @brief Initialises the s-particles feedback props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param feedback_props The properties of the feedback model.
 */
void feedback_prepare_spart(struct spart *sp,
                            const struct feedback_props *feedback_props) {}

/**
 * @brief Prepare a #spart for the feedback task.
 *
 * This is called in the stars ghost task.
 *
 * In here, we only need to add the missing coefficients.
 *
 * @param sp The particle to act upon
 * @param feedback_props The #feedback_props structure.
 * @param cosmo The current cosmological model.
 * @param us The unit system.
 * @param phys_const The #phys_const.
 * @param star_age_beg_step The age of the star at the star of the time-step in
 * internal units.
 * @param dt The time-step size of this star in internal units.
 * @param time The physical time in internal units.
 * @param ti_begin The integer time at the beginning of the step.
 * @param with_cosmology Are we running with cosmology on?
 */
void feedback_prepare_feedback(struct spart *restrict sp,
                               const struct feedback_props *feedback_props,
                               const struct cosmology *cosmo,
                               const struct unit_system *us,
                               const struct phys_const *phys_const,
                               const double star_age_beg_step, const double dt,
                               const double time, const integertime_t ti_begin,
                               const int with_cosmology) {
  /* Add missing h factor */
  const float hi_inv = 1.f / sp->h;
  const float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  sp->feedback_data.enrichment_weight *= hi_inv_dim;
}
