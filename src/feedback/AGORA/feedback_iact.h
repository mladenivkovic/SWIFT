/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yves Revaz (yves.reavz@epfl.ch)
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
#ifndef SWIFT_AGORA_FEEDBACK_IACT_H
#define SWIFT_AGORA_FEEDBACK_IACT_H

/* Local includes */
#include "chemistry.h"
#include "feedback.h"
#include "hydro.h"
#include "random.h"
#include "timestep_sync_part.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pindj Index of second particle (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const size_t pindj,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {

  /* Get the gas mass. */
  const float mj = hydro_get_mass(pindj);

  /* Get r */
  const float r = sqrtf(r2);

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  float wi;
  kernel_eval(ui, &wi);

  /* Add contribution of pj to normalisation of density weighted fraction
   * which determines how much mass to distribute to neighbouring
   * gas particles */

  /* The normalization by 1 / h^d is done in feedback.h */
  si->feedback_data.enrichment_weight += mj * wi;
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, size_t pindj, const struct cosmology *cosmo, const struct hydro_props *hydro_props, const struct feedback_props *fb_props, const integertime_t ti_current) {

  const double e_sn = si->feedback_data.energy_ejected;

  /* Do we have supernovae? */
  if (e_sn == 0) {
    return;
  }

  const float mj = hydro_get_mass(pindj);
  const float r = sqrtf(r2);

  /* Get the kernel for hi. */
  float hi_inv = 1.0f / hi;
  float hi_inv_dim = pow_dimension(hi_inv); /* 1/h^d */
  float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  wi *= hi_inv_dim;

  /* Compute inverse enrichment weight */
  const double si_inv_weight = si->feedback_data.enrichment_weight == 0
                                   ? 0.
                                   : 1. / si->feedback_data.enrichment_weight;

  /* Mass received */
  const double m_ej = si->feedback_data.mass_ejected;
  const double weight = mj * wi * si_inv_weight;
  const double dm = m_ej * weight;
  const double new_mass = mj + dm;

  /* Energy received */
  const double du = e_sn * weight / new_mass;

  struct feedback_xpart_data *restrict fdj = part_get_extra_feedback_data_p(pindj);
  fdj->delta_mass += dm;
  fdj->delta_u += du;

  /* Compute momentum received. */
  for (int i = 0; i < 3; i++) {
    fdj->delta_p[i] += dm * (si->v[i] - part_get_v_full_ind(pindj, i));
  }

  /* Add the metals */
  for (int i = 0; i < AGORA_CHEMISTRY_ELEMENT_COUNT; i++) {
   part_get_chemistry_data_p(pindj)->metal_mass[i] +=
        weight * si->feedback_data.metal_mass_ejected[i];
  }

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pindj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pindj);
}

#endif /* SWIFT_AGORA_FEEDBACK_IACT_H */
