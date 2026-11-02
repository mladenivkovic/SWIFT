/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2026 Filip Husko (husko@strw.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_IDEALIZED_AGN_JET_H
#define SWIFT_FORCING_IDEALIZED_AGN_JET_H

/* Config parameters. */
#include <config.h>

/* Standard includes. */
#include <float.h>

/* Local includes. */
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "space.h"
#include "timestep_sync_part.h"
#include "units.h"

/**
 * @brief Forcing Term Properties
 */
struct forcing_terms {

  /*! Particles with ID <= boundary_particle_max_id are treated as boundary
   * particles. */
  int boundary_particle_max_id;

  /*! Whether or not boundary particles are fixed in place (0 (false, default)
   * or 1 (true))  */
  int enable_fixed_position;

  /*! Whether or not boundary particles have hydrodynamic accelerations enabled
   * (0 (false, default) or 1 (true))  */
  int enable_hydro_acceleration;

  /*! Whether or not boundary particles have gravitational accelerations enabled
   * (0 (false, default) or 1 (true))  */
  int enable_grav_acceleration;

  /* The jet power */
  float jet_power;

  /* The jet kick velocity */
  float jet_velocity;

  /* How long the jet should be on */
  float jet_duration;
};

/**
 * @brief Computes the hydrodynamic forcing terms.
 *
 * Resets the hydrodynamic accelerations of boundary particles to zero.
 * If using fixed boundary particles, their velocities are also set to zero.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param s The #space we act on.
 * @param phys_const The physical constants in internal units.
 * @param pind Index of the particle.
 */
__attribute__((always_inline)) INLINE static void forcing_hydro_terms_apply(
    const double time, const struct forcing_terms *terms, const struct space *s,
    const struct phys_const *phys_const, size_t pind) {

  if (terms->enable_hydro_acceleration) {
    /* Skip if not resetting hydro accelerations. */
    return;
  }

  /* Only kick particles if they haven't yet been kicked, and if it is time to
     kick them; define max_id of particle that should be kicked up to now */
  const float max_id =
      (terms->jet_power * time) /
      (0.5 *part_get_mass(pind) * terms->jet_velocity * terms->jet_velocity);

  if ((part_get_id(pind) <= terms->boundary_particle_max_id) && (part_get_id(pind) > max_id)) {
    /* Reset the hydro accelerations of boundary particles. */
   part_set_a_hydro_ind(pind, 0, 0.0f);
   part_set_a_hydro_ind(pind, 1, 0.0f);
   part_set_a_hydro_ind(pind, 2, 0.0f);

#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)

    /* Some values need to be reset in the Gizmo case. */
    hydro_prepare_force(p, &c->hydro.xparts[k], cosmo, e->hydro_properties,
                        e->pressure_floor_props,
                        /*dt_alpha=*/0, /*dt_therm=*/0);
    rt_prepare_force(p);
#endif

    if (terms->enable_fixed_position) {
      /* Set velocity of fixed boundary particle to zero. */
      part_set_v_full_ind(pind, 0, 0.f);
      part_set_v_full_ind(pind, 1, 0.f);
      part_set_v_full_ind(pind, 2, 0.f);
    }
  } else if (part_get_id(pind) <= terms->boundary_particle_max_id) {
    /* Launch particles into jets */

    /* Particle position relative to the box centre */
    const double delta_x = (part_get_x_ind(pind, 0) - s->dim[0] / 2.0);
    const double delta_y = (part_get_x_ind(pind, 1) - s->dim[1] / 2.0);
    const double delta_z = (part_get_x_ind(pind, 2) - s->dim[2] / 2.0);

    const double R =
        sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

    /* Assign velocity to be given. We do a kick radially from the origin */
    const double vel_kick_vec[3] = {terms->jet_velocity * delta_x / R,
                                    terms->jet_velocity * delta_y / R,
                                    terms->jet_velocity * delta_z / R};

    part_set_v_ind(pind, 0, vel_kick_vec[0]);
    part_set_v_ind(pind, 1, vel_kick_vec[1]);
    part_set_v_ind(pind, 2, vel_kick_vec[2]);
    part_set_v_full_ind(pind, 0, vel_kick_vec[0]);
    part_set_v_full_ind(pind, 1, vel_kick_vec[1]);
    part_set_v_full_ind(pind, 2, vel_kick_vec[2]);

    /* Reset some hydro quantities */
    hydro_diffusive_feedback_reset(pind);

    /* Recompute the signal velocity of the particle */
    hydro_set_v_sig_based_on_velocity_kick(pind, s->e->cosmology,
                                           terms->jet_velocity);

    /* Synchronize the particle on the time-line */
    timestep_sync_part(pind);

    /* Increase the particle's id so it's no longer ever kicked */
    part_set_id(pind, part_get_id(pind) + terms->boundary_particle_max_id);
  }
}

/**
 * @brief Computes the gravitational forcing terms.
 *
 * Resets the gravitational accelerations of boundary particles to zero.
 * If using fixed boundary particles, their velocities are also set to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_grav_terms_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {
  /* Nothing to do here */
}

/**
 * @brief Sets the forcing of gparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param gp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_gpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct gpart *gp) {
  /* Nothing to do here */
}

/**
 * @brief Sets the forcing of parts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param pind Index of the particle.
 */
__attribute__((always_inline)) INLINE static void forcing_part_drift_apply(
    const long long id, const struct forcing_terms *terms, size_t pind) {

  if (id <= terms->boundary_particle_max_id && terms->enable_fixed_position) {
    /* Set velocity of fixed boundary particle to zero. */
    part_set_v_full_ind(pind, 0, 0.f);
    part_set_v_full_ind(pind, 1, 0.f);
    part_set_v_full_ind(pind, 2, 0.f);
  }
}

/**
 * @brief Sets the forcing of sparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_spart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct spart *sp) {
  /* Nothing to do here */
}

/**
 * @brief Sets the forcing of bparts prior to drift.
 *
 * Resets the velocities of fixed boundary particles to zero.
 *
 * @param id The particle ID.
 * @param terms The properties of the forcing terms.
 * @param bp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void forcing_bpart_drift_apply(
    const long long id, const struct forcing_terms *terms, struct bpart *bp) {
  /* Nothing to do here */
}

/**
 * @brief Computes the time-step condition due to the forcing terms.
 *
 * Nothing to do here. --> Return FLT_MAX.
 *
 * @param time The current time.
 * @param terms The properties of the forcing terms.
 * @param phys_const The physical constants in internal units.
 * @param pind Index of the particle.
 */
__attribute__((always_inline)) INLINE static float forcing_terms_timestep(
    double time, const struct forcing_terms *terms, const struct space *s,
    const struct phys_const *phys_const, const size_t pind) {

  return FLT_MAX;
}

/**
 * @brief updates the forcing terms
 *
 * Nothing to do here
 *
 * @param terms The #forcing_terms properties of the run
 * @param time_old The previous system time
 */
INLINE static void forcing_update(struct forcing_terms *terms,
                                  const double time_old) {}

/**
 * @brief Prints the properties of the forcing terms to stdout.
 *
 * @param terms The #forcing_terms properties of the run.
 */
static INLINE void forcing_terms_print(const struct forcing_terms *terms) {

  message("Forcing terms is 'Idealized AGN jet'. Max boundary particle ID: %i.",
          terms->boundary_particle_max_id);
  message(
      "Run with options: Enable fixed position: %i, Enable hydro acceleration: "
      "%i, Enable grav acceleration: %i",
      terms->enable_fixed_position, terms->enable_hydro_acceleration,
      terms->enable_grav_acceleration);
}

/**
 * @brief Initialises the forcing term properties
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space object.
 * @param terms The forcing term properties to initialize
 */
static INLINE void forcing_terms_init(struct swift_params *parameter_file,
                                      const struct phys_const *phys_const,
                                      const struct unit_system *us,
                                      const struct space *s,
                                      struct forcing_terms *terms) {

  terms->boundary_particle_max_id = parser_get_param_int(
      parameter_file, "BoundaryParticles:boundary_particle_max_id");
  terms->enable_fixed_position = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_fixed_position", 0);
  terms->enable_hydro_acceleration = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_hydro_acceleration", 0);
  terms->enable_grav_acceleration = parser_get_opt_param_int(
      parameter_file, "BoundaryParticles:enable_grav_acceleration", 0);

  if ((terms->enable_fixed_position != 0) &&
      (terms->enable_fixed_position != 1)) {
    error(
        "BoundaryParticles:enable_fixed_position must be either 0 (false) or 1 "
        "(true).");
  }

  if (terms->enable_fixed_position) {
    /* If using fixed boundary particles, both hydro and grav accelerations must
     * be reset to 0. */
    if (terms->enable_hydro_acceleration) {
      error(
          "BoundaryParticles:enable_hydro_acceleration must be 0 (false) when "
          "using fixed boundary particles.");
    }
    if (terms->enable_grav_acceleration) {
      error(
          "BoundaryParticles:enable_grav_acceleration must be 0 (false) when "
          "using fixed boundary particles.");
    }

  } else {
    /* If not using fixed boundary particles, hydro and grav accelerations can
     * be enabled independently. */
    if ((terms->enable_hydro_acceleration != 0) &&
        (terms->enable_hydro_acceleration != 1)) {
      error(
          "BoundaryParticles:enable_hydro_acceleration must be either 0 "
          "(false) or 1 (true).");
    }
    if ((terms->enable_grav_acceleration != 0) &&
        (terms->enable_grav_acceleration != 1)) {
      error(
          "BoundaryParticles:enable_grav_acceleration must be either 0 (false) "
          "or 1 (true).");
    }
  }

  terms->jet_power =
      parser_get_param_float(parameter_file, "BoundaryParticles:jet_power");
  terms->jet_velocity =
      parser_get_param_float(parameter_file, "BoundaryParticles:jet_velocity");
  terms->jet_duration =
      parser_get_param_float(parameter_file, "BoundaryParticles:jet_duration");
}

/**
 * @brief Clean-up the memory allocated for the forcing routine
 *
 * Nothing to do here
 *
 * @param terms The forcing term properties
 */
static INLINE void forcing_terms_clean(struct forcing_terms *terms) {}

#endif /* SWIFT_FORCING_IDEALIZED_AGN_JET_H */
