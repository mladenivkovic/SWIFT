/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_TRACERS_FLAMINGO_H
#define SWIFT_TRACERS_FLAMINGO_H

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "black_holes.h"
#include "cooling.h"
#include "part.h"
#include "tracers_struct.h"

/**
 * @brief Update the particle tracers just after it has been initialised at the
 * start of a step.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 */
static INLINE void tracers_after_init(
    const size_t pind, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after it has been drifted.
 *
 * Nothing to do here in the EAGLE model.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 */
static INLINE void tracers_after_drift(
    const size_t pind, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time) {}

/**
 * @brief Update the particle tracers just after its time-step has been
 * computed.
 *
 * In FLAMINGO, we record the highest temperature reached.
 *
 * @param pind Index of the particle.
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 * @param time The current time.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_part)
 */
static INLINE void tracers_after_timestep_part(
    const size_t pind, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling, const double time,
    const double time_step_length, const int *const tracers_triggers_started) {

  struct tracers_xpart_data *restrict td = part_get_extra_tracers_data_p(pind);

  /* Current temperature */
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, pind);
  /* New record? */
  if (temperature > td->maximum_temperature) {

    td->maximum_temperature = temperature;

    if (with_cosmology) {
      td->maximum_temperature_scale_factor = cosmo->a;
    } else {
      td->maximum_temperature_time = time;
    }
  }
}

static INLINE void tracers_after_recording_trigger_part(const int trigger_index, const double time_to_remove) {}

/**
 * @brief Update the star particle tracers just after its time-step has been
 * computed.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param pind Index of the particle.
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_spart)
 */
static INLINE void tracers_after_timestep_spart(
    struct spart *sp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const double time_step_length,
    const int *const tracers_triggers_started) {}

/**
 * @brief Update the black hole particle tracers just after its time-step has
 * been computed.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_bpart)
 */
static INLINE void tracers_after_timestep_bpart(
    struct bpart *bp, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const double time_step_length,
    const int *const tracers_triggers_started) {}

static INLINE void tracers_after_recording_trigger_bpart(
    struct bpart *bp, const int trigger_index, const double time_to_remove) {}

/**
 * @brief Update the sink particle tracers just after its time-step has
 * been computed.
 *
 * Nothing to do here.
 *
 * @param sink Pointer to the sink particle data.
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The current cosmological model.
 * @param time_step_length The length of the step that just finished
 * @param tracers_triggers_started Which triggers have started? (array of size
 * num_snapshot_triggers_sink)
 */
static INLINE void tracers_after_timestep_sink(
    struct sink *sink, const struct unit_system *us,
    const struct phys_const *phys_const, const int with_cosmology,
    const struct cosmology *cosmo, const double time_step_length,
    const int *const tracers_triggers_started) {}

static INLINE void tracers_after_recording_trigger_sink(
    struct sink *sink, const int trigger_index, const double time_to_remove) {}

/**
 * @brief Initialise the tracer data at the start of a calculation.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param hydro_props the hydro_props struct
 * @param cooling The #cooling_function_data used in the run.
 */
static INLINE void tracers_first_init_xpart(
    const size_t pind, const struct unit_system *us,
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct cooling_function_data *cooling) {

  struct tracers_xpart_data *restrict td = part_get_extra_tracers_data_p(pind);

  td->maximum_temperature = -1.f;
  td->maximum_temperature_time = -1.f;
  td->last_AGN_injection_scale_factor = -1.f;
  td->last_AGN_jet_feedback_scale_factor = -1.f;
}

/**
 * @brief Initialise the star tracer data at the start of a calculation.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 */
static INLINE void tracers_first_init_spart(struct spart *sp,
                                            const struct unit_system *us,
                                            const struct phys_const *phys_const,
                                            const struct cosmology *cosmo) {}

/**
 * @brief Initialise the black hole tracer data at the start of a calculation.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param pind Index of the particle.
 * struct).
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 */
static INLINE void tracers_first_init_bpart(struct bpart *bp,
                                            const struct unit_system *us,
                                            const struct phys_const *phys_const,
                                            const struct cosmology *cosmo) {}

/**
 * @brief Initialise the sink tracer data at the start of a calculation.
 *
 * Nothing to do here.
 *
 * @param sink Pointer to the sink particle data.
 * @param us The internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param cosmo The current cosmological model.
 */
static INLINE void tracers_first_init_sink(struct sink *sink,
                                           const struct unit_system *us,
                                           const struct phys_const *phys_const,
                                           const struct cosmology *cosmo) {}

/**
 * @brief Update the gas particle tracer data after it received stellar winds
 * feedback.
 *
 * Nothing to do here.
 *
 * @param delta_p_magnitude Norm of the momentum received (internal physical
 * units).
 * @param delta_u Specific internal energy received (internal physical units).
 * @param kick_velocity Norm of the velocity kick, in the same frame as
 * delta_p_magnitude (internal physical units).
 */
static INLINE void tracers_after_stellar_winds_feedback_part(const float delta_p_magnitude, const float delta_u,
    const float kick_velocity) {}

/**
 * @brief Update the gas particle tracer data after it received supernovae
 * feedback.
 *
 * Nothing to do here.
 *
 * @param delta_p_magnitude Norm of the momentum received (internal physical
 * units).
 * @param delta_u Specific internal energy received (internal physical units).
 * @param kick_velocity Norm of the velocity kick, in the same frame as
 * delta_p_magnitude (internal physical units).
 */
static INLINE void tracers_after_supernovae_feedback_part(const float delta_p_magnitude, const float delta_u,
    const float kick_velocity) {}

/**
 * @brief Update the star particle tracer data after it exploded as type II
 * supernovae.
 *
 * Nothing to do here.
 *
 * @param sp The star particle.
 * @param number_events Number of SNII in this step (may be fractional).
 * @param comoving_density Gas density around the star (comoving).
 * @param with_cosmology Are we running with cosmology?
 * @param cosmo The current cosmological model.
 * @param time The current time (used only without cosmology).
 */
static INLINE void tracers_after_snii_event_spart(struct spart *sp,
                                                  const float number_events,
                                                  const float comoving_density,
                                                  const int with_cosmology,
                                                  const struct cosmology *cosmo,
                                                  const double time) {}

/**
 * @brief Update the star particle tracer data after it exploded as type Ia
 * supernovae.
 *
 * Nothing to do here.
 *
 * @param sp The star particle.
 * @param number_events Number of SNIa in this step (may be fractional).
 * @param comoving_density Gas density around the star (comoving).
 * @param with_cosmology Are we running with cosmology?
 * @param cosmo The current cosmological model.
 * @param time The current time (used only without cosmology).
 */
static INLINE void tracers_after_snia_event_spart(struct spart *sp,
                                                  const float number_events,
                                                  const float comoving_density,
                                                  const int with_cosmology,
                                                  const struct cosmology *cosmo,
                                                  const double time) {}

/**
 * @brief Update the particles' tracer data with values before an AGN feedback
 * event. Note: this function is called in `black_holes_iact.h` before the
 * particle data are updated.
 *
 * @param pind Index of the particle.
 * (internal physical units)
 */
static INLINE void tracers_before_black_holes_feedback(
    const size_t pind, const float scale_factor) {}

/**
 * @brief Update the particles' tracer data after an AGN feedback
 * event.
 *
 * @param with_cosmology Are we running with cosmology?
 * @param scale_factor The current scale-factor (if running with cosmo)
 * @param time The current time (if running without cosmo)
 * @param delta_energy Amount of energy injected in the feedback event
 * (internal physical units)
 */
static INLINE void tracers_after_black_holes_feedback(
    const size_t pind, const int with_cosmology,
    const float scale_factor, const double time, const double delta_energy) {

  struct tracers_xpart_data *restrict td = part_get_extra_tracers_data_p(pind);

  if (with_cosmology)
    td->last_AGN_injection_scale_factor = scale_factor;
  else
    td->last_AGN_injection_time = time;
}

static INLINE void tracers_after_jet_feedback(
    const size_t pind, const int with_cosmology,
    const float scale_factor, const double time, const double delta_energy,
    const float vel_kick, const enum BH_accretion_modes accretion_mode,
    const long long id) {

  struct tracers_xpart_data *restrict td = part_get_extra_tracers_data_p(pind);

  if (with_cosmology)
    td->last_AGN_jet_feedback_scale_factor = scale_factor;
  else
    td->last_AGN_jet_feedback_time = time;
  td->hit_by_jet_feedback++;
  td->jet_feedback_energy += delta_energy;
}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param pind the particle index.
 */
static INLINE void tracers_after_snapshot_part(const size_t pind) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param sp the #spart.
 */
static INLINE void tracers_after_snapshot_spart(struct spart *sp) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * In FLAMINGO, nothing to do.
 *
 * @param bp the #bpart.
 */
static INLINE void tracers_after_snapshot_bpart(struct bpart *bp) {}

/**
 * @brief Tracer event called after a snapshot was written.
 *
 * Nothing to do here.
 *
 * @param sink the #sink.
 */
static INLINE void tracers_after_snapshot_sink(struct sink *sink) {}

/**
 * @brief Split the tracer content of a particle into n pieces
 *
 * @param pind Index of the hydro particle.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void tracers_split_part(
    size_t pind, const double n) {

  struct tracers_xpart_data *restrict td = part_get_extra_tracers_data_p(pind);
  td->jet_feedback_energy /= n;
}

#endif /* SWIFT_TRACERS_FLAMINGO_H */
