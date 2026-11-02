/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_STAR_FORMATION_H
#define SWIFT_GEAR_STAR_FORMATION_H

/* Local includes */
#include "chemistry.h"
#include "cooling.h"
#include "cosmology.h"
#include "engine.h"
#include "entropy_floor.h"
#include "error.h"
#include "feedback.h"
#include "hydro_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "random.h"
#include "star_formation_setters.h"
#include "star_formation_struct.h"
#include "stars.h"
#include "units.h"

/* Does the star formation model move the hydro particles?
   This will update the c->hydro.dx_max_part and
   c->hydro.dx_max_sort after forming a star. */
#define star_formation_need_update_dx_max 1

/**
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * Use the star formation criterion given by eq. 3 in Revaz & Jablonka 2018.
 *
 * @param pind Index of the particle.
 * @param starform the star formation law properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 * @param entropy_floor The #entropy_floor_properties.
 *
 */
INLINE static int star_formation_is_star_forming(
    size_t pind,
    const struct star_formation *starform, const struct phys_const *phys_const,
    const struct cosmology *cosmo,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cooling_function_data *restrict cooling,
    const struct entropy_floor_properties *restrict entropy_floor) {

  /* Check if collapsing particles */
  if (part_get_const_sf_data_p(pind)->div_v > 0) {
    return 0;
  }

  const float temperature = cooling_get_temperature(phys_const, hydro_props, us, cosmo, cooling, pind);

  const float density = hydro_get_physical_density(pind, cosmo);

  const float temperature_max = starform->maximal_temperature;

  const float density_threshold = starform->density_threshold;

  /* Check the temperature criterion */
  if (temperature > temperature_max) {
    return 0;
  }

  /* Check the density threshold */
  if (density < density_threshold) {
    return 0;
  }

  /* mode = 4 : density > density_threshold : accept */
  if (starform->star_formation_mode == gear_star_formation_agora) return 1;

  /* Get the required variables */
  const float n_jeans_2_3 = starform->n_jeans_2_3;

  const float h = part_get_h(pind) * kernel_gamma * cosmo->a;

  // TODO use GRACKLE */
  const float mu = hydro_props->mu_neutral;

  /* Compute the density criterion */
  const float coef =
      M_PI_4 / (phys_const->const_newton_G * n_jeans_2_3 * h * h);
  const float density_criterion = coef * hydro_gamma *
                                  phys_const->const_boltzmann_k * temperature /
                                  (mu * phys_const->const_proton_mass);

  /* Check the density criterion */
  return density > density_criterion;
}

/**
 * @brief Compute the star-formation rate of a given particle.
 *
 * Nothing to do here. Everything is done in
 * #star_formation_should_convert_to_star.
 *
 * @param pind Index of the particle.
 * @param starform the star formation law properties to use
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cosmo the cosmological parameters and properties.
 * @param dt_star The time-step of this particle.
 */
INLINE static void star_formation_compute_SFR(
    size_t pind, const struct star_formation *starform, const struct phys_const *phys_const,
    const struct hydro_props *hydro_props, const struct cosmology *cosmo,
    const double dt_star) {}

/**
 * @brief Decides whether a particle should be converted into a
 * star or not.

 * Compute the star formation rate from eq. 4 in Revaz & Jablonka 2012.
 *
 * @param pind Index of the hydro particle.
 * @param starform The properties of the star formation model.
 * @param e The #engine (for random numbers).
 * @param dt_star The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int star_formation_should_convert_to_star(
    size_t pind, const struct star_formation *starform,
    const struct engine *e, const double dt_star) {

  const struct phys_const *phys_const = e->physical_constants;
  const struct cosmology *cosmo = e->cosmology;

  /* Check that we are running a full time step */
  if (dt_star == 0.) {
    return 0;
  }

  /* Get a few variables */
  const float G = phys_const->const_newton_G;
  const float density = hydro_get_physical_density(pind, cosmo);

  /* Get the mass of the future possible star */
  const float mass_gas = hydro_get_mass(pind);

  /* Compute the probability */
  const float inv_free_fall_time =
      sqrtf(density * 32.f * G * 0.33333333f * M_1_PI);
  float prob = 1.f - expf(-starform->star_formation_efficiency *
                          inv_free_fall_time * dt_star);

  /* Add the mass factor */
  if (starform->n_stars_per_part != 1) {
    const float min_mass =
        starform->mass_stars * starform->min_mass_frac_plus_one;
    const float mass_star =
        mass_gas > min_mass ? starform->mass_stars : mass_gas;

    prob *= mass_gas / mass_star;
  }

  /* Roll the dice... */
  const float random_number =
      random_unit_interval(part_get_id(pind), e->ti_current, random_number_star_formation);

  /* Can we form a star? */
  return random_number < prob;
}

/**
 * @brief Returns the number of new star particles to create per SF event.
 *
 * @param pind Index of the hydro particle.
 * @param starform The properties of the star formation model.
 *
 * @return The number of extra star particles to generate per gas particles.
 *        (return 0 if the gas particle itself is to be converted)
 */
INLINE static int star_formation_number_spart_to_spawn(
    size_t pind, const struct star_formation *starform) {

  /* Check if we are splitting the particles or not */
  if (starform->n_stars_per_part == 1) {
    return 0;
  }

  const float mass_min =
      starform->min_mass_frac_plus_one * starform->mass_stars;
  return hydro_get_mass(pind) > mass_min;
}

/**
 * @brief Returns the number of particles to convert per SF event.
 *
 * @param pind Index of the hydro particle.
 * @param starform The properties of the star formation model.
 *
 * @return The number of particles to generate per gas particles.
 *        (This has to be 0 or 1)
 */
INLINE static int star_formation_number_spart_to_convert(
    const size_t pind,
    const struct star_formation *starform) {

  if (starform->n_stars_per_part == 1) {
    return 1;
  }

  const float mass_min =
      starform->min_mass_frac_plus_one * starform->mass_stars;
  return hydro_get_mass(pind) <= mass_min;
}

/**
 * @brief Update the SF properties of a particle that is not star forming.
 *
 * @param pind Index of the hydro particle.
 * @param e The #engine.
 * @param starform The properties of the star formation model.
 * @param with_cosmology Are we running with cosmology switched on?
 */
INLINE static void star_formation_update_part_not_SFR(
    size_t pind, const struct engine *e,
    const struct star_formation *starform, const int with_cosmology) {}

/**
 * @brief Separate the #spart and #part by randomly moving both of them.
 *
 * @param e The #engine.
 * @param pind Index of the hydro particle generating a star.
 * @param sp The new #spart.
 * @param (return) displacement The 3D displacement vector of the star with
 * respect to the star position.
 */
INLINE static void star_formation_separate_particles(const struct engine *e,
                                                     size_t pind,
                                                     struct spart *sp,
                                                     float displacement[3]) {

  const double x[3] = {part_get_x_ind(pind, 0), part_get_x_ind(pind, 1), part_get_x_ind(pind, 2)};
  const float h = part_get_h(pind);
  const long long id = part_get_id(pind);

#ifdef SWIFT_DEBUG_CHECKS
  if (x[0] != sp->x[0] || x[1] != sp->x[1] || x[2] != sp->x[2]) {
    error(
        "Moving particles that are not at the same location."
        " (%g, %g, %g) - (%g, %g, %g)",
        x[0], x[1], x[2], sp->x[0], sp->x[1], sp->x[2]);
  }
#endif

  /* Move a bit the particle in order to avoid
     division by 0.
  */
  const float max_displacement = 0.1;
  const double delta_x =
      2.f * random_unit_interval(id, e->ti_current, (enum random_number_type)0) - 1.f;
  const double delta_y =
      2.f * random_unit_interval(id, e->ti_current, (enum random_number_type)1) - 1.f;
  const double delta_z =
      2.f * random_unit_interval(id, e->ti_current, (enum random_number_type)2) - 1.f;

  /* Update the displacement */
  displacement[0] = delta_x * max_displacement * h;
  displacement[1] = delta_y * max_displacement * h;
  displacement[2] = delta_z * max_displacement * h;

  /* Move the spart */
  sp->x[0] += displacement[0];
  sp->x[1] += displacement[1];
  sp->x[2] += displacement[2];

  /* Copy the position to the gpart */
  sp->gpart->x[0] = sp->x[0];
  sp->gpart->x[1] = sp->x[1];
  sp->gpart->x[2] = sp->x[2];

  /* Do the gas particle. */
  const double mass_ratio = sp->mass / hydro_get_mass(pind);
  const double dx[3] = {mass_ratio * displacement[0],
                        mass_ratio * displacement[1],
                        mass_ratio * displacement[2]};

 part_set_x_ind(pind, 0, x[0] - dx[0]);
 part_set_x_ind(pind, 1, x[1] - dx[1]);
 part_set_x_ind(pind, 2, x[2] - dx[2]);

  /* Compute offsets since last cell construction */
  part_set_x_diff_ind(pind, 0, part_get_x_diff_ind(pind, 0) + dx[0]);
  part_set_x_diff_ind(pind, 1, part_get_x_diff_ind(pind, 1) + dx[1]);
  part_set_x_diff_ind(pind, 2, part_get_x_diff_ind(pind, 2) + dx[2]);
  part_set_x_diff_sort_ind(pind, 0, part_get_x_diff_sort_ind(pind, 0)+ dx[0]);
  part_set_x_diff_sort_ind(pind, 1, part_get_x_diff_sort_ind(pind, 1)+ dx[1]);
  part_set_x_diff_sort_ind(pind, 2, part_get_x_diff_sort_ind(pind, 2)+ dx[2]);

  /* Copy the position to the gpart */
  struct gpart* restrict gp = part_get_gpart(pind);
  gp->x[0] =part_get_x_ind(pind, 0);
  gp->x[1] =part_get_x_ind(pind, 1);
  gp->x[2] =part_get_x_ind(pind, 2);
}

/**
 * @brief Copies the properties of the gas particle over to the
 * star particle.
 *
 * @param pind Index of the particle.
 * @param sp the new created star particle with its properties.
 * @param e The #engine
 * @param starform the star formation law properties to use.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param phys_const the physical constants in internal units.
 * @param hydro_props The #hydro_props.
 * @param us The #unit_system.
 * @param cooling The #cooling_function_data.
 * @param chem_data The global properties of the chemistry scheme.
 * @param convert_part Did we convert a part (or spawned one)?
 * @param (return) displacement The 3D displacement vector of the star with
 * respect to the sink position.
 */
INLINE static void star_formation_copy_properties(
    size_t pind, struct spart *sp, const struct engine *e,
    const struct star_formation *starform, const struct cosmology *cosmo,
    const int with_cosmology, const struct phys_const *phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cooling_function_data *restrict cooling,
    const struct chemistry_global_data *chem_data, const int convert_part,
    float displacement[3]) {

  /* Initialize the feedback */
  feedback_init_after_star_formation(sp, e->feedback_props, star_population);
  struct gpart* restrict gp = part_get_gpart(pind);

  /* Store the current mass */
  const float mass_gas = hydro_get_mass(pind);
  if (!convert_part) {
    const float mass_star = starform->mass_stars;
    const float new_mass_gas = mass_gas - mass_star;

    /* Update the spart */
    sp->mass = mass_star;
    sp->gpart->mass = mass_star;

    /* Update the part */
    hydro_set_mass(pind, new_mass_gas);
    gp->mass = new_mass_gas;

    star_formation_separate_particles(e, pind, sp, displacement);
  } else {
    sp->mass = mass_gas;
  }
  sp->sf_data.birth_mass = sp->mass;

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the tracers data */
  // sp->tracers_data = xp->tracers_data;

  /* Move over the splitting data */
  sp->split_data = part_get_split_data(pind);

  /* Store the birth density and temperature in the star particle */
  const float birth_density = hydro_get_physical_density(pind, cosmo);
  const float birth_temperature = cooling_get_temperature(
      phys_const, hydro_props, us, cosmo, cooling, pind);

  star_formation_set_spart_birth_density(sp, birth_density);
  star_formation_set_spart_birth_temperature(sp, birth_temperature);

  /* Copy the progenitor id */
  star_formation_set_spart_progenitor_id(sp,part_get_id(pind));

  /* Copy the chemistry properties */
  chemistry_copy_star_formation_properties(pind, sp, chem_data, cosmo);
}

/**
 * @brief Prints the used parameters of the star formation law
 *
 * @param starform the star formation law properties.
 */
INLINE static void starformation_print_backend(
    const struct star_formation *starform) {
  message("Star formation law is 'GEAR'");
}

/**
 * @brief Return the star formation rate of a particle.
 *
 * This scheme does not store the SFR in the particles. We return 0.
 *
 * @param pind Index of the particle.
 */
INLINE static float star_formation_get_SFR(const size_t pind) {
  return 0.f;
}

/**
 * @brief Finishes the density calculation.
 *
 * Nothing to do here.
 *
 * @param pind Index of the particle to act upon
 * @param cd The global star_formation information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void star_formation_end_density(
    size_t pind, const struct star_formation *cd, const struct cosmology *cosmo) {

  struct star_formation_xpart_data *restrict sd = part_get_sf_data_p(pind);

#ifdef SPHENIX_SPH
  /* Copy the velocity divergence */
  sd->div_v = part_get_div_v(pind);
  /* SPHENIX is already including the Hubble flow */
#elif GADGET2_SPH
  /* Copy the velocity divergence */
  sd->div_v = part_get_div_v(pind);
  sd->div_v += hydro_dimension * cosmo->H;
#elif MINIMAL_SPH
  /* Copy the velocity divergence */
  sd->div_v = part_get_div_v(pind);
  sd->div_v += hydro_dimension * cosmo->H;
#elif GASOLINE_SPH
  /* Copy the velocity divergence */
  sd->div_v = (1. / 3.) * (p->viscosity.velocity_gradient[0][0] +
                                   p->viscosity.velocity_gradient[1][1] +
                                   p->viscosity.velocity_gradient[2][2]);
#elif HOPKINS_PU_SPH
  sd->div_v = part_get_div_v(pind);
#elif defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
  float dummy[3], gradvx[3], gradvy[3], gradvz[3];
  hydro_part_get_gradients(pind, dummy, gradvx, gradvy, gradvz, dummy);
  float div_v = gradvx[0] + gradvy[1] + gradvz[2];

  /* Multiply by the missing scale factors */
  div_v *= cosmo->a2_inv;

  /* Add the missing term */
  div_v += hydro_dimension * cosmo->H;
  sd->div_v = div_v;
#else
#error \
    "This scheme is not implemented. Note that Different scheme apply the Hubble flow in different places. Be careful about it."
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #part has 0 ngbs.
 *
 * Nothing to do here.
 *
 * @param pind Index of the particle to act upon
 * @param cd #star_formation containing star_formation informations.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
star_formation_part_has_no_neighbours(size_t pind,
                                      const struct star_formation *cd,
                                      const struct cosmology *cosmo) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * state to start the density loop.
 *
 * Nothing to do here.
 *
 * @param data The global star_formation information used for this run.
 * @param pind Index of the particle.
 */
__attribute__((always_inline)) INLINE static void star_formation_init_part(
    size_t pind, const struct star_formation *data) {}

/**
 * @brief Sets the star_formation properties of the (x-)particles to a valid
 * start state at the beginning of the simulation after the ICs have been read.
 *
 * Nothing to do here.
 *
 * @param phys_const The physical constant in internal units.
 * @param us The unit system.
 * @param cosmo The current cosmological model.
 * @param data The global star_formation information used for this run.
 * @param pind Index of the particle.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_part(const struct phys_const *restrict phys_const,
                               const struct unit_system *restrict us,
                               const struct cosmology *restrict cosmo,
                               const struct star_formation *data,
                               const size_t pind) {}

/**
 * @brief Split the star formation content of a particle into n pieces
 *
 *
 * @param pind Index of the hydro particle.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void star_formation_split_part(
    size_t pind, const double n) {
  error("Loic: to be implemented");
}

/**
 * @brief Deal with the case where no spart are available for star formation.
 *
 * @param e The #engine.
 * @param pind Index of the hydro particle.
 */
__attribute__((always_inline)) INLINE static void
star_formation_no_spart_available(const struct engine *e, const size_t pind) {
  error(
      "Failed to get a new particle. Please increase "
      "Scheduler:cell_extra_sparts "
      "or Scheduler:cell_extra_gparts");
}

/**
 * @brief Compute some information for the star formation model based
 * on all the particles that were read in.
 *
 * This is called once on start-up of the code.
 *
 * @param star_form The #star_formation structure.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static void
star_formation_first_init_stats(struct star_formation *star_form,
                                const struct engine *e) {

  const struct space *s = e->s;
  double avg_mass = 0;

  /* Sum the mass over all the particles. */
  for (size_t i = 0; i < s->nr_parts; i++) {
    avg_mass += hydro_get_mass(i);
  }

#ifdef WITH_MPI
  /* Compute the contribution from the other nodes. */
  MPI_Allreduce(MPI_IN_PLACE, &avg_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  star_form->mass_stars =
      avg_mass / (e->total_nr_parts * star_form->n_stars_per_part);

  if (e->nodeID == 0) {
    message("Mass new stars: %g", star_form->mass_stars);
  }
}

#endif /* SWIFT_GEAR_STAR_FORMATION_H */
