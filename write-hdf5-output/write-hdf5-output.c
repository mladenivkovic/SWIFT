/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Elijah Cann (e.a.cann72@gmail.com).
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

/* Some standard headers. */
#include <config.h>
#include <stdlib.h>

/* Includes. */
#include "swift.h"

/* Box size for the test lattice. Used by both generate_particles()
 * and write_hdf5_output_run(). */
static const double kBoxSize = 1.0;

/* Builds numPart gas particles directly in memory
 * Returns an array of numPart particles. Caller must free it. */
static struct part *generate_particles(size_t numPart) {

  struct part *parts = NULL;
  if (swift_memalign("parts", (void **)&parts, part_align,
                     numPart * sizeof(struct part)) != 0)
    error("Failed to allocate parts.");
  bzero(parts, numPart * sizeof(struct part));

  for (size_t i = 0; i < numPart; ++i) {
    struct part *p = &parts[i];

    p->id = (long long)i;

    p->x[0] = 0.5 * kBoxSize;
    p->x[1] = 0.5 * kBoxSize;
    p->x[2] = 0.5 * kBoxSize;

    p->v[0] = 0.0f;
    p->v[1] = 0.0f;
    p->v[2] = 0.0f;

    p->mass = 1.0f;
    p->h = 0.1f;
    p->u = 1.0f;
    p->rho = 1.0f;
  }

  return parts;
}

static void hdf5_output_engine_init(struct engine *e, struct space *s,
                                    struct cosmology *cosmo,
                                    struct swift_params *params,
                                    struct output_options *output,
                                    struct cooling_function_data *cooling,
                                    struct hydro_props *hydro_properties,
                                    struct ic_info *ics_metadata) {
  /* set structures */
  e->s = s;
  e->cooling_func = cooling;
  e->parameter_file = params;
  e->output_options = output;
  e->cosmology = cosmo;
  e->policy = engine_policy_hydro;
  e->hydro_properties = hydro_properties;
  e->ics_metadata = ics_metadata;

  /* initialization of threadpool */
  threadpool_init(&e->threadpool, 1);

  /* set parameters */
  e->verbose = 1;
  e->time = 0;
  e->snapshot_output_count = 0;
  e->snapshot_compression = 0;
};

static void hdf5_output_space_init(struct space *s, double *dim, int periodic,
                                   size_t Ngas, size_t Nspart, size_t Ngpart,
                                   struct part *parts, struct spart *sparts,
                                   struct gpart *gparts) {
  // zero everything first
  bzero(s, sizeof(struct space));

  s->periodic = periodic;
  for (int i = 0; i < 3; i++) {
    s->dim[i] = dim[i];
  }

  /* init space particles */
  s->nr_parts = Ngas;
  s->nr_sparts = Nspart;
  s->nr_gparts = Ngpart;

  s->parts = parts;
  s->gparts = gparts;
  s->sparts = sparts;

  /* Allocate the extra parts array for the gas particles. */
  if (posix_memalign((void **)&s->xparts, xpart_align,
                     Ngas * sizeof(struct xpart)) != 0)
    error("Failed to allocate xparts.");
  bzero(s->xparts, Ngas * sizeof(struct xpart));
};

static void hdf5_output_space_clean(struct space *s) { free(s->xparts); };

static void hdf5_output_engine_clean(struct engine *e) {
  threadpool_clean(&e->threadpool);
}

/**
 * @brief Builds gas particles in memory and writes them out with
 *        write_output_single()
 *
 * @param numberOfParticles How many particles to make.
 * @param param_filename Path to the SWIFT parameter file
 *
 * @return 0 on success, non-zero on failure.
 */
static int write_hdf5_output_run(int numberOfParticles,
                                 const char *param_filename) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  message("Number of particles requested: %d", numberOfParticles);

  size_t Ngas = (size_t)numberOfParticles;

  message("Generating %zu particles directly in memory.", Ngas);
  struct part *parts = generate_particles(Ngas);

  int periodic = 1;
  double dim[3] = {kBoxSize, kBoxSize, kBoxSize};

  /* parse parameters */
  message("Reading parameters.");
  struct swift_params param_file;
  parser_read_file(param_filename, &param_file);

  struct ic_info ics_metadata;
  ic_info_init(&ics_metadata, &param_file);

  struct output_options output_options;
  output_options_init(&param_file, 0, &output_options);

  /* Default unit system */
  message("Initialization of the unit system.");
  struct unit_system us;
  units_init_cgs(&us);

  /* Default physical constants */
  message("Initialization of the physical constants.");
  struct phys_const prog_const;
  phys_const_init(&us, &param_file, &prog_const);

  /* pseudo initialization of the space */
  message("Initialization of the space.");
  struct space s;
  hdf5_output_space_init(&s, dim, periodic, Ngas, /*Nspart=*/0, /*Ngpart=*/0,
                         parts, /*sparts=*/NULL, /*gparts=*/NULL);

  /* initialization of cosmology */
  message("Initialization of the cosmology.");
  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  /* pseudo initialization of cooling */
  message("Initialization of the cooling.");
  struct cooling_function_data cooling;

  /* pseudo initialization of hydro */
  message("Initialization of the hydro.");
  struct hydro_props hydro_properties;
  hydro_props_init(&hydro_properties, &prog_const, &us, &param_file);

  /* pseudo initialization of the engine */
  message("Initialization of the engine.");
  struct engine e;

  bzero(&e, sizeof(struct engine));
  e.physical_constants = &prog_const;
  sprintf(e.snapshot_base_name, "write-hdf5-output");
  sprintf(e.run_name, "HDF5 writing test");
  hdf5_output_engine_init(&e, &s, &cosmo, &param_file, &output_options,
                          &cooling, &hydro_properties, &ics_metadata);

  /* check output selection */
  message("Checking output parameters.");
  io_prepare_output_fields(&output_options, /*with_cosmology=*/0,
                           /*with_fof=*/0,
                           /*with_structure_finding=*/0,
                           /*verbose=*/1);

  /* write output file */
  message("Writing output.");
  write_output_single(&e, &us, &us, /*fof=*/0);

  /* Clean-up */
  message("Cleaning memory.");
  hdf5_output_engine_clean(&e);
  hdf5_output_space_clean(&s);
  cosmology_clean(&cosmo);
  swift_free("parts", parts);

  return 0;
}

int main(int argc, char *argv[]) {
  int numberOfParticles = 16;                               // default amount
  const char *param_filename = "HDF5WritingParameters.yml"; // default

  if (argc > 1) {
    numberOfParticles = atoi(argv[1]);
    if (numberOfParticles <= 0) {
      fprintf(stderr, "Error: invalid numberOfParticles '%s'\n", argv[1]);
      return 1;
    }
  }

  if (argc > 2) {
    param_filename = argv[2];
  }

  return write_hdf5_output_run(numberOfParticles, param_filename);
}
