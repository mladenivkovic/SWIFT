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
#include <hdf5.h>
#include <math.h>
#include <stdlib.h>

/* Includes. */
#include "swift.h"
#include "write-hdf5-output.h"

/* Generates a SWIFT IC file, replicating makeInput.py
 * L is the number of particles along one axis
 * filename is the path to write to */
static void generate_input_hdf5(size_t L, const char *filename) {

  const double boxSize = 1.0;
  const int periodic = 1;
  const double density = 2.0;
  const double P = 1.0;
  const double gamma_val = 5.0 / 3.0;
  const int material = 0;

  const size_t numPart = L * L * L;
  const double mass = boxSize * boxSize * boxSize * density / (double)numPart;
  const double internalEnergy = P / ((gamma_val - 1.0) * density);
  const double he_density = density * 0.24;

  /* Allocate flat buffers matching each dataset's shape */
  double *coords = malloc(numPart * 3 * sizeof(double));
  float *v = malloc(numPart * 3 * sizeof(float));
  float *m = malloc(numPart * sizeof(float));
  float *h = malloc(numPart * sizeof(float));
  float *u = malloc(numPart * sizeof(float));
  float *rho = malloc(numPart * sizeof(float));
  unsigned long *ids = malloc(numPart * sizeof(unsigned long));
  int *mat = malloc(numPart * sizeof(int));
  float *he = malloc(numPart * sizeof(float));

  for (size_t i = 0; i < L; ++i) {
    for (size_t j = 0; j < L; ++j) {
      for (size_t k = 0; k < L; ++k) {
        size_t index = i * L * L + j * L + k;

        coords[index * 3 + 0] = i * boxSize / L + boxSize / (2 * L);
        coords[index * 3 + 1] = j * boxSize / L + boxSize / (2 * L);
        coords[index * 3 + 2] = k * boxSize / L + boxSize / (2 * L);

        v[index * 3 + 0] = 0.0f;
        v[index * 3 + 1] = 0.0f;
        v[index * 3 + 2] = 0.0f;

        m[index] = (float)mass;
        h[index] = (float)(2.251 * boxSize / L);
        u[index] = (float)internalEnergy;
        rho[index] = (float)density;
        ids[index] = index;
        mat[index] = material;
        he[index] = (float)he_density;
      }
    }
  }

  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* ---- /Header ---- */
  hid_t header =
      H5Gcreate2(file_id, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hid_t scalar = H5Screate(H5S_SCALAR);
  hid_t attr;

  attr = H5Acreate2(header, "BoxSize", H5T_NATIVE_DOUBLE, scalar, H5P_DEFAULT,
                    H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_DOUBLE, &boxSize);
  H5Aclose(attr);

  double time_val = 0.0;
  attr = H5Acreate2(header, "Time", H5T_NATIVE_DOUBLE, scalar, H5P_DEFAULT,
                    H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_DOUBLE, &time_val);
  H5Aclose(attr);

  int num_files = 1;
  attr = H5Acreate2(header, "NumFilesPerSnapshot", H5T_NATIVE_INT, scalar,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, &num_files);
  H5Aclose(attr);
  H5Sclose(scalar);

  hsize_t six = 6;
  hid_t vec6 = H5Screate_simple(1, &six, NULL);

  long long numpart_total[6] = {(long long)numPart, 0, 0, 0, 0, 0};
  long long numpart_highword[6] = {0, 0, 0, 0, 0, 0};
  long long numpart_thisfile[6] = {(long long)numPart, 0, 0, 0, 0, 0};
  double mass_table[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int flag_entropy[6] = {0, 0, 0, 0, 0, 0};

  attr = H5Acreate2(header, "NumPart_Total", H5T_NATIVE_LLONG, vec6,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_LLONG, numpart_total);
  H5Aclose(attr);

  attr = H5Acreate2(header, "NumPart_Total_HighWord", H5T_NATIVE_LLONG, vec6,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_LLONG, numpart_highword);
  H5Aclose(attr);

  attr = H5Acreate2(header, "NumPart_ThisFile", H5T_NATIVE_LLONG, vec6,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_LLONG, numpart_thisfile);
  H5Aclose(attr);

  attr = H5Acreate2(header, "MassTable", H5T_NATIVE_DOUBLE, vec6, H5P_DEFAULT,
                    H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_DOUBLE, mass_table);
  H5Aclose(attr);

  attr = H5Acreate2(header, "Flag_Entropy_ICs", H5T_NATIVE_INT, vec6,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, flag_entropy);
  H5Aclose(attr);

  H5Sclose(vec6);
  H5Gclose(header);

  /* ---- /RuntimePars ---- */
  hid_t runtime = H5Gcreate2(file_id, "/RuntimePars", H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
  scalar = H5Screate(H5S_SCALAR);
  attr = H5Acreate2(runtime, "PeriodicBoundariesOn", H5T_NATIVE_INT, scalar,
                    H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, &periodic);
  H5Aclose(attr);
  H5Sclose(scalar);
  H5Gclose(runtime);

  /* ---- /Units ---- */
  hid_t units =
      H5Gcreate2(file_id, "/Units", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  scalar = H5Screate(H5S_SCALAR);
  double one = 1.0;
  const char *unit_names[5] = {
      "Unit length in cgs (U_L)", "Unit mass in cgs (U_M)",
      "Unit time in cgs (U_t)", "Unit current in cgs (U_I)",
      "Unit temperature in cgs (U_T)"};
  for (int n = 0; n < 5; ++n) {
    attr = H5Acreate2(units, unit_names[n], H5T_NATIVE_DOUBLE, scalar,
                      H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &one);
    H5Aclose(attr);
  }
  H5Sclose(scalar);
  H5Gclose(units);

  /* ---- /PartType0 ---- */
  hid_t part0 =
      H5Gcreate2(file_id, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  hsize_t dims3[2] = {numPart, 3};
  hsize_t dims1[1] = {numPart};

  hid_t space3 = H5Screate_simple(2, dims3, NULL);
  hid_t space1 = H5Screate_simple(1, dims1, NULL);

  hid_t ds;

  ds = H5Dcreate2(part0, "Coordinates", H5T_NATIVE_DOUBLE, space3, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "Velocities", H5T_NATIVE_FLOAT, space3, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "Masses", H5T_NATIVE_FLOAT, space1, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "SmoothingLength", H5T_NATIVE_FLOAT, space1,
                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, h);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "InternalEnergy", H5T_NATIVE_FLOAT, space1,
                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "Density", H5T_NATIVE_FLOAT, space1, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "ParticleIDs", H5T_NATIVE_ULONG, space1, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ids);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "MaterialIDs", H5T_NATIVE_INT, space1, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat);
  H5Dclose(ds);

  ds = H5Dcreate2(part0, "HeDensity", H5T_NATIVE_FLOAT, space1, H5P_DEFAULT,
                  H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, he);
  H5Dclose(ds);

  H5Sclose(space3);
  H5Sclose(space1);
  H5Gclose(part0);
  H5Fclose(file_id);

  free(coords);
  free(v);
  free(m);
  free(h);
  free(u);
  free(rho);
  free(ids);
  free(mat);
  free(he);
}

static void select_output_engine_init(struct engine *e, struct space *s,
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

static void select_output_space_init(struct space *s, double *dim,
                                     int periodic, size_t Ngas,
                                     size_t Nspart, size_t Ngpart,
                                     struct part *parts, struct spart *sparts,
                                     struct gpart *gparts) {
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

static void select_output_space_clean(struct space *s) { free(s->xparts); };

static void select_output_engine_clean(struct engine *e) {
  threadpool_clean(&e->threadpool);
}

int write_hdf5_output_run(int numberOfParticles, const char *param_filename) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  message("Number of particles requested: %d", numberOfParticles);

 // get the approx cube root of number of particles for side length of cube
  size_t L = (size_t)ceil(cbrt((double)numberOfParticles));

  message("Generating IC with L=%zu (%zu particles).", L, L * L * L);
  generate_input_hdf5(L, "write-hdf5-output.hdf5");

  size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nspart = 0, Nbpart = 0,
         Nsink = 0, Nnupart = 0;
  int flag_entropy_ICs = -1;
  int periodic = 1;
  double dim[3];
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;
  struct sink *sinks = NULL;
  // struct ic_info ics_metadata;
  // strcpy(ics_metadata.group_name, "NoSUCH");

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

  /* Read data */
  message("Reading initial conditions.");
  read_ic_single("write-hdf5-output.hdf5", &us, dim, &parts, &gparts, &sinks,
                 &sparts, &bparts, &Ngas, &Ngpart, &Ngpart_background,
                 &Nnupart, &Nsink, &Nspart, &Nbpart, &flag_entropy_ICs,
                 /*with_hydro=*/1,
                 /*with_gravity=*/0,
                 /*with_sink=*/0,
                 /*with_stars=*/0,
                 /*with_black_holes=*/0,
                 /*with_cosmology=*/0,
                 /*cleanup_h=*/0,
                 /*cleanup_sqrt_a=*/0,
                 /*h=*/1., /*a=*/1., /*n_threads=*/1, /*dry_run=*/0,
                 /*remap_ids=*/0, &ics_metadata);

  /* pseudo initialization of the space */
  message("Initialization of the space.");
  struct space s;
  select_output_space_init(&s, dim, periodic, Ngas, Nspart, Ngpart, parts,
                           sparts, gparts);

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
  e.physical_constants = &prog_const;
  sprintf(e.snapshot_base_name, "write-hdf5-output");
  sprintf(e.run_name, "HDF5 writing test");
  select_output_engine_init(&e, &s, &cosmo, &param_file, &output_options,
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
  select_output_engine_clean(&e);
  select_output_space_clean(&s);
  cosmology_clean(&cosmo);
  free(parts);
  free(gparts);

  return 0;
}

int main(int argc, char *argv[]) {

  int numberOfParticles = 10; // default amount

  if (argc > 1) {
    FILE *file = fopen(argv[1], "r");
    if (file == NULL) {
      fprintf(stderr, "Error: Could not open file %s\n", argv[1]);
      return 1;
    }
    fscanf(file, "numberOfParticles = %d", &numberOfParticles);
    fclose(file);
  }

  return write_hdf5_output_run(numberOfParticles, "HDF5WritingParameters.yml");
}
