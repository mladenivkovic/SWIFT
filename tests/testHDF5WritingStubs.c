#include <stdio.h>
#include <stdint.h>

struct engine;
struct space;
struct neutrino_model;
struct gpart;
struct cosmology;
struct lightcone_props;
struct lightcone_map;

void engine_struct_dump(struct engine *e, FILE *stream) {
  (void)e;
  (void)stream;
}
void engine_struct_restore(struct engine *e, FILE *stream) {
  (void)e;
  (void)stream;
}
void gather_neutrino_consts(const struct space *s, struct neutrino_model *nm) {
  (void)s;
  (void)nm;
}

void gpart_neutrino_weight_mesh_only(const struct gpart *gp,
                                     const struct neutrino_model *nm,
                                     double *weight) {
  (void)gp;
  (void) nm;
  (void) weight;
}

double lightcone_map_neutrino_baseline_value(
    const struct cosmology *c, const struct lightcone_props *lightcone_props,
    const struct lightcone_map *map) {
(void)c;
(void)lightcone_props;
(void)map;
return 0;
}

double neutrino_seed_to_fermi_dirac(uint64_t seed) {
(void)seed;
return 0;
}

void gpart_neutrino_mass_weight(const struct gpart *gp,
                                const struct neutrino_model *nm, double *mass,
                                double *weight) {
(void)gp;
(void)nm;
(void)mass;
(void)weight;
}

double lightcone_map_neutrino_mass_get_value(
    const struct engine *e, const struct lightcone_props *lightcone_props,
    const struct gpart *gp, const double a_cross, const double x_cross[3]) {
(void)e;
(void)lightcone_props;
(void)gp;
(void)a_cross;
(void)x_cross;
return 0;
}

int lightcone_map_neutrino_mass_type_contributes(int ptype) { return 0; }

int engine_rank = 0;
const char *engine_policy_names[] = {"none",
                                     "rand",
                                     "steal",
                                     "keep",
                                     "block",
                                     "cpu tight",
                                     "mpi",
                                     "numa affinity",
                                     "hydro",
                                     "self gravity",
                                     "external gravity",
                                     "cosmological integration",
                                     "drift everything",
                                     "reconstruct multi-poles",
                                     "temperature",
                                     "cooling",
                                     "stars",
                                     "structure finding",
                                     "star formation",
                                     "feedback",
                                     "black holes",
                                     "fof search",
                                     "time-step limiter",
                                     "time-step sync",
                                     "csds",
                                     "line of sight",
                                     "sink",
                                     "rt",
                                     "power spectra",
                                     "moving mesh",
                                     "moving mesh hydro",
                                     "no_io"};
