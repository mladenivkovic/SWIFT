#ifndef SWIFT_HYDRO_PART_SPHENIX_H
#define SWIFT_HYDRO_PART_SPHENIX_H

 /*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) & Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#include "adaptive_softening_struct.h"
#include "align.h"
#include "black_holes_struct.h"
#include "chemistry_struct.h"
#include "cooling_struct.h"
#include "csds.h"
#include "error.h"
#include "feedback_struct.h"
#include "fvpm_geometry_struct.h"
#include "hydro_part_arrays_struct.h"
#include "mhd_struct.h"
#include "particle_splitting_struct.h"
#include "pressure_floor_struct.h"
#include "rt_struct.h"
#include "sink_struct.h"
#include "star_formation_struct.h"
#include "timestep_limiter_struct.h"
#include "tracers_struct.h"

#include <float.h>
#include <limits.h>
#include <stddef.h>


extern struct part_arrays global_part_arrays;



struct part {
  /*! the particle position */
  double _x[3];

  /*! Particle smoothing length */
  float _h;

  /*! Particle predicted velocity */
  float _v[3];

  /*! Particle mass */
  float _mass;

  /*! Particle acceleration */
  float _a_hydro[3];

  /*! Time derivative of the internal energy */
  float _u_dt;

  /*! Particle internal energy */
  float _u;

  /*! Particle density */
  float _rho;

  /*! Store density/force specific stuff */
  union {

    /*! Structure for the variables only used in the density loop over neighbours Quantities in this sub-structure should only be accessed in the density loop over neighbours and the ghost task */
    struct {

      /*! Derivative of density with respect to h */
      float _rho_dh;

      /*! Neighbour number count */
      float _wcount;

      /*! Derivative of the neighbour number with respect to h */
      float _wcount_dh;

      /*! Particle velocity curl */
      float _rot_v[3];

     }_density;

    /*! Structure for the variables only used in the force loop over neighbours Quantities in this sub-structure should only be accessed in the force loop over neighbours and the ghost, drift and kick tasks */
    struct {

      /*! 'Grad h' term -- only partial in P-U */
      float _f_gradh;

      /*! Particle pressure */
      float _pressure;

      /*! Balsara switch */
      float _balsara;

      /*! Particle soundspeed */
      float _soundspeed;

      /*! Maximal alpha (viscosity) over neighbours */
      float _alpha_visc_max_ngb;

      /*! Time derivative of smoothing length */
      float _h_dt;

     }_force;

   };

  /*! store viscosity information in a separate struct */
  struct {

    /*! Particle velocity divergence */
    float _div_v;

    /*! Time differential of velocity divergence */
    float _div_v_dt;

    /*! Particle velocity divergence from previous step */
    float _div_v_previous_step;

    /*! Artificial viscosity parameter */
    float _alpha_av;

    /*! Signal velocity */
    float _v_sig;

   }_viscosity;

  /*! Store thermal diffusion information in a separate struct */
  struct {

    /*! del^2 u, a smoothed quantity */
    float _laplace_u;

    /*! Thermal diffusion coefficient */
    float _alpha_diff;

   }_diffusion;

  /*! Particle unique ID */
  long long _id;

  /*! Pointer to corresponding gravity part */
  struct gpart* _gpart;

  /*! Additional data used for adaptive softening */
  struct adaptive_softening_part_data _adaptive_softening_data;

  /*! Additional data used by the MHD scheme */
  struct mhd_part_data _mhd_data;

  /*! Chemistry information */
  struct chemistry_part_data _chemistry_data;

  /*! Cooling information */
  struct cooling_part_data _cooling_data;

  /*! Additional data used by the feedback */
  struct feedback_part_data _feedback_data;

  /*! Black holes information (eg swallowing ID) */
  struct black_holes_part_data _black_holes_data;

  /*! Sink information (eg swallowing ID) */
  struct sink_part_data _sink_data;

  /*! Additional data used by the pressure floor */
  struct pressure_floor_part_data _pressure_floor_data;

  /*! Additional Radiative Transfer Data */
  struct rt_part_data _rt_data;

  /*! RT sub-cycling time stepping data */
  struct rt_timestepping_data _rt_time_data;

  /*! Tree-depth at which size / 2 <= h * gamma < size */
  char _depth_h;

  /*! Time-step length */
  timebin_t _time_bin;

  /*! Time-step limiter information */
  struct timestep_limiter_data _limiter_data;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Time of the last drift */
  integertime_t _ti_drift;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /*! Time of the last drift */
  integertime_t _ti_kick;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Integer number of neighbours in the density loop */
  int _N_density;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact integer number of neighbours in the density loop */
  int _N_density_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Integer number of neighbours in the gradient loop */
  int _N_gradient;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact integer number of neighbours in the gradient loop */
  int _N_gradient_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Integer number of neighbours in the force loop */
  int _N_force;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact integer number of neighbours in the force loop */
  int _N_force_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact value of the density field obtained via brute-force loop */
  float _rho_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Weighted number of neighbours in the density loop */
  float _n_density;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact value of the weighted number of neighbours in the density loop */
  float _n_density_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Weighted number of neighbours in the gradient loop */
  float _n_gradient;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact value of the weighted number of neighbours in the gradient loop */
  float _n_gradient_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Weighted number of neighbours in the force loop */
  float _n_force;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Exact value of the weighted number of neighbours in the force loop */
  float _n_force_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Has this particle interacted with any unhibited neighbour? */
  char _inhibited_exact;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /*! Has this particle been woken up by the limiter? */
  char _limited_part;
#endif

  /*! Geometrical quantities used for Finite Volume Particle Method RT */
  struct fvpm_geometry_struct _geometry;

#ifdef SWIFT_DEBUG_CHECKS
  /*! This particle's accessor ID, identical for all structs associated with this particle. */
  long long _accessor_id;
#endif

} SWIFT_STRUCT_ALIGN_PART;

struct xpart {
  /*! Offset between current position and position at last tree rebuild. */
  float _x_diff[3];

  /*! Offset between the current position and position at the last sort. */
  float _x_diff_sort[3];

  /*! Velocity at the last full step. */
  float _v_full[3];

  /*! Gravitational acceleration at the end of the last step */
  float _a_grav[3];

  /*! Internal energy at the last full step. */
  float _u_full;

  /*! Additional data used to record particle splits */
  struct particle_splitting_data _split_data;

  /*! Additional data used to record cooling information */
  struct cooling_xpart_data _extra_cooling_data;

  /*! Additional data used by the tracers */
  struct tracers_xpart_data _extra_tracers_data;

  /*! Additional data used by the tracers */
  struct star_formation_xpart_data _extra_sf_data;

  /*! Additional data used by the feedback */
  struct feedback_xpart_data _extra_feedback_data;

  /*! Additional data used by the MHD scheme */
  struct mhd_xpart_data _extra_mhd_data;

#ifdef WITH_CSDS
  /*! Additional data for the particle csds */
  struct csds_part_data _csds_data;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /*! This particle's accessor ID, identical for all structs associated with this particle. */
  long long _accessor_id;
#endif

} SWIFT_STRUCT_ALIGN_XPART;


/**
 * @brief get x, the particle position,
 * for read and write access. For read-only access, use
 * part_get_const_x() instead.
 */
static __attribute__((always_inline)) INLINE double*
  part_get_x(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_x;
}

/**
 * @brief get x, the particle position, for read-only access.
 */
static __attribute__((always_inline)) INLINE const double*
  part_get_const_x(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_x;
}

/**
 * @brief get x, the particle position, by index.
 */
static __attribute__((always_inline)) INLINE double
  part_get_x_ind(const size_t pind, const int i) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_x[i];
}

/**
 * @brief set all values of x, the particle position,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x(const size_t pind, const double x[3]) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_x[0] = x[0];
  part_s->_x[1] = x[1];
  part_s->_x[2] = x[2];
}

/**
 * @brief set the value of x, the particle position, by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x_ind(const size_t pind, const int i, const double x) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_x[i] = x;
}




/**
 * @brief get h, Particle smoothing length.
 */
static __attribute__((always_inline)) INLINE float
  part_get_h(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_h;
}


/**
 * @brief get a pointer to h, Particle smoothing length.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to h. If you need read-only access to h, use part_get_const_h_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_h_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_h;
}


/**
 * @brief get read-only access to pointer to h,
 * Particle smoothing length.
 * If you need write access to h, use part_get_h_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_h_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_h;
}


/**
 * @brief set the value of h, Particle smoothing length.
 */
static __attribute__((always_inline)) INLINE void
  part_set_h(const size_t pind, const float h) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_h = h;
}





/**
 * @brief get v, Particle predicted velocity,
 * for read and write access. For read-only access, use
 * part_get_const_v() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_v(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_v;
}

/**
 * @brief get v, Particle predicted velocity, for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_v(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_v;
}

/**
 * @brief get v, Particle predicted velocity, by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_v_ind(const size_t pind, const int i) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_v[i];
}

/**
 * @brief set all values of v, Particle predicted velocity,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_v(const size_t pind, const float v[3]) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_v[0] = v[0];
  part_s->_v[1] = v[1];
  part_s->_v[2] = v[2];
}

/**
 * @brief set the value of v, Particle predicted velocity, by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_v_ind(const size_t pind, const int i, const float v) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_v[i] = v;
}




/**
 * @brief get mass, Particle mass.
 */
static __attribute__((always_inline)) INLINE float
  part_get_mass(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_mass;
}


/**
 * @brief get a pointer to mass, Particle mass.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to mass. If you need read-only access to mass, use part_get_const_mass_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_mass_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_mass;
}


/**
 * @brief get read-only access to pointer to mass,
 * Particle mass.
 * If you need write access to mass, use part_get_mass_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_mass_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_mass;
}


/**
 * @brief set the value of mass, Particle mass.
 */
static __attribute__((always_inline)) INLINE void
  part_set_mass(const size_t pind, const float mass) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_mass = mass;
}





/**
 * @brief get a_hydro, Particle acceleration,
 * for read and write access. For read-only access, use
 * part_get_const_a_hydro() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_a_hydro(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_a_hydro;
}

/**
 * @brief get a_hydro, Particle acceleration, for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_a_hydro(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_a_hydro;
}

/**
 * @brief get a_hydro, Particle acceleration, by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_a_hydro_ind(const size_t pind, const int i) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_a_hydro[i];
}

/**
 * @brief set all values of a_hydro, Particle acceleration,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_a_hydro(const size_t pind, const float a_hydro[3]) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_a_hydro[0] = a_hydro[0];
  part_s->_a_hydro[1] = a_hydro[1];
  part_s->_a_hydro[2] = a_hydro[2];
}

/**
 * @brief set the value of a_hydro, Particle acceleration, by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_a_hydro_ind(const size_t pind, const int i, const float a_hydro) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_a_hydro[i] = a_hydro;
}




/**
 * @brief get u_dt, Time derivative of the internal energy.
 */
static __attribute__((always_inline)) INLINE float
  part_get_u_dt(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_u_dt;
}


/**
 * @brief get a pointer to u_dt, Time derivative of the internal energy.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to u_dt. If you need read-only access to u_dt, use part_get_const_u_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_u_dt_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_u_dt;
}


/**
 * @brief get read-only access to pointer to u_dt,
 * Time derivative of the internal energy.
 * If you need write access to u_dt, use part_get_u_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_u_dt_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_u_dt;
}


/**
 * @brief set the value of u_dt, Time derivative of the internal energy.
 */
static __attribute__((always_inline)) INLINE void
  part_set_u_dt(const size_t pind, const float u_dt) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_u_dt = u_dt;
}





/**
 * @brief get u, Particle internal energy.
 */
static __attribute__((always_inline)) INLINE float
  part_get_u(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_u;
}


/**
 * @brief get a pointer to u, Particle internal energy.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to u. If you need read-only access to u, use part_get_const_u_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_u_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_u;
}


/**
 * @brief get read-only access to pointer to u,
 * Particle internal energy.
 * If you need write access to u, use part_get_u_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_u_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_u;
}


/**
 * @brief set the value of u, Particle internal energy.
 */
static __attribute__((always_inline)) INLINE void
  part_set_u(const size_t pind, const float u) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_u = u;
}





/**
 * @brief get rho, Particle density.
 */
static __attribute__((always_inline)) INLINE float
  part_get_rho(size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_rho;
}


/**
 * @brief get a pointer to rho, Particle density.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to rho. If you need read-only access to rho, use part_get_const_rho_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_rho_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rho;
}


/**
 * @brief get read-only access to pointer to rho,
 * Particle density.
 * If you need write access to rho, use part_get_rho_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_rho_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rho;
}


/**
 * @brief set the value of rho, Particle density.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rho(const size_t pind, const float rho) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_rho = rho;
}





/**
 * @brief get rho_dh, Derivative of density with respect to h.
 */
static __attribute__((always_inline)) INLINE float
  part_get_rho_dh(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._rho_dh;
}


/**
 * @brief get a pointer to rho_dh, Derivative of density with respect to h.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to rho_dh. If you need read-only access to rho_dh, use part_get_const_rho_dh_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_rho_dh_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._rho_dh;
}


/**
 * @brief get read-only access to pointer to rho_dh,
 * Derivative of density with respect to h.
 * If you need write access to rho_dh, use part_get_rho_dh_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_rho_dh_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._rho_dh;
}


/**
 * @brief set the value of rho_dh, Derivative of density with respect to h.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rho_dh(const size_t pind, const float rho_dh) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_density._rho_dh = rho_dh;
}



/**
 * @brief get wcount, Neighbour number count.
 */
static __attribute__((always_inline)) INLINE float
  part_get_wcount(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._wcount;
}


/**
 * @brief get a pointer to wcount, Neighbour number count.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to wcount. If you need read-only access to wcount, use part_get_const_wcount_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_wcount_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._wcount;
}


/**
 * @brief get read-only access to pointer to wcount,
 * Neighbour number count.
 * If you need write access to wcount, use part_get_wcount_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_wcount_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._wcount;
}


/**
 * @brief set the value of wcount, Neighbour number count.
 */
static __attribute__((always_inline)) INLINE void
  part_set_wcount(const size_t pind, const float wcount) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_density._wcount = wcount;
}



/**
 * @brief get wcount_dh, Derivative of the neighbour number with respect to h.
 */
static __attribute__((always_inline)) INLINE float
  part_get_wcount_dh(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._wcount_dh;
}


/**
 * @brief get a pointer to wcount_dh, Derivative of the neighbour number with respect to h.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to wcount_dh. If you need read-only access to wcount_dh, use part_get_const_wcount_dh_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_wcount_dh_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._wcount_dh;
}


/**
 * @brief get read-only access to pointer to wcount_dh,
 * Derivative of the neighbour number with respect to h.
 * If you need write access to wcount_dh, use part_get_wcount_dh_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_wcount_dh_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_density._wcount_dh;
}


/**
 * @brief set the value of wcount_dh, Derivative of the neighbour number with respect to h.
 */
static __attribute__((always_inline)) INLINE void
  part_set_wcount_dh(const size_t pind, const float wcount_dh) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_density._wcount_dh = wcount_dh;
}



/**
 * @brief get rot_v, Particle velocity curl,
 * for read and write access. For read-only access, use
 * part_get_const_rot_v() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_rot_v(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._rot_v;
}

/**
 * @brief get rot_v, Particle velocity curl, for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_rot_v(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._rot_v;
}

/**
 * @brief get rot_v, Particle velocity curl, by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_rot_v_ind(const size_t pind, const int i) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_density._rot_v[i];
}

/**
 * @brief set all values of rot_v, Particle velocity curl,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rot_v(const size_t pind, const float rot_v[3]) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_density._rot_v[0] = rot_v[0];
  part_s->_density._rot_v[1] = rot_v[1];
  part_s->_density._rot_v[2] = rot_v[2];
}

/**
 * @brief set the value of rot_v, Particle velocity curl, by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rot_v_ind(const size_t pind, const int i, const float rot_v) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_density._rot_v[i] = rot_v;
}


/**
 * @brief get f_gradh, 'Grad h' term -- only partial in P-U.
 */
static __attribute__((always_inline)) INLINE float
  part_get_f_gradh(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._f_gradh;
}


/**
 * @brief get a pointer to f_gradh, 'Grad h' term -- only partial in P-U.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to f_gradh. If you need read-only access to f_gradh, use part_get_const_f_gradh_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_f_gradh_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._f_gradh;
}


/**
 * @brief get read-only access to pointer to f_gradh,
 * 'Grad h' term -- only partial in P-U.
 * If you need write access to f_gradh, use part_get_f_gradh_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_f_gradh_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._f_gradh;
}


/**
 * @brief set the value of f_gradh, 'Grad h' term -- only partial in P-U.
 */
static __attribute__((always_inline)) INLINE void
  part_set_f_gradh(const size_t pind, const float f_gradh) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._f_gradh = f_gradh;
}



/**
 * @brief get pressure, Particle pressure.
 */
static __attribute__((always_inline)) INLINE float
  part_get_pressure(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._pressure;
}


/**
 * @brief get a pointer to pressure, Particle pressure.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to pressure. If you need read-only access to pressure, use part_get_const_pressure_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_pressure_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._pressure;
}


/**
 * @brief get read-only access to pointer to pressure,
 * Particle pressure.
 * If you need write access to pressure, use part_get_pressure_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_pressure_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._pressure;
}


/**
 * @brief set the value of pressure, Particle pressure.
 */
static __attribute__((always_inline)) INLINE void
  part_set_pressure(const size_t pind, const float pressure) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._pressure = pressure;
}



/**
 * @brief get balsara, Balsara switch.
 */
static __attribute__((always_inline)) INLINE float
  part_get_balsara(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._balsara;
}


/**
 * @brief get a pointer to balsara, Balsara switch.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to balsara. If you need read-only access to balsara, use part_get_const_balsara_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_balsara_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._balsara;
}


/**
 * @brief get read-only access to pointer to balsara,
 * Balsara switch.
 * If you need write access to balsara, use part_get_balsara_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_balsara_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._balsara;
}


/**
 * @brief set the value of balsara, Balsara switch.
 */
static __attribute__((always_inline)) INLINE void
  part_set_balsara(const size_t pind, const float balsara) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._balsara = balsara;
}



/**
 * @brief get soundspeed, Particle soundspeed.
 */
static __attribute__((always_inline)) INLINE float
  part_get_soundspeed(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._soundspeed;
}


/**
 * @brief get a pointer to soundspeed, Particle soundspeed.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to soundspeed. If you need read-only access to soundspeed, use part_get_const_soundspeed_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_soundspeed_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._soundspeed;
}


/**
 * @brief get read-only access to pointer to soundspeed,
 * Particle soundspeed.
 * If you need write access to soundspeed, use part_get_soundspeed_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_soundspeed_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._soundspeed;
}


/**
 * @brief set the value of soundspeed, Particle soundspeed.
 */
static __attribute__((always_inline)) INLINE void
  part_set_soundspeed(const size_t pind, const float soundspeed) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._soundspeed = soundspeed;
}



/**
 * @brief get alpha_visc_max_ngb, Maximal alpha (viscosity) over neighbours.
 */
static __attribute__((always_inline)) INLINE float
  part_get_alpha_visc_max_ngb(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._alpha_visc_max_ngb;
}


/**
 * @brief get a pointer to alpha_visc_max_ngb, Maximal alpha (viscosity) over neighbours.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to alpha_visc_max_ngb. If you need read-only access to alpha_visc_max_ngb, use part_get_const_alpha_visc_max_ngb_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_alpha_visc_max_ngb_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._alpha_visc_max_ngb;
}


/**
 * @brief get read-only access to pointer to alpha_visc_max_ngb,
 * Maximal alpha (viscosity) over neighbours.
 * If you need write access to alpha_visc_max_ngb, use part_get_alpha_visc_max_ngb_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_alpha_visc_max_ngb_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._alpha_visc_max_ngb;
}


/**
 * @brief set the value of alpha_visc_max_ngb, Maximal alpha (viscosity) over neighbours.
 */
static __attribute__((always_inline)) INLINE void
  part_set_alpha_visc_max_ngb(const size_t pind, const float alpha_visc_max_ngb) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._alpha_visc_max_ngb = alpha_visc_max_ngb;
}



/**
 * @brief get h_dt, Time derivative of smoothing length.
 */
static __attribute__((always_inline)) INLINE float
  part_get_h_dt(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_force._h_dt;
}


/**
 * @brief get a pointer to h_dt, Time derivative of smoothing length.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to h_dt. If you need read-only access to h_dt, use part_get_const_h_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_h_dt_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._h_dt;
}


/**
 * @brief get read-only access to pointer to h_dt,
 * Time derivative of smoothing length.
 * If you need write access to h_dt, use part_get_h_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_h_dt_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_force._h_dt;
}


/**
 * @brief set the value of h_dt, Time derivative of smoothing length.
 */
static __attribute__((always_inline)) INLINE void
  part_set_h_dt(const size_t pind, const float h_dt) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_force._h_dt = h_dt;
}





/**
 * @brief get div_v, Particle velocity divergence.
 */
static __attribute__((always_inline)) INLINE float
  part_get_div_v(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_viscosity._div_v;
}


/**
 * @brief get a pointer to div_v, Particle velocity divergence.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to div_v. If you need read-only access to div_v, use part_get_const_div_v_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_div_v_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v;
}


/**
 * @brief get read-only access to pointer to div_v,
 * Particle velocity divergence.
 * If you need write access to div_v, use part_get_div_v_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_div_v_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v;
}


/**
 * @brief set the value of div_v, Particle velocity divergence.
 */
static __attribute__((always_inline)) INLINE void
  part_set_div_v(const size_t pind, const float div_v) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_viscosity._div_v = div_v;
}



/**
 * @brief get div_v_dt, Time differential of velocity divergence.
 */
static __attribute__((always_inline)) INLINE float
  part_get_div_v_dt(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_viscosity._div_v_dt;
}


/**
 * @brief get a pointer to div_v_dt, Time differential of velocity divergence.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to div_v_dt. If you need read-only access to div_v_dt, use part_get_const_div_v_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_div_v_dt_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v_dt;
}


/**
 * @brief get read-only access to pointer to div_v_dt,
 * Time differential of velocity divergence.
 * If you need write access to div_v_dt, use part_get_div_v_dt_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_div_v_dt_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v_dt;
}


/**
 * @brief set the value of div_v_dt, Time differential of velocity divergence.
 */
static __attribute__((always_inline)) INLINE void
  part_set_div_v_dt(const size_t pind, const float div_v_dt) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_viscosity._div_v_dt = div_v_dt;
}



/**
 * @brief get div_v_previous_step, Particle velocity divergence from previous step.
 */
static __attribute__((always_inline)) INLINE float
  part_get_div_v_previous_step(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_viscosity._div_v_previous_step;
}


/**
 * @brief get a pointer to div_v_previous_step, Particle velocity divergence from previous step.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to div_v_previous_step. If you need read-only access to div_v_previous_step, use part_get_const_div_v_previous_step_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_div_v_previous_step_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v_previous_step;
}


/**
 * @brief get read-only access to pointer to div_v_previous_step,
 * Particle velocity divergence from previous step.
 * If you need write access to div_v_previous_step, use part_get_div_v_previous_step_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_div_v_previous_step_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._div_v_previous_step;
}


/**
 * @brief set the value of div_v_previous_step, Particle velocity divergence from previous step.
 */
static __attribute__((always_inline)) INLINE void
  part_set_div_v_previous_step(const size_t pind, const float div_v_previous_step) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_viscosity._div_v_previous_step = div_v_previous_step;
}



/**
 * @brief get alpha_av, Artificial viscosity parameter.
 */
static __attribute__((always_inline)) INLINE float
  part_get_alpha_av(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_viscosity._alpha_av;
}


/**
 * @brief get a pointer to alpha_av, Artificial viscosity parameter.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to alpha_av. If you need read-only access to alpha_av, use part_get_const_alpha_av_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_alpha_av_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._alpha_av;
}


/**
 * @brief get read-only access to pointer to alpha_av,
 * Artificial viscosity parameter.
 * If you need write access to alpha_av, use part_get_alpha_av_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_alpha_av_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._alpha_av;
}


/**
 * @brief set the value of alpha_av, Artificial viscosity parameter.
 */
static __attribute__((always_inline)) INLINE void
  part_set_alpha_av(const size_t pind, const float alpha_av) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_viscosity._alpha_av = alpha_av;
}



/**
 * @brief get v_sig, Signal velocity.
 */
static __attribute__((always_inline)) INLINE float
  part_get_v_sig(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_viscosity._v_sig;
}


/**
 * @brief get a pointer to v_sig, Signal velocity.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to v_sig. If you need read-only access to v_sig, use part_get_const_v_sig_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_v_sig_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._v_sig;
}


/**
 * @brief get read-only access to pointer to v_sig,
 * Signal velocity.
 * If you need write access to v_sig, use part_get_v_sig_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_v_sig_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_viscosity._v_sig;
}


/**
 * @brief set the value of v_sig, Signal velocity.
 */
static __attribute__((always_inline)) INLINE void
  part_set_v_sig(const size_t pind, const float v_sig) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_viscosity._v_sig = v_sig;
}





/**
 * @brief get laplace_u, del^2 u, a smoothed quantity.
 */
static __attribute__((always_inline)) INLINE float
  part_get_laplace_u(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_diffusion._laplace_u;
}


/**
 * @brief get a pointer to laplace_u, del^2 u, a smoothed quantity.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to laplace_u. If you need read-only access to laplace_u, use part_get_const_laplace_u_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_laplace_u_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_diffusion._laplace_u;
}


/**
 * @brief get read-only access to pointer to laplace_u,
 * del^2 u, a smoothed quantity.
 * If you need write access to laplace_u, use part_get_laplace_u_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_laplace_u_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_diffusion._laplace_u;
}


/**
 * @brief set the value of laplace_u, del^2 u, a smoothed quantity.
 */
static __attribute__((always_inline)) INLINE void
  part_set_laplace_u(const size_t pind, const float laplace_u) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_diffusion._laplace_u = laplace_u;
}



/**
 * @brief get alpha_diff, Thermal diffusion coefficient.
 */
static __attribute__((always_inline)) INLINE float
  part_get_alpha_diff(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_diffusion._alpha_diff;
}


/**
 * @brief get a pointer to alpha_diff, Thermal diffusion coefficient.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to alpha_diff. If you need read-only access to alpha_diff, use part_get_const_alpha_diff_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_alpha_diff_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_diffusion._alpha_diff;
}


/**
 * @brief get read-only access to pointer to alpha_diff,
 * Thermal diffusion coefficient.
 * If you need write access to alpha_diff, use part_get_alpha_diff_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_alpha_diff_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_diffusion._alpha_diff;
}


/**
 * @brief set the value of alpha_diff, Thermal diffusion coefficient.
 */
static __attribute__((always_inline)) INLINE void
  part_set_alpha_diff(const size_t pind, const float alpha_diff) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_diffusion._alpha_diff = alpha_diff;
}





/**
 * @brief get id, Particle unique ID.
 */
static __attribute__((always_inline)) INLINE long long
  part_get_id(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_id;
}


/**
 * @brief get a pointer to id, Particle unique ID.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to id. If you need read-only access to id, use part_get_const_id_p() instead.
 */
static __attribute__((always_inline)) INLINE long long*
  part_get_id_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_id;
}


/**
 * @brief get read-only access to pointer to id,
 * Particle unique ID.
 * If you need write access to id, use part_get_id_p() instead.
 */
static __attribute__((always_inline)) INLINE const long long*
  part_get_const_id_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_id;
}


/**
 * @brief set the value of id, Particle unique ID.
 */
static __attribute__((always_inline)) INLINE void
  part_set_id(const size_t pind, const long long id) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_id = id;
}





/**
 * @brief get gpart, Pointer to corresponding gravity part.
 */
static __attribute__((always_inline)) INLINE struct gpart*
  part_get_gpart(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_gpart;
}


/**
 * @brief get a pointer to gpart, Pointer to corresponding gravity part.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to gpart. If you need read-only access to gpart, use part_get_const_gpart_p() instead.
 */
static __attribute__((always_inline)) INLINE struct gpart**
  part_get_gpart_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_gpart;
}
/**
 * @brief set the value of gpart, Pointer to corresponding gravity part.
 */
static __attribute__((always_inline)) INLINE void
  part_set_gpart(const size_t pind,  struct gpart* gpart) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_gpart = gpart;
}





/**
 * @brief get adaptive_softening_data, Additional data used for adaptive softening.
 */
static __attribute__((always_inline)) INLINE struct adaptive_softening_part_data
  part_get_adaptive_softening_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_adaptive_softening_data;
}


/**
 * @brief get a pointer to adaptive_softening_data, Additional data used for adaptive softening.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to adaptive_softening_data. If you need read-only access to adaptive_softening_data, use part_get_const_adaptive_softening_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct adaptive_softening_part_data*
  part_get_adaptive_softening_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_adaptive_softening_data;
}


/**
 * @brief get read-only access to pointer to adaptive_softening_data,
 * Additional data used for adaptive softening.
 * If you need write access to adaptive_softening_data, use part_get_adaptive_softening_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct adaptive_softening_part_data*
  part_get_const_adaptive_softening_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_adaptive_softening_data;
}


/**
 * @brief set the value of adaptive_softening_data, Additional data used for adaptive softening.
 */
static __attribute__((always_inline)) INLINE void
  part_set_adaptive_softening_data(const size_t pind, const struct adaptive_softening_part_data adaptive_softening_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_adaptive_softening_data = adaptive_softening_data;
}





/**
 * @brief get mhd_data, Additional data used by the MHD scheme.
 */
static __attribute__((always_inline)) INLINE struct mhd_part_data
  part_get_mhd_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_mhd_data;
}


/**
 * @brief get a pointer to mhd_data, Additional data used by the MHD scheme.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to mhd_data. If you need read-only access to mhd_data, use part_get_const_mhd_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct mhd_part_data*
  part_get_mhd_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_mhd_data;
}


/**
 * @brief get read-only access to pointer to mhd_data,
 * Additional data used by the MHD scheme.
 * If you need write access to mhd_data, use part_get_mhd_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct mhd_part_data*
  part_get_const_mhd_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_mhd_data;
}


/**
 * @brief set the value of mhd_data, Additional data used by the MHD scheme.
 */
static __attribute__((always_inline)) INLINE void
  part_set_mhd_data(const size_t pind, const struct mhd_part_data mhd_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_mhd_data = mhd_data;
}





/**
 * @brief get chemistry_data, Chemistry information.
 */
static __attribute__((always_inline)) INLINE struct chemistry_part_data
  part_get_chemistry_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_chemistry_data;
}


/**
 * @brief get a pointer to chemistry_data, Chemistry information.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to chemistry_data. If you need read-only access to chemistry_data, use part_get_const_chemistry_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct chemistry_part_data*
  part_get_chemistry_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_chemistry_data;
}


/**
 * @brief get read-only access to pointer to chemistry_data,
 * Chemistry information.
 * If you need write access to chemistry_data, use part_get_chemistry_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct chemistry_part_data*
  part_get_const_chemistry_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_chemistry_data;
}


/**
 * @brief set the value of chemistry_data, Chemistry information.
 */
static __attribute__((always_inline)) INLINE void
  part_set_chemistry_data(const size_t pind, const struct chemistry_part_data chemistry_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_chemistry_data = chemistry_data;
}





/**
 * @brief get cooling_data, Cooling information.
 */
static __attribute__((always_inline)) INLINE struct cooling_part_data
  part_get_cooling_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_cooling_data;
}


/**
 * @brief get a pointer to cooling_data, Cooling information.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to cooling_data. If you need read-only access to cooling_data, use part_get_const_cooling_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct cooling_part_data*
  part_get_cooling_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_cooling_data;
}


/**
 * @brief get read-only access to pointer to cooling_data,
 * Cooling information.
 * If you need write access to cooling_data, use part_get_cooling_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct cooling_part_data*
  part_get_const_cooling_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_cooling_data;
}


/**
 * @brief set the value of cooling_data, Cooling information.
 */
static __attribute__((always_inline)) INLINE void
  part_set_cooling_data(const size_t pind, const struct cooling_part_data cooling_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_cooling_data = cooling_data;
}





/**
 * @brief get feedback_data, Additional data used by the feedback.
 */
static __attribute__((always_inline)) INLINE struct feedback_part_data
  part_get_feedback_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_feedback_data;
}


/**
 * @brief get a pointer to feedback_data, Additional data used by the feedback.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to feedback_data. If you need read-only access to feedback_data, use part_get_const_feedback_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct feedback_part_data*
  part_get_feedback_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_feedback_data;
}


/**
 * @brief get read-only access to pointer to feedback_data,
 * Additional data used by the feedback.
 * If you need write access to feedback_data, use part_get_feedback_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct feedback_part_data*
  part_get_const_feedback_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_feedback_data;
}


/**
 * @brief set the value of feedback_data, Additional data used by the feedback.
 */
static __attribute__((always_inline)) INLINE void
  part_set_feedback_data(const size_t pind, const struct feedback_part_data feedback_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_feedback_data = feedback_data;
}





/**
 * @brief get black_holes_data, Black holes information (eg swallowing ID).
 */
static __attribute__((always_inline)) INLINE struct black_holes_part_data
  part_get_black_holes_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_black_holes_data;
}


/**
 * @brief get a pointer to black_holes_data, Black holes information (eg swallowing ID).
 * Use this only if you need to modify the value, i.e. if you need write access
 * to black_holes_data. If you need read-only access to black_holes_data, use part_get_const_black_holes_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct black_holes_part_data*
  part_get_black_holes_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_black_holes_data;
}


/**
 * @brief get read-only access to pointer to black_holes_data,
 * Black holes information (eg swallowing ID).
 * If you need write access to black_holes_data, use part_get_black_holes_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct black_holes_part_data*
  part_get_const_black_holes_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_black_holes_data;
}


/**
 * @brief set the value of black_holes_data, Black holes information (eg swallowing ID).
 */
static __attribute__((always_inline)) INLINE void
  part_set_black_holes_data(const size_t pind, const struct black_holes_part_data black_holes_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_black_holes_data = black_holes_data;
}





/**
 * @brief get sink_data, Sink information (eg swallowing ID).
 */
static __attribute__((always_inline)) INLINE struct sink_part_data
  part_get_sink_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_sink_data;
}


/**
 * @brief get a pointer to sink_data, Sink information (eg swallowing ID).
 * Use this only if you need to modify the value, i.e. if you need write access
 * to sink_data. If you need read-only access to sink_data, use part_get_const_sink_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct sink_part_data*
  part_get_sink_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_sink_data;
}


/**
 * @brief get read-only access to pointer to sink_data,
 * Sink information (eg swallowing ID).
 * If you need write access to sink_data, use part_get_sink_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct sink_part_data*
  part_get_const_sink_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_sink_data;
}


/**
 * @brief set the value of sink_data, Sink information (eg swallowing ID).
 */
static __attribute__((always_inline)) INLINE void
  part_set_sink_data(const size_t pind, const struct sink_part_data sink_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_sink_data = sink_data;
}





/**
 * @brief get pressure_floor_data, Additional data used by the pressure floor.
 */
static __attribute__((always_inline)) INLINE struct pressure_floor_part_data
  part_get_pressure_floor_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_pressure_floor_data;
}


/**
 * @brief get a pointer to pressure_floor_data, Additional data used by the pressure floor.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to pressure_floor_data. If you need read-only access to pressure_floor_data, use part_get_const_pressure_floor_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct pressure_floor_part_data*
  part_get_pressure_floor_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_pressure_floor_data;
}


/**
 * @brief get read-only access to pointer to pressure_floor_data,
 * Additional data used by the pressure floor.
 * If you need write access to pressure_floor_data, use part_get_pressure_floor_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct pressure_floor_part_data*
  part_get_const_pressure_floor_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_pressure_floor_data;
}


/**
 * @brief set the value of pressure_floor_data, Additional data used by the pressure floor.
 */
static __attribute__((always_inline)) INLINE void
  part_set_pressure_floor_data(const size_t pind, const struct pressure_floor_part_data pressure_floor_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_pressure_floor_data = pressure_floor_data;
}





/**
 * @brief get rt_data, Additional Radiative Transfer Data.
 */
static __attribute__((always_inline)) INLINE struct rt_part_data
  part_get_rt_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_rt_data;
}


/**
 * @brief get a pointer to rt_data, Additional Radiative Transfer Data.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to rt_data. If you need read-only access to rt_data, use part_get_const_rt_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct rt_part_data*
  part_get_rt_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rt_data;
}


/**
 * @brief get read-only access to pointer to rt_data,
 * Additional Radiative Transfer Data.
 * If you need write access to rt_data, use part_get_rt_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct rt_part_data*
  part_get_const_rt_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rt_data;
}


/**
 * @brief set the value of rt_data, Additional Radiative Transfer Data.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rt_data(const size_t pind, const struct rt_part_data rt_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_rt_data = rt_data;
}





/**
 * @brief get rt_time_data, RT sub-cycling time stepping data.
 */
static __attribute__((always_inline)) INLINE struct rt_timestepping_data
  part_get_rt_time_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_rt_time_data;
}


/**
 * @brief get a pointer to rt_time_data, RT sub-cycling time stepping data.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to rt_time_data. If you need read-only access to rt_time_data, use part_get_const_rt_time_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct rt_timestepping_data*
  part_get_rt_time_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rt_time_data;
}


/**
 * @brief get read-only access to pointer to rt_time_data,
 * RT sub-cycling time stepping data.
 * If you need write access to rt_time_data, use part_get_rt_time_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct rt_timestepping_data*
  part_get_const_rt_time_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rt_time_data;
}


/**
 * @brief set the value of rt_time_data, RT sub-cycling time stepping data.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rt_time_data(const size_t pind, const struct rt_timestepping_data rt_time_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_rt_time_data = rt_time_data;
}





/**
 * @brief get depth_h, Tree-depth at which size / 2 <= h * gamma < size.
 */
static __attribute__((always_inline)) INLINE char
  part_get_depth_h(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_depth_h;
}


/**
 * @brief get a pointer to depth_h, Tree-depth at which size / 2 <= h * gamma < size.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to depth_h. If you need read-only access to depth_h, use part_get_const_depth_h_p() instead.
 */
static __attribute__((always_inline)) INLINE char*
  part_get_depth_h_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_depth_h;
}


/**
 * @brief get read-only access to pointer to depth_h,
 * Tree-depth at which size / 2 <= h * gamma < size.
 * If you need write access to depth_h, use part_get_depth_h_p() instead.
 */
static __attribute__((always_inline)) INLINE const char*
  part_get_const_depth_h_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_depth_h;
}


/**
 * @brief set the value of depth_h, Tree-depth at which size / 2 <= h * gamma < size.
 */
static __attribute__((always_inline)) INLINE void
  part_set_depth_h(const size_t pind, const char depth_h) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_depth_h = depth_h;
}





/**
 * @brief get time_bin, Time-step length.
 */
static __attribute__((always_inline)) INLINE timebin_t
  part_get_time_bin(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_time_bin;
}


/**
 * @brief get a pointer to time_bin, Time-step length.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to time_bin. If you need read-only access to time_bin, use part_get_const_time_bin_p() instead.
 */
static __attribute__((always_inline)) INLINE timebin_t*
  part_get_time_bin_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_time_bin;
}


/**
 * @brief get read-only access to pointer to time_bin,
 * Time-step length.
 * If you need write access to time_bin, use part_get_time_bin_p() instead.
 */
static __attribute__((always_inline)) INLINE const timebin_t*
  part_get_const_time_bin_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_time_bin;
}


/**
 * @brief set the value of time_bin, Time-step length.
 */
static __attribute__((always_inline)) INLINE void
  part_set_time_bin(const size_t pind, const timebin_t time_bin) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_time_bin = time_bin;
}





/**
 * @brief get limiter_data, Time-step limiter information.
 */
static __attribute__((always_inline)) INLINE struct timestep_limiter_data
  part_get_limiter_data(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_limiter_data;
}


/**
 * @brief get a pointer to limiter_data, Time-step limiter information.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to limiter_data. If you need read-only access to limiter_data, use part_get_const_limiter_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct timestep_limiter_data*
  part_get_limiter_data_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_limiter_data;
}


/**
 * @brief get read-only access to pointer to limiter_data,
 * Time-step limiter information.
 * If you need write access to limiter_data, use part_get_limiter_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct timestep_limiter_data*
  part_get_const_limiter_data_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_limiter_data;
}


/**
 * @brief set the value of limiter_data, Time-step limiter information.
 */
static __attribute__((always_inline)) INLINE void
  part_set_limiter_data(const size_t pind, const struct timestep_limiter_data limiter_data) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_limiter_data = limiter_data;
}





/**
 * @brief get ti_drift, Time of the last drift.
 */
static __attribute__((always_inline)) INLINE integertime_t
  part_get_ti_drift(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return part_s->_ti_drift;
#else
  return LLONG_MAX;
#endif
}


/**
 * @brief get a pointer to ti_drift, Time of the last drift.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to ti_drift. If you need read-only access to ti_drift, use part_get_const_ti_drift_p() instead.
 */
static __attribute__((always_inline)) INLINE integertime_t*
  part_get_ti_drift_p(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return &part_s->_ti_drift;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to ti_drift,
 * Time of the last drift.
 * If you need write access to ti_drift, use part_get_ti_drift_p() instead.
 */
static __attribute__((always_inline)) INLINE const integertime_t*
  part_get_const_ti_drift_p(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return &part_s->_ti_drift;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of ti_drift, Time of the last drift.
 */
static __attribute__((always_inline)) INLINE void
  part_set_ti_drift(const size_t pind, const integertime_t ti_drift) {
#ifdef SWIFT_DEBUG_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  part_s->_ti_drift = ti_drift;
#endif
}





/**
 * @brief get ti_kick, Time of the last drift.
 */
static __attribute__((always_inline)) INLINE integertime_t
  part_get_ti_kick(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return part_s->_ti_kick;
#else
  return LLONG_MAX;
#endif
}


/**
 * @brief get a pointer to ti_kick, Time of the last drift.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to ti_kick. If you need read-only access to ti_kick, use part_get_const_ti_kick_p() instead.
 */
static __attribute__((always_inline)) INLINE integertime_t*
  part_get_ti_kick_p(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return &part_s->_ti_kick;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to ti_kick,
 * Time of the last drift.
 * If you need write access to ti_kick, use part_get_ti_kick_p() instead.
 */
static __attribute__((always_inline)) INLINE const integertime_t*
  part_get_const_ti_kick_p(const size_t pind) {
#ifdef SWIFT_DEBUG_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  return &part_s->_ti_kick;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of ti_kick, Time of the last drift.
 */
static __attribute__((always_inline)) INLINE void
  part_set_ti_kick(const size_t pind, const integertime_t ti_kick) {
#ifdef SWIFT_DEBUG_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;

  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);

  part_s->_ti_kick = ti_kick;
#endif
}





/**
 * @brief get N_density, Integer number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_density(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_density;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_density, Integer number of neighbours in the density loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_density. If you need read-only access to N_density, use part_get_const_N_density_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_density_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_density;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_density,
 * Integer number of neighbours in the density loop.
 * If you need write access to N_density, use part_get_N_density_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_density_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_density;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_density, Integer number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_density(const size_t pind, const int N_density) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_density = N_density;
#endif
}





/**
 * @brief get N_density_exact, Exact integer number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_density_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_density_exact;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_density_exact, Exact integer number of neighbours in the density loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_density_exact. If you need read-only access to N_density_exact, use part_get_const_N_density_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_density_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_density_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_density_exact,
 * Exact integer number of neighbours in the density loop.
 * If you need write access to N_density_exact, use part_get_N_density_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_density_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_density_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_density_exact, Exact integer number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_density_exact(const size_t pind, const int N_density_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_density_exact = N_density_exact;
#endif
}





/**
 * @brief get N_gradient, Integer number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_gradient(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_gradient;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_gradient, Integer number of neighbours in the gradient loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_gradient. If you need read-only access to N_gradient, use part_get_const_N_gradient_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_gradient_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_gradient;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_gradient,
 * Integer number of neighbours in the gradient loop.
 * If you need write access to N_gradient, use part_get_N_gradient_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_gradient_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_gradient;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_gradient, Integer number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_gradient(const size_t pind, const int N_gradient) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_gradient = N_gradient;
#endif
}





/**
 * @brief get N_gradient_exact, Exact integer number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_gradient_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_gradient_exact;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_gradient_exact, Exact integer number of neighbours in the gradient loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_gradient_exact. If you need read-only access to N_gradient_exact, use part_get_const_N_gradient_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_gradient_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_gradient_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_gradient_exact,
 * Exact integer number of neighbours in the gradient loop.
 * If you need write access to N_gradient_exact, use part_get_N_gradient_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_gradient_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_gradient_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_gradient_exact, Exact integer number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_gradient_exact(const size_t pind, const int N_gradient_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_gradient_exact = N_gradient_exact;
#endif
}





/**
 * @brief get N_force, Integer number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_force(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_force;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_force, Integer number of neighbours in the force loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_force. If you need read-only access to N_force, use part_get_const_N_force_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_force_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_force;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_force,
 * Integer number of neighbours in the force loop.
 * If you need write access to N_force, use part_get_N_force_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_force_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_force;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_force, Integer number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_force(const size_t pind, const int N_force) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_force = N_force;
#endif
}





/**
 * @brief get N_force_exact, Exact integer number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE int
  part_get_N_force_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_N_force_exact;
#else
  return INT_MAX;
#endif
}


/**
 * @brief get a pointer to N_force_exact, Exact integer number of neighbours in the force loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to N_force_exact. If you need read-only access to N_force_exact, use part_get_const_N_force_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE int*
  part_get_N_force_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_force_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to N_force_exact,
 * Exact integer number of neighbours in the force loop.
 * If you need write access to N_force_exact, use part_get_N_force_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const int*
  part_get_const_N_force_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_N_force_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of N_force_exact, Exact integer number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_N_force_exact(const size_t pind, const int N_force_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_N_force_exact = N_force_exact;
#endif
}





/**
 * @brief get rho_exact, Exact value of the density field obtained via brute-force loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_rho_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_rho_exact;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to rho_exact, Exact value of the density field obtained via brute-force loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to rho_exact. If you need read-only access to rho_exact, use part_get_const_rho_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_rho_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rho_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to rho_exact,
 * Exact value of the density field obtained via brute-force loop.
 * If you need write access to rho_exact, use part_get_rho_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_rho_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_rho_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of rho_exact, Exact value of the density field obtained via brute-force loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_rho_exact(const size_t pind, const float rho_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_rho_exact = rho_exact;
#endif
}





/**
 * @brief get n_density, Weighted number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_density(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_density;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_density, Weighted number of neighbours in the density loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_density. If you need read-only access to n_density, use part_get_const_n_density_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_density_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_density;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_density,
 * Weighted number of neighbours in the density loop.
 * If you need write access to n_density, use part_get_n_density_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_density_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_density;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_density, Weighted number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_density(const size_t pind, const float n_density) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_density = n_density;
#endif
}





/**
 * @brief get n_density_exact, Exact value of the weighted number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_density_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_density_exact;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_density_exact, Exact value of the weighted number of neighbours in the density loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_density_exact. If you need read-only access to n_density_exact, use part_get_const_n_density_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_density_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_density_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_density_exact,
 * Exact value of the weighted number of neighbours in the density loop.
 * If you need write access to n_density_exact, use part_get_n_density_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_density_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_density_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_density_exact, Exact value of the weighted number of neighbours in the density loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_density_exact(const size_t pind, const float n_density_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_density_exact = n_density_exact;
#endif
}





/**
 * @brief get n_gradient, Weighted number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_gradient(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_gradient;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_gradient, Weighted number of neighbours in the gradient loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_gradient. If you need read-only access to n_gradient, use part_get_const_n_gradient_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_gradient_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_gradient;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_gradient,
 * Weighted number of neighbours in the gradient loop.
 * If you need write access to n_gradient, use part_get_n_gradient_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_gradient_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_gradient;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_gradient, Weighted number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_gradient(const size_t pind, const float n_gradient) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_gradient = n_gradient;
#endif
}





/**
 * @brief get n_gradient_exact, Exact value of the weighted number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_gradient_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_gradient_exact;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_gradient_exact, Exact value of the weighted number of neighbours in the gradient loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_gradient_exact. If you need read-only access to n_gradient_exact, use part_get_const_n_gradient_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_gradient_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_gradient_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_gradient_exact,
 * Exact value of the weighted number of neighbours in the gradient loop.
 * If you need write access to n_gradient_exact, use part_get_n_gradient_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_gradient_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_gradient_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_gradient_exact, Exact value of the weighted number of neighbours in the gradient loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_gradient_exact(const size_t pind, const float n_gradient_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_gradient_exact = n_gradient_exact;
#endif
}





/**
 * @brief get n_force, Weighted number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_force(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_force;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_force, Weighted number of neighbours in the force loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_force. If you need read-only access to n_force, use part_get_const_n_force_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_force_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_force;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_force,
 * Weighted number of neighbours in the force loop.
 * If you need write access to n_force, use part_get_n_force_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_force_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_force;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_force, Weighted number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_force(const size_t pind, const float n_force) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_force = n_force;
#endif
}





/**
 * @brief get n_force_exact, Exact value of the weighted number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE float
  part_get_n_force_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_n_force_exact;
#else
  return FLT_MAX;
#endif
}


/**
 * @brief get a pointer to n_force_exact, Exact value of the weighted number of neighbours in the force loop.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to n_force_exact. If you need read-only access to n_force_exact, use part_get_const_n_force_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_n_force_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_force_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to n_force_exact,
 * Exact value of the weighted number of neighbours in the force loop.
 * If you need write access to n_force_exact, use part_get_n_force_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_n_force_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_n_force_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of n_force_exact, Exact value of the weighted number of neighbours in the force loop.
 */
static __attribute__((always_inline)) INLINE void
  part_set_n_force_exact(const size_t pind, const float n_force_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_n_force_exact = n_force_exact;
#endif
}





/**
 * @brief get inhibited_exact, Has this particle interacted with any unhibited neighbour?.
 */
static __attribute__((always_inline)) INLINE char
  part_get_inhibited_exact(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_inhibited_exact;
#else
  return CHAR_MAX;
#endif
}


/**
 * @brief get a pointer to inhibited_exact, Has this particle interacted with any unhibited neighbour?.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to inhibited_exact. If you need read-only access to inhibited_exact, use part_get_const_inhibited_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE char*
  part_get_inhibited_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_inhibited_exact;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to inhibited_exact,
 * Has this particle interacted with any unhibited neighbour?.
 * If you need write access to inhibited_exact, use part_get_inhibited_exact_p() instead.
 */
static __attribute__((always_inline)) INLINE const char*
  part_get_const_inhibited_exact_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_inhibited_exact;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of inhibited_exact, Has this particle interacted with any unhibited neighbour?.
 */
static __attribute__((always_inline)) INLINE void
  part_set_inhibited_exact(const size_t pind, const char inhibited_exact) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_inhibited_exact = inhibited_exact;
#endif
}





/**
 * @brief get limited_part, Has this particle been woken up by the limiter?.
 */
static __attribute__((always_inline)) INLINE char
  part_get_limited_part(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_limited_part;
#else
  return CHAR_MAX;
#endif
}


/**
 * @brief get a pointer to limited_part, Has this particle been woken up by the limiter?.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to limited_part. If you need read-only access to limited_part, use part_get_const_limited_part_p() instead.
 */
static __attribute__((always_inline)) INLINE char*
  part_get_limited_part_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_limited_part;
#else
  return NULL;
#endif
}


/**
 * @brief get read-only access to pointer to limited_part,
 * Has this particle been woken up by the limiter?.
 * If you need write access to limited_part, use part_get_limited_part_p() instead.
 */
static __attribute__((always_inline)) INLINE const char*
  part_get_const_limited_part_p(const size_t pind) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_limited_part;
#else
  return NULL;
#endif
}


/**
 * @brief set the value of limited_part, Has this particle been woken up by the limiter?.
 */
static __attribute__((always_inline)) INLINE void
  part_set_limited_part(const size_t pind, const char limited_part) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_limited_part = limited_part;
#endif
}





/**
 * @brief get geometry, Geometrical quantities used for Finite Volume Particle Method RT.
 */
static __attribute__((always_inline)) INLINE struct fvpm_geometry_struct
  part_get_geometry(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return part_s->_geometry;
}


/**
 * @brief get a pointer to geometry, Geometrical quantities used for Finite Volume Particle Method RT.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to geometry. If you need read-only access to geometry, use part_get_const_geometry_p() instead.
 */
static __attribute__((always_inline)) INLINE struct fvpm_geometry_struct*
  part_get_geometry_p(const size_t pind) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_geometry;
}


/**
 * @brief get read-only access to pointer to geometry,
 * Geometrical quantities used for Finite Volume Particle Method RT.
 * If you need write access to geometry, use part_get_geometry_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct fvpm_geometry_struct*
  part_get_const_geometry_p(const size_t pind) {

  const struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  return &part_s->_geometry;
}


/**
 * @brief set the value of geometry, Geometrical quantities used for Finite Volume Particle Method RT.
 */
static __attribute__((always_inline)) INLINE void
  part_set_geometry(const size_t pind, const struct fvpm_geometry_struct geometry) {

  struct part* restrict part_s = global_part_arrays._part + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(part_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", part_s->_accessor_id, p->_accessor_id);
#endif
  part_s->_geometry = geometry;
}







/**
 * @brief get x_diff, Offset between current position and position at last tree rebuild.,
 * for read and write access. For read-only access, use
 * part_get_const_x_diff() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_x_diff(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff;
}

/**
 * @brief get x_diff, Offset between current position and position at last tree rebuild., for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_x_diff(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff;
}

/**
 * @brief get x_diff, Offset between current position and position at last tree rebuild., by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_x_diff_ind(const size_t pind, const int i) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff[i];
}

/**
 * @brief set all values of x_diff, Offset between current position and position at last tree rebuild.,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x_diff(const size_t pind, const float x_diff[3]) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_x_diff[0] = x_diff[0];
  xpart_s->_x_diff[1] = x_diff[1];
  xpart_s->_x_diff[2] = x_diff[2];
}

/**
 * @brief set the value of x_diff, Offset between current position and position at last tree rebuild., by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x_diff_ind(const size_t pind, const int i, const float x_diff) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_x_diff[i] = x_diff;
}




/**
 * @brief get x_diff_sort, Offset between the current position and position at the last sort.,
 * for read and write access. For read-only access, use
 * part_get_const_x_diff_sort() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_x_diff_sort(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff_sort;
}

/**
 * @brief get x_diff_sort, Offset between the current position and position at the last sort., for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_x_diff_sort(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff_sort;
}

/**
 * @brief get x_diff_sort, Offset between the current position and position at the last sort., by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_x_diff_sort_ind(const size_t pind, const int i) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_x_diff_sort[i];
}

/**
 * @brief set all values of x_diff_sort, Offset between the current position and position at the last sort.,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x_diff_sort(const size_t pind, const float x_diff_sort[3]) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_x_diff_sort[0] = x_diff_sort[0];
  xpart_s->_x_diff_sort[1] = x_diff_sort[1];
  xpart_s->_x_diff_sort[2] = x_diff_sort[2];
}

/**
 * @brief set the value of x_diff_sort, Offset between the current position and position at the last sort., by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_x_diff_sort_ind(const size_t pind, const int i, const float x_diff_sort) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_x_diff_sort[i] = x_diff_sort;
}




/**
 * @brief get v_full, Velocity at the last full step.,
 * for read and write access. For read-only access, use
 * part_get_const_v_full() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_v_full(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_v_full;
}

/**
 * @brief get v_full, Velocity at the last full step., for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_v_full(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_v_full;
}

/**
 * @brief get v_full, Velocity at the last full step., by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_v_full_ind(const size_t pind, const int i) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_v_full[i];
}

/**
 * @brief set all values of v_full, Velocity at the last full step.,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_v_full(const size_t pind, const float v_full[3]) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_v_full[0] = v_full[0];
  xpart_s->_v_full[1] = v_full[1];
  xpart_s->_v_full[2] = v_full[2];
}

/**
 * @brief set the value of v_full, Velocity at the last full step., by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_v_full_ind(const size_t pind, const int i, const float v_full) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_v_full[i] = v_full;
}




/**
 * @brief get a_grav, Gravitational acceleration at the end of the last step,
 * for read and write access. For read-only access, use
 * part_get_const_a_grav() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_a_grav(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_a_grav;
}

/**
 * @brief get a_grav, Gravitational acceleration at the end of the last step, for read-only access.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_a_grav(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_a_grav;
}

/**
 * @brief get a_grav, Gravitational acceleration at the end of the last step, by index.
 */
static __attribute__((always_inline)) INLINE float
  part_get_a_grav_ind(const size_t pind, const int i) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_a_grav[i];
}

/**
 * @brief set all values of a_grav, Gravitational acceleration at the end of the last step,
 * from an array.
 */
static __attribute__((always_inline)) INLINE void
  part_set_a_grav(const size_t pind, const float a_grav[3]) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_a_grav[0] = a_grav[0];
  xpart_s->_a_grav[1] = a_grav[1];
  xpart_s->_a_grav[2] = a_grav[2];
}

/**
 * @brief set the value of a_grav, Gravitational acceleration at the end of the last step, by index i.
 */
static __attribute__((always_inline)) INLINE void
  part_set_a_grav_ind(const size_t pind, const int i, const float a_grav) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_a_grav[i] = a_grav;
}




/**
 * @brief get u_full, Internal energy at the last full step..
 */
static __attribute__((always_inline)) INLINE float
  part_get_u_full(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_u_full;
}


/**
 * @brief get a pointer to u_full, Internal energy at the last full step..
 * Use this only if you need to modify the value, i.e. if you need write access
 * to u_full. If you need read-only access to u_full, use part_get_const_u_full_p() instead.
 */
static __attribute__((always_inline)) INLINE float*
  part_get_u_full_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_u_full;
}


/**
 * @brief get read-only access to pointer to u_full,
 * Internal energy at the last full step..
 * If you need write access to u_full, use part_get_u_full_p() instead.
 */
static __attribute__((always_inline)) INLINE const float*
  part_get_const_u_full_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_u_full;
}


/**
 * @brief set the value of u_full, Internal energy at the last full step..
 */
static __attribute__((always_inline)) INLINE void
  part_set_u_full(const size_t pind, const float u_full) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_u_full = u_full;
}





/**
 * @brief get split_data, Additional data used to record particle splits.
 */
static __attribute__((always_inline)) INLINE struct particle_splitting_data
  part_get_split_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_split_data;
}


/**
 * @brief get a pointer to split_data, Additional data used to record particle splits.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to split_data. If you need read-only access to split_data, use part_get_const_split_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct particle_splitting_data*
  part_get_split_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_split_data;
}


/**
 * @brief get read-only access to pointer to split_data,
 * Additional data used to record particle splits.
 * If you need write access to split_data, use part_get_split_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct particle_splitting_data*
  part_get_const_split_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_split_data;
}


/**
 * @brief set the value of split_data, Additional data used to record particle splits.
 */
static __attribute__((always_inline)) INLINE void
  part_set_split_data(const size_t pind, const struct particle_splitting_data split_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_split_data = split_data;
}





/**
 * @brief get extra_cooling_data, Additional data used to record cooling information.
 */
static __attribute__((always_inline)) INLINE struct cooling_xpart_data
  part_get_extra_cooling_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_extra_cooling_data;
}


/**
 * @brief get a pointer to extra_cooling_data, Additional data used to record cooling information.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to extra_cooling_data. If you need read-only access to extra_cooling_data, use part_get_const_extra_cooling_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct cooling_xpart_data*
  part_get_extra_cooling_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_cooling_data;
}


/**
 * @brief get read-only access to pointer to extra_cooling_data,
 * Additional data used to record cooling information.
 * If you need write access to extra_cooling_data, use part_get_extra_cooling_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct cooling_xpart_data*
  part_get_const_extra_cooling_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_cooling_data;
}


/**
 * @brief set the value of extra_cooling_data, Additional data used to record cooling information.
 */
static __attribute__((always_inline)) INLINE void
  part_set_extra_cooling_data(const size_t pind, const struct cooling_xpart_data extra_cooling_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_extra_cooling_data = extra_cooling_data;
}





/**
 * @brief get extra_tracers_data, Additional data used by the tracers.
 */
static __attribute__((always_inline)) INLINE struct tracers_xpart_data
  part_get_extra_tracers_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_extra_tracers_data;
}


/**
 * @brief get a pointer to extra_tracers_data, Additional data used by the tracers.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to extra_tracers_data. If you need read-only access to extra_tracers_data, use part_get_const_extra_tracers_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct tracers_xpart_data*
  part_get_extra_tracers_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_tracers_data;
}


/**
 * @brief get read-only access to pointer to extra_tracers_data,
 * Additional data used by the tracers.
 * If you need write access to extra_tracers_data, use part_get_extra_tracers_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct tracers_xpart_data*
  part_get_const_extra_tracers_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_tracers_data;
}


/**
 * @brief set the value of extra_tracers_data, Additional data used by the tracers.
 */
static __attribute__((always_inline)) INLINE void
  part_set_extra_tracers_data(const size_t pind, const struct tracers_xpart_data extra_tracers_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_extra_tracers_data = extra_tracers_data;
}





/**
 * @brief get extra_sf_data, Additional data used by the tracers.
 */
static __attribute__((always_inline)) INLINE struct star_formation_xpart_data
  part_get_extra_sf_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_extra_sf_data;
}


/**
 * @brief get a pointer to extra_sf_data, Additional data used by the tracers.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to extra_sf_data. If you need read-only access to extra_sf_data, use part_get_const_extra_sf_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct star_formation_xpart_data*
  part_get_extra_sf_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_sf_data;
}


/**
 * @brief get read-only access to pointer to extra_sf_data,
 * Additional data used by the tracers.
 * If you need write access to extra_sf_data, use part_get_extra_sf_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct star_formation_xpart_data*
  part_get_const_extra_sf_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_sf_data;
}


/**
 * @brief set the value of extra_sf_data, Additional data used by the tracers.
 */
static __attribute__((always_inline)) INLINE void
  part_set_extra_sf_data(const size_t pind, const struct star_formation_xpart_data extra_sf_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_extra_sf_data = extra_sf_data;
}





/**
 * @brief get extra_feedback_data, Additional data used by the feedback.
 */
static __attribute__((always_inline)) INLINE struct feedback_xpart_data
  part_get_extra_feedback_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_extra_feedback_data;
}


/**
 * @brief get a pointer to extra_feedback_data, Additional data used by the feedback.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to extra_feedback_data. If you need read-only access to extra_feedback_data, use part_get_const_extra_feedback_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct feedback_xpart_data*
  part_get_extra_feedback_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_feedback_data;
}


/**
 * @brief get read-only access to pointer to extra_feedback_data,
 * Additional data used by the feedback.
 * If you need write access to extra_feedback_data, use part_get_extra_feedback_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct feedback_xpart_data*
  part_get_const_extra_feedback_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_feedback_data;
}


/**
 * @brief set the value of extra_feedback_data, Additional data used by the feedback.
 */
static __attribute__((always_inline)) INLINE void
  part_set_extra_feedback_data(const size_t pind, const struct feedback_xpart_data extra_feedback_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_extra_feedback_data = extra_feedback_data;
}





/**
 * @brief get extra_mhd_data, Additional data used by the MHD scheme.
 */
static __attribute__((always_inline)) INLINE struct mhd_xpart_data
  part_get_extra_mhd_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_extra_mhd_data;
}


/**
 * @brief get a pointer to extra_mhd_data, Additional data used by the MHD scheme.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to extra_mhd_data. If you need read-only access to extra_mhd_data, use part_get_const_extra_mhd_data_p() instead.
 */
static __attribute__((always_inline)) INLINE struct mhd_xpart_data*
  part_get_extra_mhd_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_mhd_data;
}


/**
 * @brief get read-only access to pointer to extra_mhd_data,
 * Additional data used by the MHD scheme.
 * If you need write access to extra_mhd_data, use part_get_extra_mhd_data_p() instead.
 */
static __attribute__((always_inline)) INLINE const struct mhd_xpart_data*
  part_get_const_extra_mhd_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_extra_mhd_data;
}


/**
 * @brief set the value of extra_mhd_data, Additional data used by the MHD scheme.
 */
static __attribute__((always_inline)) INLINE void
  part_set_extra_mhd_data(const size_t pind, const struct mhd_xpart_data extra_mhd_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_extra_mhd_data = extra_mhd_data;
}





/**
 * @brief get csds_data, Additional data for the particle csds.
 */
#ifdef WITH_CSDS
static __attribute__((always_inline)) INLINE struct csds_part_data
  part_get_csds_data(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return xpart_s->_csds_data;
}
#endif

/**
 * @brief get a pointer to csds_data, Additional data for the particle csds.
 * Use this only if you need to modify the value, i.e. if you need write access
 * to csds_data. If you need read-only access to csds_data, use part_get_const_csds_data_p() instead.
 */
#ifdef WITH_CSDS
static __attribute__((always_inline)) INLINE struct csds_part_data*
  part_get_csds_data_p(const size_t pind) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_csds_data;
}
#endif

/**
 * @brief get read-only access to pointer to csds_data,
 * Additional data for the particle csds.
 * If you need write access to csds_data, use part_get_csds_data_p() instead.
 */
#ifdef WITH_CSDS
static __attribute__((always_inline)) INLINE const struct csds_part_data*
  part_get_const_csds_data_p(const size_t pind) {

  const struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  return &xpart_s->_csds_data;
}
#endif

/**
 * @brief set the value of csds_data, Additional data for the particle csds.
 */
#ifdef WITH_CSDS
static __attribute__((always_inline)) INLINE void
  part_set_csds_data(const size_t pind, const struct csds_part_data csds_data) {

  struct xpart* restrict xpart_s = global_part_arrays._xpart + pind;
#ifdef SWIFT_DEBUG_CHECKS
  /* Forbid ID = 0 to prevent false positives by forgotten initialisation */
  const struct part* restrict p = global_part_arrays._part + pind;
  swift_assert(p->_accessor_id != 0);
  /* Make sure we're accessing the correct data */
  if(xpart_s->_accessor_id != p->_accessor_id)
    error("Accessor IDs not equal: %lld %lld", xpart_s->_accessor_id, p->_accessor_id);
#endif
  xpart_s->_csds_data = csds_data;
}
#endif







#endif /* SWIFT_HYDRO_PART_SPHENIX_H */
