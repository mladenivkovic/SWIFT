/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

/* This object's header. */
#include "tools.h"

/* Local includes. */
#include "active.h"
#include "adaptive_softening.h"
#include "cell.h"
#include "chemistry.h"
#include "cosmology.h"
#include "error.h"
#include "feedback.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_data.h"
#include "mhd.h"
#include "part.h"
#include "periodic.h"
#include "pressure_floor_iact.h"
#include "runner.h"
#include "sink_iact.h"
#include "sink_properties.h"
#include "space_getsid.h"
#include "star_formation_iact.h"
#include "stars.h"

/**
 *  Factorize a given integer, attempts to keep larger pair of factors.
 */
void factor(int value, int *f1, int *f2) {

  int j = (int)sqrt(value);
  for (int i = j; i > 0; i--) {
    if ((value % i) == 0) {
      *f1 = i;
      *f2 = value / i;
      break;
    }
  }
}

/**
 * @brief Compute the average number of pairs per particle using
 *      a brute-force O(N^2) computation.
 *
 * @param dim The space dimensions.
 * @param N The number of parts.
 * @param periodic Periodic boundary conditions flag.
 */
void pairs_n2(double *dim, size_t first, size_t N, int periodic) {
  size_t i;
  size_t j;
  size_t k;
  int count = 0;
  // int mj, mk;
  // double maxratio = 1.0;
  double r2;
  double rho = 0.0;
  double rho_max = 0.0;
  double rho_min = 100;
  float a = 1.f;
  float H = 0.f;
  float dx[3];

  /* Loop over all particle pairs. */
  for (j = first; j < first + N; j++) {
    if (j % 1000 == 0) {
      printf("pairs_n2: j=%zu.\n", j);
      fflush(stdout);
    }
    const double *const xj = part_get_const_x(j);
    const float hj = part_get_h(j);
    for (k = j + 1; k < first + N; k++) {
      const double *const xk = part_get_const_x(k);
      const float hk = part_get_h(k);
      for (i = 0; i < 3; i++) {
        dx[i] = xj[i] - xk[i];
        if (periodic) {
          if (dx[i] < -dim[i] / 2)
            dx[i] += dim[i];
          else if (dx[i] > dim[i] / 2)
            dx[i] -= dim[i];
        }
      }
      r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (r2 < hj * hj || r2 < hk * hk) {
        runner_iact_density(r2, dx, hj, hk, j, k, a, H);
        /* if ( hj / hk > maxratio )
            {
            maxratio = hj / hk;
            mj = j; mk = k;
            }
        else if ( hk / hj > maxratio )
            {
            maxratio = hk / hj;
            mj = j; mk = k;
            } */
      }
    }
  }

  /* Aggregate the results. */
  for (k = first; k < first + N; k++) {
    // count += parts[k].icount;
    const float wcount = part_get_wcount(k);
    rho += wcount;
    rho_min = fmin(wcount, rho_min);
    rho_max = fmax(wcount, rho_max);
  }

  /* Dump the result. */
  printf("pairs_n2: avg. density per part is %.3f (nr. pairs %.3f).\n",
         rho / N + 32.0 / 3, ((double)count) / N);
  printf("pairs_n2: densities are in [ %e , %e ].\n", rho_min / N + 32.0 / 3,
         rho_max / N + 32.0 / 3);
  /* printf( "pairs_n2: maximum ratio between parts %i [%e,%e,%e] and %i
     [%e,%e,%e] is %.3f/%.3f\n" ,
      mj , parts[mj].x[0] , parts[mj].x[1] , parts[mj].x[2] ,
      mk , parts[mk].x[0] , parts[mk].x[1] , parts[mk].x[2] ,
      parts[mj].h , parts[mk].h ); fflush(stdout); */
  fflush(stdout);
}

void pairs_single_density(double *dim, long long int pid, size_t pind, size_t first, size_t N, int periodic) {
  const float a = 1.f;
  const float H = 0.f;

  /* Find "our" part. */
  size_t k;
  for (k = first; k < first + N && part_get_id(k) != pid; k++) {
    /* Nothing to do here */
  }

  /* Clear accumulators. */
  if (k == first + N) error("Part not found.");
  printf("pairs_single: part[%zu].id == %lli.\n", k, pid);

  hydro_init_part(pind, NULL);
  adaptive_softening_init_part(pind);
  mhd_init_part(pind);
  const double *xp = part_get_const_x(pind);
  const float hp = part_get_h(pind);

  /* Loop over all particle pairs. */
  for (k = 0; k < N; k++) {
    if (part_get_id(k) == part_get_id(pind)) continue;

    double dx[3];
    float fdx[3];

    const double *xk = part_get_const_x(k);
    for (int i = 0; i < 3; i++) {
      dx[i] = xp[i] - xk[i];
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    const float r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    if (r2 < hp * hp) {
      runner_iact_nonsym_density(r2, fdx, hp, part_get_h(k), pind, k, a, H);
      /* printf( "pairs_simple: interacting particles %lli [%i,%i,%i] and %lli
         [%i,%i,%i], r=%e.\n" ,
          pid , (int)(xp[0]*ih) , (int)(xp[1]*ih) , (int)(xp[2]*ih) ,
          part_get_id(&parts[k]) , (int)(part_get_x_ind(&parts[k], 0)*ih),
          (int)(part_get_x_ind(&parts[k], 1)*ih),
          (int)(part_get_x_ind(&parts[k], 2)*ih),
          sqrtf(r2) ); */
    }
  }

  /* Dump the result. */
  printf("pairs_single: wcount of part %lli (h=%e) is %f.\n", part_get_id(pind),
         part_get_h(pind), part_get_wcount(pind) + 32.0 / 3);
  fflush(stdout);
}

void pairs_all_density(struct runner *r, struct cell *ci, struct cell *cj) {

  size_t pindi;
  size_t pindj;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;
  double shift[3] = {0.0, 0.0, 0.0};
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);
  const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                             cj->loc[2] + shift[2]};
  const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    const float hi = part_get_h(pindi);
    const float hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pindi, e)) continue;

    const double *xi = part_get_const_x(pindi);
    const float pix = xi[0] - shift_i[0];
    const float piy = xi[1] - shift_i[1];
    const float piz = xi[2] - shift_i[2];

    for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

      const double *xj = part_get_const_x(pindj);
      const float pjx = xj[0] - shift_j[0];
      const float pjy = xj[1] - shift_j[1];
      const float pjz = xj[2] - shift_j[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pix - pjx, dim[0]),
                           nearest(piy - pjy, dim[1]),
                           nearest(piz - pjz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      const float hj = part_get_h(pindj);

      /* Hit or miss? */
      if (r2 < hig2 && !part_is_inhibited(j, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, i, j, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

    const float hj = part_get_h(pindj);
    const float hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(pindj, e)) continue;

    const double *xj = part_get_const_x(pindj);
    const float pjx = xj[0] - shift_j[0];
    const float pjy = xj[1] - shift_j[1];
    const float pjz = xj[2] - shift_j[2];

    for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

      const double *xi = part_get_const_x(pindi);
      const float pix = xi[0] - shift_i[0];
      const float piy = xi[1] - shift_i[1];
      const float piz = xi[2] - shift_i[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pjx - pix, dim[0]),
                           nearest(pjy - piy, dim[1]),
                           nearest(pjz - piz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      const float hi = part_get_h(pindi);

      /* Hit or miss? */
      if (r2 < hjg2 && !part_is_inhibited(i, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, j, i, a, H);
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP
void pairs_all_gradient(struct runner *r, struct cell *ci, struct cell *cj) {

  // Keepig hig2, hjg2 here to disable -Wunused-variable warnings
  float hig2 = 0.f;
  float hjg2 = 0.f;
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;
  double shift[3] = {0.0, 0.0, 0.0};
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);
  const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                             cj->loc[2] + shift[2]};
  const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    const float hi = part_get_h(i);
    hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(i, e)) continue;

    const double *xi = part_get_const_x(i);
    const float pix = xi[0] - shift_i[0];
    const float piy = xi[1] - shift_i[1];
    const float piz = xi[2] - shift_i[2];

    for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

      const float hj = part_get_h(j);
      hjg2 = hj * hj * kernel_gamma2;

      const double *xj = part_get_const_x(j);
      const float pjx = xj[0] - shift_j[0];
      const float pjy = xj[1] - shift_j[1];
      const float pjz = xj[2] - shift_j[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pix - pjx, dim[0]),
                           nearest(piy - pjy, dim[1]),
                           nearest(piz - pjz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (!part_is_inhibited(j, e)) {
#ifdef EXTRA_HYDRO_LOOP_TYPE2
        if (r2 < hig2 || r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hi, hj, i, j, a, H);
        }
#else
        if (r2 < hig2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hi, hj, i, j, a, H);
        }
#endif /* EXTRA_HYDRO_LOOP_TYPE2 */
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

    const float hj = part_get_h(j);
    hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(j, e)) continue;

    const double *xj = part_get_const_x(j);
    const float pjx = xj[0] - shift_j[0];
    const float pjy = xj[1] - shift_j[1];
    const float pjz = xj[2] - shift_j[2];

    for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

      const float hi = part_get_h(i);
      hig2 = hi * hi * kernel_gamma2;

      const double *xi = part_get_const_x(i);
      const float pix = xi[0] - shift_i[0];
      const float piy = xi[1] - shift_i[1];
      const float piz = xi[2] - shift_i[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pjx - pix, dim[0]),
                           nearest(pjy - piy, dim[1]),
                           nearest(pjz - piz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (!part_is_inhibited(i, e)) {
#ifdef EXTRA_HYDRO_LOOP_TYPE2
        if (r2 < hig2 || r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hj, hi, j, i, a, H);
        }
#else
        if (r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hj, hi, j, i, a, H);
        }
#endif /* EXTRA_HYDRO_LOOP_TYPE2 */
      }
    }
  }
}
#endif /* EXTRA_HDYRO_LOOP */

void pairs_all_force(struct runner *r, struct cell *ci, struct cell *cj) {

  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;
  double shift[3] = {0.0, 0.0, 0.0};
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);
  const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                             cj->loc[2] + shift[2]};
  const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    const float hi = part_get_h(i);
    const float hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(i, e)) continue;

    const double *xi = part_get_const_x(i);
    const float pix = xi[0] - shift_i[0];
    const float piy = xi[1] - shift_i[1];
    const float piz = xi[2] - shift_i[2];

    for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

      const float hj = part_get_h(j);
      const float hjg2 = hj * hj * kernel_gamma2;

      const double *xj = part_get_const_x(j);
      const float pjx = xj[0] - shift_j[0];
      const float pjy = xj[1] - shift_j[1];
      const float pjz = xj[2] - shift_j[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pix - pjx, dim[0]),
                           nearest(piy - pjy, dim[1]),
                           nearest(piz - pjz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hi, hj, i, j, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

    const float hj = part_get_h(j);
    const float hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!part_is_active(j, e)) continue;

    const double *xj = part_get_const_x(j);
    const float pjx = xj[0] - shift_j[0];
    const float pjy = xj[1] - shift_j[1];
    const float pjz = xj[2] - shift_j[2];

    for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

      const float hi = part_get_h(i);
      const float hig2 = hi * hi * kernel_gamma2;

      const double *xi = part_get_const_x(i);
      const float pix = xi[0] - shift_i[0];
      const float piy = xi[1] - shift_i[1];
      const float piz = xi[2] - shift_i[2];

      /* Pairwise distance */
      const float dx[3] = {nearest(pjx - pix, dim[0]),
                           nearest(pjy - piy, dim[1]),
                           nearest(pjz - piz, dim[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hjg2 || r2 < hig2) {

        /* Interact */
        runner_iact_nonsym_force(r2, dx, hj, hi, j, i, a, H);
      }
    }
  }
}

void pairs_all_stars_density(struct runner *r, struct cell *ci,
                             struct cell *cj) {

  float r2;
  float dx[3];
  const double dim[3] = {r->e->s->dim[0], r->e->s->dim[1], r->e->s->dim[2]};
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->stars.count; ++i) {
    struct spart *spi = &ci->stars.parts[i];

    float hi = spi->h;
    float hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles. */
    if (!spart_is_active(spi, e)) continue;
    if (!feedback_is_active(spi, e)) continue;

    for (size_t j = cj->hydro.first; j < cj->hydro.first + cj->hydro.count; ++j) {

      /* Early abort? */
      if (part_is_inhibited(j, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = spi->x[k] - part_get_x_ind(j, k);
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dx, hi, part_get_h(j), spi, j, a, H);
      }
    }
  }

  /* Reverse double-for loop and checks every interaction */
  for (int j = 0; j < cj->stars.count; ++j) {

    struct spart *spj = &cj->stars.parts[j];
    float hj = spj->h;
    float hjg2 = hj * hj * kernel_gamma2;

    /* Skip inactive particles. */
    if (!spart_is_active(spj, e)) continue;
    if (!feedback_is_active(spj, e)) continue;

    for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

      /* Early abort? */
      if (part_is_inhibited(i, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dx[k] = spj->x[k] - part_get_x_ind(i, k);
        dx[k] = nearest(dx[k], dim[k]);
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dx, hj, part_get_h(i), spj, i, a, H);
      }
    }
  }
}

void self_all_density(struct runner *r, struct cell *ci) {
  float hi;
  float hj;
  float hig2;
  float hjg2;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    hi = part_get_h(i);
    hig2 = hi * hi * kernel_gamma2;

    const double *pix = part_get_const_x(i);

    for (size_t j = i + 1; j < ci->hydro.first + ci->hydro.count; ++j) {

      hj = part_get_h(j);
      hjg2 = hj * hj * kernel_gamma2;

      if (i == j) continue;

      const double *pjx = part_get_const_x(j);

      /* Pairwise distance */
      float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                     (float)(pix[2] - pjx[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hig2 && part_is_active(i, e) && !part_is_inhibited(j, e)) {

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hi, hj, i, j, a, H);
        runner_iact_nonsym_sink(r2, dx, hi, hj, i, j, a, H);
      }

      /* Hit or miss? */
      if (r2 < hjg2 && part_is_active(j, e) && !part_is_inhibited(i, e)) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        /* Interact */
        runner_iact_nonsym_density(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_pressure_floor(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_star_formation(r2, dx, hj, hi, j, i, a, H);
        runner_iact_nonsym_sink(r2, dx, hj, hi, j, i, a, H);
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP
void self_all_gradient(struct runner *r, struct cell *ci) {
  float hi;
  float hj;
  float hig2;
  float hjg2;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    hi = part_get_h(i);
    hig2 = hi * hi * kernel_gamma2;

    const double *pix = part_get_const_x(i);

    for (size_t j = i + 1; j < ci->hydro.first + ci->hydro.count; ++j) {

      hj = part_get_h(j);
      hjg2 = hj * hj * kernel_gamma2;

      if (i == j) continue;

      const double *pjx = part_get_const_x(j);

      /* Pairwise distance */
      float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                     (float)(pix[2] - pjx[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (part_is_active(i, e) && !part_is_inhibited(j, e)) {
#ifdef EXTRA_HYDRO_LOOP_TYPE2
        if (r2 < hig2 || r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hi, hj, i, j, a, H);
        }
#else
        if (r2 < hig2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hi, hj, i, j, a, H);
        }
#endif /* EXTRA_HYDRO_LOOP_TYPE2 */
      }

      /* Hit or miss? */
      if (part_is_active(j, e) && !part_is_inhibited(i, e)) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

#ifdef EXTRA_HYDRO_LOOP_TYPE2
        if (r2 < hig2 || r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hj, hi, j, i, a, H);
        }
#else
        if (r2 < hjg2) {
          /* Interact */
          runner_iact_nonsym_gradient(r2, dx, hj, hi, j, i, a, H);
        }
#endif /* EXTRA_HYDRO_LOOP_TYPE2 */
      }
    }
  }
}
#endif /* EXTRA_HYDRO_LOOP */

void self_all_force(struct runner *r, struct cell *ci) {
  float hi;
  float hj;
  float hig2;
  float hjg2;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (size_t i = ci->hydro.first; i < ci->hydro.first + ci->hydro.count; ++i) {

    hi = part_get_h(i);
    hig2 = hi * hi * kernel_gamma2;

    const double *pix = part_get_const_x(i);

    for (size_t j = i + 1; j < ci->hydro.first + ci->hydro.count; ++j) {

      hj = part_get_h(j);
      hjg2 = hj * hj * kernel_gamma2;

      if (i == j) continue;

      const double *pjx = part_get_const_x(j);

      /* Pairwise distance */
      float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                     (float)(pix[2] - pjx[2])};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hig2 || r2 < hjg2) {

        /* Interact */
        runner_iact_force(r2, dx, hi, hj, i, j, a, H);
      }
    }
  }
}

void self_all_stars_density(struct runner *r, struct cell *ci) {

  float r2;
  float hi;
  float hj;
  float hig2;
  float dxi[3];
  struct spart *spi;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Implements a double-for loop and checks every interaction */
  for (int i = 0; i < ci->stars.count; ++i) {

    spi = &ci->stars.parts[i];
    hi = spi->h;
    hig2 = hi * hi * kernel_gamma2;

    if (!spart_is_active(spi, e)) continue;
    if (!feedback_is_active(spi, e)) continue;

    for (size_t j = ci->hydro.first; j < ci->hydro.first + ci->hydro.count; ++j) {

      hj = part_get_h(j);

      /* Early abort? */
      if (part_is_inhibited(j, e)) continue;

      /* Pairwise distance */
      r2 = 0.0f;
      for (int k = 0; k < 3; k++) {
        dxi[k] = spi->x[k] - part_get_x_ind(j, k);
        r2 += dxi[k] * dxi[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {
        /* Interact */
        runner_iact_nonsym_stars_density(r2, dxi, hi, hj, spi, j, a, H);
      }
    }
  }
}

/**
 * @brief Compute the force on a single particle brute-force.
 */
void engine_single_density(const double dim[3], const long long int pid,
                           size_t first, const size_t N,
                           const int periodic, const struct cosmology *cosmo,
                           const struct gravity_props *grav_props) {
  const float a = 1.f;
  const float H = 0.f;

  /* Find "our" part. */
  size_t k;
  for (k = first; k < first + N && part_get_id(k) != pid; k++) {
    /* Nothing to do here */
  }

  if (k == first + N) error("Part not found.");
  size_t pind = k;

  /* Clear accumulators. */
  hydro_init_part(pind, NULL);
  adaptive_softening_init_part(pind);
  mhd_init_part(pind);

  /* Loop over all particle pairs (density). */
  for (k = first; k < first + N; k++) {
    if (k == pind) continue;

    double dx[3];
    float fdx[3];

    for (int i = 0; i < 3; i++) {
      dx[i] = part_get_x_ind(pind, i) - part_get_x_ind(k, i);
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }

    const float r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    const float hp = part_get_h(pind);
    if (r2 < hp * hp * kernel_gamma2) {
      runner_iact_nonsym_density(r2, fdx, hp, part_get_h(k), pind, k, a, H);
    }
  }

  /* Dump the result. */
  hydro_end_density(pind, cosmo);
  adaptive_softening_end_density(pind, grav_props);
  mhd_end_density(pind, cosmo);
  message("part %lli (h=%e) has wcount=%e, rho=%e.", part_get_id(pind),
          part_get_h(pind), part_get_wcount(pind), hydro_get_comoving_density(pind));
  fflush(stdout);
}

void engine_single_force(const double *dim, long long int pid,
                         size_t first, size_t N, int periodic) {
  size_t i;
  size_t k;
  double r2;
  double dx[3];
  float fdx[3];
  float a = 1.f;
  float H = 0.f;

  /* Find "our" part. */
  for (k = first; k < first + N && part_get_id(k) != pid; k++) {
    /* Nothing to do here */
  }

  if (k == first + N) error("Part not found.");
  size_t pind = k;

  /* Clear accumulators. */
  hydro_reset_acceleration(pind);
  mhd_reset_acceleration(pind);

  /* Loop over all particle pairs (force). */
  for (k = first; k < first+N; k++) {
    // for ( k = N-1 ; k >= 0 ; k-- ) {
    if (k == pind) continue;
    for (i = 0; i < 3; i++) {
      dx[i] = part_get_x_ind(pind, i) - part_get_x_ind(k, i);
      if (periodic) {
        if (dx[i] < -dim[i] / 2)
          dx[i] += dim[i];
        else if (dx[i] > dim[i] / 2)
          dx[i] -= dim[i];
      }
      fdx[i] = dx[i];
    }
    r2 = fdx[0] * fdx[0] + fdx[1] * fdx[1] + fdx[2] * fdx[2];
    const float hp = part_get_h(pind);
    const float hk = part_get_h(k);
    if (r2 < hp * hp * kernel_gamma2 || r2 < hk * hk * kernel_gamma2) {
      hydro_reset_acceleration(pind);
      mhd_reset_acceleration(pind);
      runner_iact_nonsym_force(r2, fdx, hp, hk, pind, k, a, H);
    }
  }

  /* Dump the result. */
  message("part %lli (h=%e) has a=[%.3e,%.3e,%.3e]", part_get_id(pind),
          part_get_h(pind), part_get_a_hydro_ind(pind, 0),
          part_get_a_hydro_ind(pind, 1), part_get_a_hydro_ind(pind, 2));
  fflush(stdout);
}

/**
 * Returns a random number (uniformly distributed) in [a,b)
 *
 * This function is *not* thread-safe.
 */
double random_uniform(const double a, const double b) {
  return (rand() / (((double)RAND_MAX) + 1.0)) * (b - a) + a;
}

/**
 * @brief Randomly shuffle an array of particles.
 */
void shuffle_particles(size_t first, const size_t count) {
  if (count > 1) {
    for (size_t i = 0; i < first + count - 1; i++) {
      size_t j = i + random_uniform(0., (double)(count - 1 - i));
      hydro_data_swap_particles(i, j);
    }
  }
}

/**
 * @brief Randomly shuffle an array of sparticles.
 */
void shuffle_sparticles(struct spart *sparts, const int scount) {
  if (scount > 1) {
    for (int i = 0; i < scount - 1; i++) {
      int j = i + random_uniform(0., (double)(scount - 1 - i));

      struct spart sparticle = sparts[j];

      sparts[j] = sparts[i];

      sparts[i] = sparticle;
    }
  }
}

/**
 * @brief Compares two values based on their relative difference:
 * |a - b|/|a + b|
 *
 * @param a Value a
 * @param b Value b
 * @param threshold The limit on the relative difference between the two values
 * @param absDiff (return) Absolute difference: |a - b|
 * @param absSum (return) Absolute sum: |a + b|
 * @param relDiff (return) Relative difference: |a - b|/|a + b|
 *
 * @return 1 if difference found, 0 otherwise
 */
int compare_values(double a, double b, double threshold, double *absDiff,
                   double *absSum, double *relDiff) {

  int result = 0;
  *absDiff = 0.0, *absSum = 0.0, *relDiff = 0.0;

  *absDiff = fabs(a - b);
  *absSum = fabs(a + b);
  if (*absSum > 0.f) {
    *relDiff = *absDiff / *absSum;
  }

  if (*relDiff > threshold) {
    result = 1;
  }

  return result;
}

/**
 * @brief Compares two particles' properties using the relative difference and a
 * threshold.
 *
 * @param pindi Index of particle i
 * @param pindj Index of particle j
 * @param threshold The limit on the relative difference between the two values
 *
 * @return 1 if difference found, 0 otherwise
 */
int compare_particles(size_t pindi, size_t pindj, double threshold) {

#ifdef GADGET2_SPH

  int result = 0;
  double absDiff = 0.0;
  double absSum = 0.0;
  double relDiff = 0.0;

  for (int k = 0; k < 3; k++) {
    if (compare_values(part_get_x_ind(pindi, k), part_get_x_ind(pindj, k), threshold, &absDiff, &absSum,
                       &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for x[%d] of "
          "particle %lld.",
          relDiff, threshold, k, part_get_id(pindi));
      message("a = %e, b = %e", part_get_x_ind(pindi, k), part_get_x_ind(pindj, k));
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(part_get_v_ind(pindi,k), part_get_v_ind(pindi, k),
          threshold, &absDiff, &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for v[%d] of "
          "particle %lld.",
          relDiff, threshold, k, part_get_id(pindi));
      message("a = %e, b = %e", part_get_v_ind(pindi, k), part_get_v_ind(pindj, k));
      result = 1;
    }
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(part_get_a_hydro_ind(pindi, k), part_get_a_hydro_ind(pindj,k), threshold, &absDiff,
                       &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for a_hydro[%d] "
          "of particle %lld.",
          relDiff, threshold, k, part_get_id(pindi));
      message("a = %e, b = %e", part_get_a_hydro_ind(pindi, k), part_get_a_hydro_ind(pindj, k));
      result = 1;
    }
  }
  if (compare_values(part_get_rho(pindi), part_get_rho(pindj), threshold, &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_rho(pindi), part_get_rho(pindj));
    result = 1;
  }
  if (compare_values(part_get_rho_dh(pindi), part_get_rho_dh(pindj), threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for rho_dh of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_rho_dh(pindi), part_get_rho_dh(pindj));
    result = 1;
  }
  if (compare_values(part_get_wcount(pindi), part_get_wcount(pindj), threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_wcount(pindi), part_get_wcount(pindj));
    result = 1;
  }
  if (compare_values(part_get_wcount_dh(pindi), part_get_wcount_dh(pindj), threshold,
                     &absDiff, &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for wcount_dh of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_wcount_dh(pindi), part_get_wcount_dh(pindj));
    result = 1;
  }
  if (compare_values(part_get_h_dt(pindi), part_get_h_dt(pindj), threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for h_dt of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_h_dt(pindi), part_get_h_dt(pindj));
    result = 1;
  }
  if (compare_values(part_get_v_sig(pindi), part_get_v_sig(pindj), threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for v_sig of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_v_sig(pindi), part_get_v_sig(pindj));
    result = 1;
  }
  if (compare_values(part_get_entropy_dt(pindi), part_get_entropy_dt(pindj), threshold, &absDiff, &absSum,
                     &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for entropy_dt of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_entropy_dt(pindi), part_get_entropy_dt(pindj));
    result = 1;
  }
  if (compare_values(part_get_div_v(pindi), part_get_div_v(pindj), threshold, &absDiff,
                     &absSum, &relDiff)) {
    message(
        "Relative difference (%e) larger than tolerance (%e) for div_v of "
        "particle %lld.",
        relDiff, threshold, part_get_id(pindi));
    message("a = %e, b = %e", part_get_div_v(pindi), part_get_div_v(pindj));
    result = 1;
  }
  for (int k = 0; k < 3; k++) {
    if (compare_values(part_get_rot_v_ind(pindi, k), part_get_rot_v_ind(pindj, k), threshold,
                       &absDiff, &absSum, &relDiff)) {
      message(
          "Relative difference (%e) larger than tolerance (%e) for rot_v[%d] "
          "of particle %lld.",
        relDiff, threshold, k,  part_get_id(pindi));
      message("a = %e, b = %e", part_get_rot_v_ind(pindi,k), part_get_rot_v_ind(pindj, k));
      result = 1;
    }
  }

  return result;

#else

  error("Function not supported for this flavour of SPH");
  return 0;

#endif
}

/**
 * @brief return the resident memory use of the process and its children.
 *
 * @result memory use in Kb.
 */
long get_maxrss(void) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss;
}

/**
 * @brief trim leading white space from a string.
 *
 * Returns pointer to first character.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_leading(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  while (isspace(*s)) s++;
  return s;
}

/**
 * @brief trim trailing white space from a string.
 *
 * Modifies the string by adding a NULL to the end.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_trailing(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  char *end = s + strlen(s) - 1;
  while (isspace(*end)) end--;
  *(end + 1) = '\0';
  return s;
}

/**
 * @brief trim leading and trailing white space from a string.
 *
 * Can modify the string by adding a NULL to the end.
 *
 * @param s the string.
 * @result the result.
 */
char *trim_both(char *s) {
  if (s == NULL || strlen(s) < 2) return s;
  return trim_trailing(trim_leading(s));
}

/**
 * @brief safe check directory existence and creation.
 *
 * Doesn't try to re-create an existing directory and makes sure that it can
 * have files created in it. If the directory doesn't exist and we want to
 * create it then we try to do that with an appropriate mask that should
 * honour the users umask. If the directory is something odd like a broken
 * link then that is dereferenced so should also result in an error, since
 * it doesn't exist and cannot be created. Finally we check that an existing
 * file with this name is indeed a directory.
 *
 * @param dir the directory we need to exist.
 * @param create create if not already exists, otherwise non-existence is an
 *               error.
 */
void safe_checkdir(const char *dir, int create) {

  if (access(dir, W_OK | X_OK) != 0) {

    /* A file with this name and the minimally correct permissions does not
     * exist. */
    if (create) {
      if (mkdir(dir, 0777) != 0)
        error("Failed to create directory %s (%s)", dir, strerror(errno));
    } else {
      error("directory %s does not exist or cannot be used (%s)", dir,
            strerror(errno));
    }
  } else {

    /* Exists and has the minimally correct permissions, is this a
     * directory? */
    struct stat sb;
    if (stat(dir, &sb) != 0) {
      /* Weird error, can stat all existing files. */
      error("Failed to stat directory %s (%s)", dir, strerror(errno));
    } else if ((sb.st_mode & S_IFMT) != S_IFDIR) {
      error("%s is not a directory, cannot use or create.", dir);
    }
  }

  /* Success. */
  return;
}
