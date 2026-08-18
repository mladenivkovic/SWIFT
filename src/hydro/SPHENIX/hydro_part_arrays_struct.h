#pragma once

 
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Mladen Ivkovic (mladen.ivkovic@durham.ac.uk)
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



/* Forward declare all part structs */

struct part;
struct xpart;

/**
 * Struct holding handles to all particle data arrays
 */
struct part_arrays {

  /*! @brief Particle fields for the SPH particles

The density and force substructures are used to contain variables only used
within the density and force loops over neighbours. All more permanent
variables should be declared in the main part of the part structure,
 */ 
  struct part *_part;

  /*! @brief Particle fields not needed during the SPH loops over neighbours.

This structure contains the particle fields that are not used in the
density or force loops. Quantities should be used in the kick, drift and
potentially ghost tasks only.
 */ 
  struct xpart *_xpart;

};

