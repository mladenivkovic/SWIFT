/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leideuniv.nl)
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
#ifndef SWIFT_NONE_HYDRO_DEBUG_H
#define SWIFT_NONE_HYDRO_DEBUG_H

/**
 * @file Shadowswift/hydro_debug.h
 * @brief Empty implementation.
 */
#include "error.h"
#include "hydro.h"
#include "inline.h"

#include <stddef.h>

__attribute__((always_inline)) INLINE static void hydro_debug_particle(
    const size_t pind) {
  error("Empty implementation");
}

#endif /* SWIFT_NONE_HYDRO_DEBUG_H */
